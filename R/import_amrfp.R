#' Import and Process AMRFinderPlus Results
#'
#' This function imports and processes AMRFinderPlus results, extracting antimicrobial resistance (AMR) elements and mapping them to standardised antibiotic names and drug classes. The function also converts gene symbols to a harmonised format and ensures compatibility with the AMR package.
#' @param input_table A character string specifying the path to the AMRFinderPlus results table (TSV format).
#' @param sample_col A character string specifying the column that identifies samples in the dataset.
#' @param amrfp_drugs A tibble containing a reference table mapping AMRFinderPlus subclasses (`AFP_Subclass`) to standardised drug agents (`drug_agent`) and drug classes (`drug_class`). Defaults to `amrfp_drugs_table`, which is provided internally.
#' @importFrom AMR as.ab
#' @importFrom dplyr all_of everything filter left_join mutate select
#' @importFrom tibble add_column
#' @importFrom tidyr separate_longer_delim, separate
#' @return A tibble containing the processed AMR elements, with harmonised gene names, mapped drug agents, and drug classes. The output retains the original columns from the AMRFinderPlus table along with the newly mapped variables.
#' @details
#' The function performs the following steps:
#' - Reads the AMRFinderPlus output table.
#' - Filters the data to only include AMR elements.
#' - Converts gene symbols to a harmonised format.
#' - Splits multiple subclass annotations into separate rows.
#' - Maps AMRFinderPlus subclasses to standardised drug agent and drug class names using `amrfp_drugs`.
#' - Converts drug agent names to the `"ab"` class from the AMR package.
#' This processing ensures compatibility with downstream AMR analysis workflows.
#' @export
#' @examples
#' \dontrun{
#' # small example E. coli AMRfinderplus data
#' ecoli_geno_raw
#'
#' # import first few rows of this data frame and parse it as AMRfp data
#' geno <- import_amrfp(ecoli_geno_raw %>% head(n = 10), "Name")
#' geno
#' }
import_amrfp <- function(input_table, sample_col, amrfp_drugs = amrfp_drugs_table) {
  in_table <- process_input(input_table)

  # filter to only include AMR elements
  if ("Element type" %in% colnames(in_table)) {
    in_table <- in_table %>% filter(`Element type` == "AMR")
  } else {
    print("No `Element type` column found, assuming all rows report AMR markers")
  }

  # create the gene class for the gene column
  if ("Gene symbol" %in% colnames(in_table)) {
    gene_col <- as.gene(in_table$`Gene symbol`)
  } else {
    stop("Expected column `Gene symbol` not found")
  }

  # Find the position of the "Gene symbol" column and insert gene after it (could replace it really)
  position <- which(names(in_table) == "Gene symbol")
  in_table_gene <- add_column(in_table, marker = gene_col, .after = position)
  
  # separate into gene and mutation, with mutation in HGVS format
  if ("Element subtype" %in% colnames(in_table)) {
    in_table_mutation <- in_table_gene %>% separate(marker, into = c("gene", "mutation"), sep = "_", remove=F, fill="right") %>%
      mutate(gene=if_else(`Element subtype`=="POINT", gene, as.character(marker))) %>%
      mutate(mutation=if_else(`Element subtype`=="POINT", convert_mutation(mutation), NA))
  } else {
    print("No `Element subtype` column found, cannot parse mutations")
    in_table_mutation <- in_table_gene %>% mutate(gene=NA, mutation=NA)
  }
  
  # create AMRrules style label with node:mutation
  if ("Hierarchy node" %in% colnames(in_table_mutation)) {
    in_table_label <- in_table_mutation %>% 
      mutate(marker.label=case_when(`Element subtype`=="POINT" & !is.na(`Hierarchy node`) ~ paste0(`Hierarchy node`,mutation,collapse=":"),
                                    `Element subtype`=="POINT" & is.na(`Hierarchy node`) ~ paste0(gene,mutation,collapse=":"),
                                    `Element subtype`!="POINT" & !is.na(`Hierarchy node`) ~ `Hierarchy node`,
                                    `Element subtype`!="POINT" & is.na(`Hierarchy node`) ~ gene,
                                    TRUE ~ NA))
    
    in_table_label <- in_table_mutation %>%
      mutate(node=if_else(is.na(`Hierarchy node`), gene, `Hierarchy node`)) %>%
      mutate(marker.label=if_else(`Element subtype`=="POINT", 
                                  paste0(node,":",mutation),
                                  node))
  }
  else {mutate(marker.label=NA, node=NA)}

  # now split the Subclass column on the "/" to make them one per row, to make adding the ab names easier
  in_table_subclass_split <- in_table_label %>% separate_longer_delim(Subclass, "/")

  # make two new columns - drug_class and drug_agent, where we control the vocab for the AMRFinderPlus Subclass column
  # into something that is comparable with the drugs in the AMR package
  in_table_ab <- in_table_subclass_split %>%
    left_join(., amrfp_drugs[, c("AFP_Subclass", "drug_agent", "drug_class")], by = c("Subclass" = "AFP_Subclass")) %>%
    relocate(any_of(c("marker","gene","mutation","node","marker.label", "drug_agent", "drug_class")), .before=`Gene symbol`)

  # convert drug_agent into the "ab" class (will leave NAs as is)
  in_table_ab <- in_table_ab %>% mutate(drug_agent = as.ab(drug_agent))

  return(in_table_ab)
}

#' @export
convert_mutation <- function(mut) {
  
  # to consider
  # AMRfp encodes protein mutations as adeS_T153M or adeS_Q163STOP or amrR_ASLD153del
  #    promoter mutations have -X position ampC_C-42G
  #    some protein substitutions might have multiple AA on either side of the coordinate
  # weird one: cirA_Y253CfsTer5
  
  aa_mapping <- c(
    A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", E = "Glu",
    Q = "Gln", G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys",
    M = "Met", F = "Phe", P = "Pro", S = "Ser", T = "Thr", W = "Trp",
    Y = "Tyr", V = "Val"
  )
  aa1 <- substr(mut, 1, 1)
  position <- substr(mut, 2, nchar(mut) - 1)
  aa2 <- substr(mut, nchar(mut), nchar(mut))
  return(paste0("p.", aa_mapping[aa1], position, aa_mapping[aa2]))
}
