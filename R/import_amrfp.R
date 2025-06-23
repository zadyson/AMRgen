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

  # detect variation type
  if ("Method" %in% colnames(in_table)) {
    in_table_mutation <- in_table %>% 
      mutate(marker=`Gene symbol`) %>%
      mutate(`variation type`=case_when(Method=="INTERNAL_STOP" ~ "Inactivating mutation detected",
                                      grepl("PARTIAL", Method) ~ "Inactivating mutation detected",
                                      Method=="POINTN" ~ "Nucleotide variant detected",
                                      Method %in% c("POINTX", "POINTP") ~ "Protein variant detected",
                                      Method %in% c("ALLELEP", "ALLELEX", "BLASTP", "BLASTX", "EXACTP", "EXACTX") ~ "Gene presence detected",
                                      TRUE ~ NA)) %>% 
      separate(marker, into = c("gene", "mutation"), sep = "_", remove=F, fill="right") %>%
      mutate(gene=if_else(startsWith(Method,"POINT"), gene, marker)) %>%
      mutate(mutation=if_else(startsWith(Method,"POINT"),convert_mutation(marker, Method), mutation))
  } else {
    print("Need Method columns to assign to parse mutations and assign variation type")
    in_table_mutation <- marker %>% mutate(`variation type`=NA, gene=NA, mutation=NA)
  }
  
  # create AMRrules style label with node:mutation
  if ("Hierarchy node" %in% colnames(in_table_mutation) & "Element subtype" %in% colnames(in_table_mutation)) {
    in_table_label <- in_table_mutation %>%
      mutate(node=if_else(is.na(`Hierarchy node`), gene, `Hierarchy node`)) %>%
      mutate(marker.label=if_else(`Element subtype`=="POINT", 
                                  if_else(!is.na(mutation), paste0(node,":",mutation), gsub("_",":",marker)),
                                  node)) %>%
      mutate(marker.label=if_else(`variation type`=="Inactivating mutation detected", 
                                  paste0(node,":-"),
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


#' Convert mutation string based on method
#'
#' This function takes a mutation string (e.g., "gene_REF123ALT") and a
#' mutation method, then extracts and converts parts of the mutation string.
#' Specifically designed for use within `dplyr::mutate()`.
#'
#' @param gene_symbol_col A character vector representing the 'Gene symbol'
#'                        column (the mutation string).
#' @param method_col A character vector representing the 'Method' column.
#' @return A character vector containing the formatted mutation strings
#'         (e.g., "Ala123Trp") or NA if not applicable/match.
#' @examples
#' data_for_mutation <- tibble(
#'   `Gene symbol` = c("gyrA_S83I","ompK36_D135DGD","cirA_W427STOP", "cirA_Y253CfsTer5", "blaSHV_C-112A"),
#'   Method = c("POINTX", "POINTP", "POINTX", "POINTP", "POINTN")
#' )
#'
#' data_for_mutation %>%
#'   mutate(mutation = convert_mutation(`Gene symbol`, Method))
convert_mutation <- function(gene_symbol_col, method_col) {
  
  # Ensure inputs are treated as vectors. `mutate` will pass them as such.
  # The regex pattern for splitting the mutation string
  regex_pattern <- "^([^_]+)_([A-Za-z]+)(-?)(\\d+)([A-Za-z]+)$"
  
  # Apply str_match to the entire gene_symbol_col vector.
  # str_match is already vectorized, so it handles all rows at once.
  result_matrix <- str_match(gene_symbol_col, regex_pattern)
  
  # Convert the result matrix to a tibble immediately
  extracted_data <- tibble(
    ref = result_matrix[,3],
    minus = result_matrix[,4],
    position = result_matrix[,5],
    alt = result_matrix[,6]
  )
  
  # Perform the conditional conversion based on Method and apply convert_aa_code
  new_ref <- ifelse(method_col %in% c("POINTX", "POINTP"), convert_aa_code(extracted_data$ref), extracted_data$ref)
  new_alt <- ifelse(method_col %in% c("POINTX", "POINTP"), convert_aa_code(extracted_data$alt), extracted_data$alt)
  
  # Construct the final mutation string
  # `paste0` is vectorized. We need to handle NA carefully if new_ref, position, or new_alt are NA.
  new_mutation_string <- ifelse(
    method_col %in% c("POINTX", "POINTP") ,
    ifelse(!is.na(new_ref) & !is.na(new_alt) & !is.na(extracted_data$position),
           paste0(new_ref, extracted_data$minus, extracted_data$position, new_alt),
           NA_character_ # Return NA if any of the components are NA
    ),
    paste0(extracted_data$minus, extracted_data$position, new_ref, ">", new_alt)
  )
  
  return(new_mutation_string)
}