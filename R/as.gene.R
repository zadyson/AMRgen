#' @title Gene Class and AMR Parsing Functions
#' @description Functions to create a custom "gene" class and parse AMR data.
#' @param x A character vector to be converted to a "gene" class.
#' @return For `as.gene`, an object of class "gene". For [import_amrfp()], a data frame with parsed AMR data.
#' @importFrom AMR as.ab
#' @importFrom readr read_tsv
#' @importFrom dplyr filter left_join select mutate
#' @importFrom tidyr separate_longer_delim
#' @importFrom tibble add_column
#' @export
#' @examples
#' \dontrun{
#' # Create a gene object
#' gene <- as.gene(c("gene1", "gene2"))
#' print(gene)
#'
#' # Parse AMR data
#' parsed_data <- import_amrfp("path/to/input_table.tsv", "SampleID")
#' }
as.gene <- function(x) {
  # Check x is a character vector
  if (!is.character(x)) {
    stop("Input must be a character vector.")
  }
  # Create the object with class "gene"
  structure(x, class = c("gene", class(x)))
}

#' @noRd
#' @export
print.gene <- function(x, ...) {
  cat("An object of class 'gene':\n")
  print(unclass(x), ...)
}

#' Import and process AMRFinderPlus results
#'
#' This function imports and processes AMRFinderPlus results, extracting antimicrobial resistance (AMR) elements
#' and mapping them to standardised antibiotic names and drug classes. The function also converts gene symbols
#' to a harmonised format and ensures compatibility with the AMR package.
#'
#' @param input_table A character string specifying the path to the AMRFinderPlus results table (TSV format).
#' @param sample_col A character string specifying the column that identifies samples in the dataset.
#' @param amrfp_drugs A tibble containing a reference table mapping AMRFinderPlus subclasses (`AFP_Subclass`) to
#'   standardised drug agents (`drug_agent`) and drug classes (`drug_class`). Defaults to `amrfp_drugs_table`,
#'   which is provided internally.
#'
#' @return A tibble containing the processed AMR elements, with harmonised gene names, mapped drug agents, and
#'   drug classes. The output retains the original columns from the AMRFinderPlus table along with the newly mapped
#'   variables.
#'
#' @details
#' The function performs the following steps:
#' - Reads the AMRFinderPlus output table.
#' - Filters the data to only include AMR elements.
#' - Converts gene symbols to a harmonised format.
#' - Splits multiple subclass annotations into separate rows.
#' - Maps AMRFinderPlus subclasses to standardised drug agent and drug class names using `amrfp_drugs`.
#' - Converts drug agent names to the `"ab"` class from the AMR package.
#'
#' This processing ensures compatibility with downstream AMR analysis workflows.
#'
#' @importFrom AMR as.ab
#' @importFrom dplyr all_of everything filter left_join mutate select
#' @importFrom tibble add_column
#' @importFrom tidyr separate_longer_delim
#'
#' @examples
#' # Example usage:
#' \dontrun{
#' # small example E. coli AMRfinderplus data
#' ecoli_geno_raw
#'
#' # import first few rows of this data frame and parse it as AMRfp data
#' geno <- import_amrfp(ecoli_geno_raw %>% head(n = 10), "Name")
#' geno
#' }
#' @export
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

  # now split the Subclass column on the "/" to make them one per row, to make adding the ab names easier
  in_table_subclass_split <- in_table_gene %>% separate_longer_delim(Subclass, "/")

  # make two new columns - drug_class and drug_agent, where we control the vocab for the AMRFinderPlus Subclass column
  # into something that is comparable with the drugs in the AMR package
  in_table_ab <- in_table_subclass_split %>%
    left_join(., amrfp_drugs[, c("AFP_Subclass", "drug_agent", "drug_class")], by = c("Subclass" = "AFP_Subclass")) %>%
    select(all_of(sample_col), marker, drug_agent, drug_class, everything())

  # convert drug_agent into the "ab" class (will leave NAs as is)
  in_table_ab <- in_table_ab %>% mutate(drug_agent = as.ab(drug_agent))

  return(in_table_ab)
}
