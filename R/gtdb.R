#' Get Microorganism from GTDB Species Name
#'
#' Parse a character vector containing species names from GTDB output to get a valid microorganism code.
#' @param species Name of species, coerced with [AMR::as.mo()].
#' @importFrom AMR as.mo mo_cleaning_regex
#' @return A character vector of valid microorganism codes as determined by [AMR::as.mo()].
#' @export
#' @examples
#' gtdb.mo("Escherichia_A coli_BC")
gtdb.mo <- function(species) {
  as.mo(species, cleaning_regex = paste0(mo_cleaning_regex(), "|_[A-Z]+($| )"), keep_synonyms = TRUE)
}

#' Import GTDB Output
#'
#' Import GTDB output (from file or data frame) and parse the species name into a microorganism recognised by AMR package. AMR mo code and species name will be added to the data frame.
#' @param file File path to GTDB output file (tab-separated). If not given, `tbl` must be provided.
#' @param tbl Data frame containing GTDB output. If not given, `file` must be provided.
#' @param species_column Name of column containing the species call (default "Species").
#' @importFrom AMR mo_name
#' @importFrom dplyr mutate
#' @importFrom readr read_tsv
#' @return A data frame containing GTDB output with AMR-parsed microorganism code (`gtdb.mo`) and species name (`gtdb.species`) appended.
#' @export
#' @examples
#' \dontrun{
#' import_gtdb(tbl = data.frame(Species = c(
#'   "Pseudomonas_E piscis",
#'   "Haemophilus_D parainfluenzae_A",
#'   "Acinetobacter calcoaceticus_C"
#' )))
#' }
import_gtdb <- function(file = NULL, tbl = NULL, species_column = "Species") {
  if (!is.null(file)) {
    tbl <- read_tsv(file)
  }
  tbl %>%
    mutate(gtdb.mo = gtdb.mo(tbl[[species_column]]),
           gtdb.species = mo_name(gtdb.mo, keep_synonyms = TRUE))
}
