# ===================================================================== #
#  Licensed as GPL-v3.0.                                                #
#                                                                       #
#  Developed as part of the AMRverse (https://github.com/AMRverse):     #
#  https://github.com/AMRverse/AMRgen                                   #
#                                                                       #
#  We created this package for both routine data analysis and academic  #
#  research and it was publicly released in the hope that it will be    #
#  useful, but it comes WITHOUT ANY WARRANTY OR LIABILITY.              #
#                                                                       #
#  This R package is free software; you can freely use and distribute   #
#  it for both personal and commercial purposes under the terms of the  #
#  GNU General Public License version 3.0 (GNU GPL-3), as published by  #
#  the Free Software Foundation.                                        #
# ===================================================================== #

#' Gene Class and AMR Parsing Functions
#'
#' Functions to create a custom "gene" class and parse AMR data.
#' @param x A character vector to be converted to a "gene" class.
#' @importFrom AMR as.ab
#' @importFrom readr read_tsv
#' @importFrom dplyr filter left_join select mutate
#' @importFrom tidyr separate_longer_delim
#' @importFrom tibble add_column
#' @return For `as.gene`, an object of class `"gene"`.
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
