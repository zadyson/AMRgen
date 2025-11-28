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

AMRgen_env <- new.env()

utils::globalVariables(c(
  ".",
  "..count..",
  "2.5 %",
  "97.5 %",
  "ab",
  "AMR::clinical_breakpoints",
  "Antibiotic",
  "antibiotics",
  "binary_comb",
  "breakpoint_R",
  "breakpoint_S",
  "category",
  "ci.lower",
  "ci.lower.est",
  "ci.lower.ppv",
  "ci.upper",
  "ci.upper.est",
  "ci.upper.ppv",
  "coalesce",
  "combination_id",
  "count",
  "Disk diffusio",
  "Disk diffusion (mm)",
  "disk",
  "disk_dose",
  "drug_agent",
  "drug_class",
  "ecoff_disk",
  "ecoff_mic",
  "ecoff_MIC",
  "Element type",
  "Element subtype",
  "est",
  "Estimate",
  "eucast",
  "fill_col",
  "Freq",
  "Gene symbol",
  "gene",
  "get(assay)",
  "get(marker_col)",
  "group",
  "guideline",
  "Hierarchy node",
  "id",
  "marker",
  "marker.label",
  "marker_count",
  "marker_list",
  "Measuremen",
  "Measurement sign",
  "median",
  "Method",
  "method",
  "MIC (mg/L)",
  "mic",
  "microorganism",
  "microorganism_code",
  "mics",
  "mo",
  "mutation",
  "na.omit",
  "Name",
  "name",
  "node",
  "NWT",
  "p",
  "perc",
  "pheno",
  "pheno_disk",
  "pheno_mic",
  "pheno_MIC",
  "point_size",
  "ppv",
  "Pr(>|z|)",
  "pval",
  "py",
  "py_run_string",
  "R",
  "r_to_py",
  "Scientifi",
  "Scientific name",
  "se",
  "setNames",
  "sig_binary",
  "solo",
  "Source",
  "spp_pheno",
  "Subclass",
  "subtype",
  "symbol",
  "Testin",
  "Testing standard",
  "type",
  "type",
  "u",
  "user",
  "value",
  "variation type",
  "x"
))


#' @importFrom readr read_tsv
process_input <- function(input) {
  if (is.character(input) && file.exists(input)) {
    # Get the full file extension (e.g., "csv.gz" or "tsv")
    full_extension <- tolower(tools::file_ext(input))
    if (full_extension %in% paste(rep(c("tsv", "txt"), 5), c(rep("", 2), rep(".gz", 2), rep(".bz2", 2), rep(".xz", 2), rep(".zip", 2)), sep = "")) {
      data <- readr::read_tsv(input)
    } else if (full_extension %in% paste(rep("csv", 5), c("", ".gz", ".bz2", ".xz", ".zip"), sep = "")) {
      data <- readr::read_csv(input)
    }
  } else if (is.data.frame(input)) {
    # Check if the input is already a dataframe
    data <- input
  } else {
    # If the input is neither a file nor a dataframe, stop with an error
    stop("Input must be either a valid file path or a dataframe.")
  }
  # Return the dataframe
  return(data)
}

#' @importFrom cli ansi_has_hyperlink_support
font_url <- function(url, txt = url) {
  if (tryCatch(isTRUE(ansi_has_hyperlink_support()), error = function(e) FALSE)) {
    paste0("\033]8;;", url, "\a", txt, "\033]8;;\a")
  } else {
    url
  }
}

#' @importFrom cli ansi_has_hyperlink_support
font_italic <- function(..., collapse = " ") {
  txt <- paste0(c(...), collapse = collapse)
  if (tryCatch(isTRUE(ansi_has_hyperlink_support()), error = function(e) FALSE)) {
    if (is.null(collapse)) {
      paste0("\033[3m", txt, "\033[23m", collapse = NULL)
    } else {
      paste0("\033[3m", txt, "\033[23m", collapse = "")
    }
  } else {
    txt
  }
}
