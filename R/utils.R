AMRgen_env <- new.env()

utils::globalVariables(c(
  ".",
  "..count..",
  "ab",
  "Antibiotic",
  "antibiotics",
  "breakpoint_R",
  "breakpoint_S",
  "category",
  "ci.lower",
  "ci.upper",
  "AMR::clinical_breakpoints",
  "combination_id",
  "count",
  "Disk diffusio",
  "Disk diffusion (mm)",
  "disk",
  "disk_dose",
  "drug_agent",
  "drug_class",
  "ecoff_disk",
  "ecoff_MIC",
  "Element type",
  "est",
  "eucast",
  "fill_col",
  "Freq",
  "group",
  "guideline",
  "marker",
  "Measuremen",
  "Measurement sign",
  "median",
  "method",
  "MIC (mg/L)",
  "mic",
  "microorganism",
  "microorganism_code",
  "mics",
  "mo",
  "na.omit",
  "Name",
  "name",
  "NWT",
  "p",
  "perc",
  "pheno",
  "pheno_disk",
  "pheno_MIC",
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
  "Testin",
  "Testing standard",
  "type",
  "u",
  "user",
  "value",
  "x"
))


#' @importFrom readr read_tsv
process_input <- function(input) {
  # Check if the input is a file path (string)
  if (is.character(input) && file.exists(input)) {
    message("Input is a file path, reading the file.")
    data <- read_tsv(input)
  }

  # Check if the input is already a dataframe
  else if (is.data.frame(input)) {
    message("Input is already a dataframe.")
    data <- input
  }

  # If the input is neither a file nor a dataframe, stop with an error
  else {
    stop("Input must be either a valid file path or a dataframe.")
  }

  # Return the dataframe
  return(data)
}
