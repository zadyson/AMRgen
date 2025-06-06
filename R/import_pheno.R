#' Import and process AST data from an NCBI file
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes
#' the data, and optionally interprets the results based on MIC or disk diffusion data.
#' It assumes that the input file is a tab-delimited text file (e.g., TSV) and
#' parses relevant columns (antibiotic names, species names, MIC or disk data) into
#' suitable classes using the AMR package. It optionally can use the AMR package to
#' determine susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines (human
#' breakpoints and/or ECOFF). If expected columns are not found warnings will be given,
#' and interpretation may not be possible.
#'
#' @param input A string representing a dataframe, or a path to a tab-delimited file,
#'    containing the AST data in NCBI antibiogram format. These files can be downloaded
#'    fromNCBI AST browser, e.g. https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa
#'
#' @param sample_col A string indicating the name of the column with sample identifiers.
#'    If `NULL`, assume this is '#BioSample'.
#'
#' @param interpret A logical value (default is FALSE). If `TRUE`, the function will interpret
#'   the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion
#'   values, against human breakpoints from either EUCAST or CLSI testing standard (as indicated
#'   in the `Testing standard` column of the input file, if blank the value of the default_guideline
#'   parameter will be used by default). If `FALSE`, no interpretation is performed.
#'
#' @param ecoff A logical value (default is FALSE). If `TRUE`, the function will interpret
#'   the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion
#'   values, against epidemiological cut-off (ECOFF) values. These will be reported in
#'   a new column 'ecoff', coded as 'NWT' (nonwildtype) or 'WT' (wildtype). If `FALSE`,
#'   no ECOFF interpretation is performed.
#'
#' @param default_guideline A string (default is "EUCAST"). Default guideline to use for
#'   interpretation via as.sir. Allowed values are 'EUCAST' or 'CLSI'. If the input file
#'   contains a column `Testing standard`, or if interpret or ecoff are set to `TRUE`, a
#'   new column `guideline` will be created to use in the interpretation step. Values are
#'   populated from those in `Testing standard`, however rows with missing/NA values or
#'   non-allowed values will be coerced to the value specified by 'default_guideline'.
#'   If there is no `Testing standard` column, all rows will be interpreted using the
#'   default_guideline.
#'
#' @return A data frame with the processed AST data, including additional columns:
#'   \itemize{
#'     \item `id`: The biological sample identifier (renamed from `#BioSample` or specified column).
#'     \item `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#'     \item `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#'     \item `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#'     \item `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#'     \item `guideline`: The guideline used for interpretation (either EUCAST or CLSI; taken from input column otherwise forced to parameter default_guideline).
#'     \item `pheno`: The phenotype interpreted against the specified breakpoint standard (as S/I/R), based on the MIC or disk diffusion data.
#'     \item `ecoff`: The wildtype/nonwildtype status interpreted against the ECOFF (as WT/NWT), based on the MIC or disk diffusion data.
#'   }
#'
#'
#' @examples
#' # Example usage
#' \dontrun{
#' # small example E. coli AST data from NCBI
#' ecoli_ast_raw
#'
#' # import without re-interpreting resistance
#' pheno <- import_ncbi_ast(ecoli_ast_raw)
#' head(pheno)
#'
#' # import and re-interpret resistance (S/I/R) and ECOFF (WT/NWT) using AMR package
#' pheno <- import_ncbi_ast(ecoli_ast_raw, interpret = T, ecoff = T)
#' head(pheno)
#' }
#' @importFrom AMR as.ab as.disk as.mic as.mo as.sir
#' @importFrom dplyr any_of case_when if_else mutate relocate rename rowwise
#' @importFrom tidyr unite
#' @export
import_ncbi_ast <- function(input, sample_col = "#BioSample", interpret = F, ecoff = F, default_guideline = "EUCAST") {
  ast <- process_input(input)

  # find id column
  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(ast)) {
      ast <- ast %>% rename(id = any_of(sample_col))
    } else {
      stop(paste("Invalid column name:", sample_col))
    }
  } else {
    stop("Please specify the column containing sample identifiers, via parameter 'sample_col'")
  }

  # parse guideline column
  if ("Testing standard" %in% colnames(ast)) {
    ast <- ast %>% mutate(guideline = if_else(`Testing standard` %in% c("CLSI", "EUCAST"),
      `Testing standard`, default_guideline
    ), .after = id)
  } else {
    print("Warning: Expected column 'Testing standard' not found in input")
    if (interpret | ecoff) {
      print(paste("Specified default standard", default_guideline, "will be used for interpretation"))
      ast <- ast %>% mutate(guideline = default_guideline, .after = id)
    }
  }

  # parse disk column
  if ("Disk diffusion (mm)" %in% colnames(ast)) {
    ast <- ast %>% mutate(disk = as.disk(`Disk diffusion (mm)`), .after = id)
  } else {
    print("Warning: Expected column 'Disk diffusion (mm)' not found in input")
  }

  # parse mic column
  if ("MIC (mg/L)" %in% colnames(ast)) {
    if ("Measurement sign" %in% colnames(ast)) {
      ast <- ast %>%
        mutate(mic = paste0(`Measurement sign`, `MIC (mg/L)`), .after = id) %>%
        mutate(mic = gsub("==", "", mic))
    } else {
      ast <- ast %>% mutate(mic = `MIC (mg/L)`, .after = id)
      print("Warning: Expected column 'Measurement sign' not found in input, be careful interpreting MIC")
    }
    ast <- ast %>% mutate(mic = as.mic(mic))
  } else {
    print("Warning: Expected column 'MIC (mg/L)' not found in input")
  }

  # parse antibiotic column
  if ("Antibiotic" %in% colnames(ast)) {
    ast <- ast %>% mutate(drug_agent = as.ab(Antibiotic), .after = id)
  } else {
    stop("Expected column 'Antibiotic' not found in input.")
  }

  # parse species column
  if ("Scientific name" %in% colnames(ast)) {
    ast <- ast %>% mutate(spp_pheno = as.mo(`Scientific name`), .after = id)
  } else {
    print("Warning: Expected column 'Scientific name' not found in input")
  }

  if (interpret | ecoff) {
    if (all(c("spp_pheno", "drug_agent") %in% colnames(ast))) {
      if (!("mic" %in% colnames(ast)) & !("disk" %in% colnames(ast))) {
        print("Could not interpret phenotypes, no mic or disk columns found")
      } else {
        if (!("mic" %in% colnames(ast))) {
          ast <- ast %>% mutate(mic = NA)
        }
        if (!("disk" %in% colnames(ast))) {
          ast <- ast %>% mutate(disk = NA)
        }
        if (interpret) {
          ast <- ast %>%
            rowwise() %>%
            mutate(pheno_MIC = if_else(!is.na(mic),
              as.sir(mic, mo = spp_pheno, ab = drug_agent, guideline = guideline),
              NA
            )) %>%
            mutate(pheno_disk = if_else(!is.na(disk),
              as.sir(disk, mo = spp_pheno, ab = drug_agent, guideline = guideline),
              NA
            )) %>%
            unite(pheno, pheno_MIC:pheno_disk, na.rm = T) %>%
            mutate(pheno = as.sir(pheno)) %>%
            relocate(pheno, .after = drug_agent)
        }
        if (ecoff) {
          ast <- ast %>%
            rowwise() %>%
            mutate(ecoff_MIC = if_else(!is.na(mic),
              as.sir(mic, mo = spp_pheno, ab = drug_agent, guideline = guideline, breakpoint_type = "ECOFF"),
              NA
            )) %>%
            mutate(ecoff_disk = if_else(!is.na(disk),
              as.sir(disk, mo = spp_pheno, ab = drug_agent, guideline = guideline, breakpoint_type = "ECOFF"),
              NA
            )) %>%
            unite(ecoff, ecoff_MIC:ecoff_disk, na.rm = T) %>%
            mutate(ecoff = case_when(ecoff == "R" ~ "NWT", ecoff == "S" ~ "WT", TRUE ~ NA)) %>%
            relocate(ecoff, .after = drug_agent)
        }
      }
    }
  }

  return(ast)
}
