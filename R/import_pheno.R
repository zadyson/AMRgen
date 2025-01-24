#' Import and process AST data from an NCBI file
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes
#' the data, and optionally interprets the results based on MIC or disk diffusion data.
#' It assumes that the input file is a delimited text file (e.g., CSV, TSV) and
#' parses relevant columns (antibiotic names, species names, MIC or disk data) into 
#' suitable classes using the AMR package. It optionally can use the AMR package to  
#' determine susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines. If expected
#' columns are not found warnings will be given, and interpretation may not be possible.
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
#'   values, and either EUCAST or CLSI testing standard (as indicated in the `Testing standard`
#'   column of the input file, if blank EUCAST will be used by default). If `FALSE`, no 
#'   interpretation is performed.
#'
#' @return A data frame with the processed AST data, including additional columns:
#'   \itemize{
#'     \item `id`: The biological sample identifier (renamed from `#BioSample` or specified column).
#'     \item `spp_pheno`: The species phenotype, formatted using the `as.mo` function.
#'     \item `drug_agent`: The antibiotic used in the test, formatted using the `as.ab` function.
#'     \item `mic`: The minimum inhibitory concentration (MIC) value, formatted using the `as.mic` function.
#'     \item `disk`: The disk diffusion measurement (in mm), formatted using the `as.disk` function.
#'     \item `guideline`: The guideline used for interpretation (either EUCAST or CLSI).
#'     \item `pheno`: The interpreted phenotype (SIR) based on the MIC or disk diffusion data.
#'   }
#'
#'
#' @examples
#' # Example usage
#' ast_data <- import_ncbi_ast("path/to/ast_data.tsv")
#' head(ast_data)
#'
#' @export
import_ncbi_ast <- function(input, sample_col="#BioSample", interpret=F) {
  
  ast <- process_input(input)
  
  # find id column
  if (!is.null(sample_col)) {
    if (sample_col %in% colnames(ast)) {ast <- ast %>% rename(id=any_of(sample_col))}
    else {stop(paste("Invalid column name:", sample_col))}
  }
  else {stop("Please specify the column containing sample identifiers, via parameter 'sample_col'")}
  
  # parse guideline column
  if ("Testing standard" %in% colnames(ast)) {
    ast <- ast %>% mutate(guideline=if_else(is.na(`Testing standard`), "EUCAST", `Testing standard`), .after=id) # interpret as EUCAST if not specified as CLSI
  }
  else {print("Warning: Expected column 'Testing standard' not found in input")}
 
  # parse disk column
  if ("Disk diffusion (mm)" %in% colnames(ast)) {
    ast <- ast %>% mutate(disk=as.disk(`Disk diffusion (mm)`), .after=id)
  } 
  else {print("Warning: Expected column 'Disk diffusion (mm)' not found in input")}

  # parse mic column
  if ("MIC (mg/L)" %in% colnames(ast)) {
    if ("Measurement sign" %in% colnames(ast)) {
      ast <- ast %>% 
        mutate(mic=paste0(`Measurement sign`, `MIC (mg/L)`), .after=id) %>%
        mutate(mic=gsub("==","",mic))
    }
    else {
      ast <- ast %>% mutate(mic=`MIC (mg/L)`, .after=id)
      print("Warning: Expected column 'Measurement sign' not found in input, be careful interpreting MIC")
    }
    ast <- ast %>% mutate(mic=as.mic(mic))
  } 
  else {print("Warning: Expected column 'MIC (mg/L)' not found in input")}

  # parse antibiotic column
  if ("Antibiotic" %in% colnames(ast)) {
    ast <- ast %>% mutate(drug_agent=as.ab(Antibiotic), .after=id)
  }
  else {stop("Expected column 'Antibiotic' not found in input.")}

  # parse species column
  if ("Scientific name" %in% colnames(ast)) {
    ast <- ast %>% mutate(spp_pheno=as.mo(`Scientific name`), .after=id)
  }
  else {print("Warning: Expected column 'Scientific name' not found in input")}
 
  if (interpret) {
    if (all(c("spp_pheno", "drug_agent", "guideline") %in% colnames(ast))) {
      if (!("mic" %in% colnames(ast)) & !("disk" %in% colnames(ast))) {
        print("Could not interpret phenotypes, no mic or disk columns found")
      }
      else {
        if (!("mic" %in% colnames(ast))) {ast <- ast %>% mutate(mic=NA)}
        if (!("disk" %in% colnames(ast))) {ast <- ast %>% mutate(disk=NA)}
        ast <- ast %>% rowwise() %>% 
          mutate(pheno_MIC=if_else(!is.na(mic), 
                               as.sir(mic, mo=spp_pheno, ab=drug_agent, guideline=guideline),
                               NA)) %>%
          mutate(pheno_disk=if_else(!is.na(disk), 
                                   as.sir(disk, mo=spp_pheno, ab=drug_agent, guideline=guideline),
                                   NA)) %>%
          unite(pheno, pheno_MIC:pheno_disk, na.rm=T)
      }
    }
  }
  
  return (ast)
  
}
