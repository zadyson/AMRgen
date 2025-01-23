#' Import and process AST data from an NCBI file
#'
#' This function imports an antibiotic susceptibility testing (AST) dataset, processes
#' the data, and optionally interprets the results based on MIC or disk diffusion data.
#' It assumes that the input file is a delimited text file (e.g., CSV, TSV) and
#' parses relevant columns (antibiotic names, species names, MIC or disk data) into 
#' suitable classes using the AMR package. It also uses the AMR package to determine 
#' susceptibility phenotype (SIR) based on EUCAST or CLSI guidelines
#'
#' @param file A string representing the path to the delimited file containing the AST data 
#'    in NCBI antibiogram format. These files can be downloaded from the NCBI AST browser, 
#'    e.g. https://www.ncbi.nlm.nih.gov/pathogens/ast#Pseudomonas%20aeruginosa
#' @param interpret A logical value (default is TRUE). If `TRUE`, the function will interpret 
#'   the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion 
#'   values, and either EUCAST or CLSI testing standard (as indicated in the `Testing standard`
#'   column of the input file, if blank EUCAST will be used by default). If `FALSE`, no 
#'   interpretation is performed.
#'
#' @return A data frame with the processed AST data, including additional columns:
#'   \itemize{
#'     \item `biosample`: The biological sample identifier (renamed from `#BioSample`).
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
#' ast_data <- import_ncbi_ast("path/to/ast_data.tsv", interpret = TRUE)
#' head(ast_data)
#'
#' @export
import_ncbi_ast <- function(file, interpret=T) {
  
  ast <- read_delim(file) %>% rename(biosample=`#BioSample`) %>% 
    mutate(spp_pheno=as.mo(`Scientific name`), .after=biosample) %>%
    mutate(drug_agent=as.ab(Antibiotic), .after=spp_pheno) %>%
    mutate(mic=as.mic(`MIC (mg/L)`), .after=drug_agent) %>%
    mutate(disk=as.disk(`Disk diffusion (mm)`), .after=mic) %>%
    mutate(guideline=if_else(is.na(`Testing standard`), "EUCAST", `Testing standard`), .after=disk) # interpret as EUCAST if not specified as CLSI

  if (interpret) {
    ast <- ast %>% rowwise() %>% 
      mutate(pheno=if_else(!is.na(mic), 
                           as.sir(mic, mo=spp_pheno, ab=drug_agent, guideline=guideline), 
                           as.sir(disk, mo=spp_pheno, ab=drug_agent, guideline=guideline)),
             .after=guideline)

  }
  
  return (ast)
  
}
