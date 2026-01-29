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

# # Download EBI data
# # Zoe A. Dyson (zoe.dyson@lshtm.ac.uk)
# # Last updated 28/01/2026
#
# # List required packages
# packages <- c("arrow", "RCurl)
#
# # Install missing packages
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }
#
# # Load packages
# invisible(lapply(packages, library, character.only = TRUE))

#' Download antimicrobial genotype data from the EBI AMR Portal
#'
#' This function will retrieve genotype data from the EBI AMR Portal, https://www.ebi.ac.uk/amr. The portal uses AMRfinderplus to identify AMR-associated genotypes, but the results are processed and not all fields returned by AMRfinderplus are included. See https://www.ebi.ac.uk/amr/about/#AMR-Genotypes for more information.
#' @param user_genus String specifying a bacterial genus to download data for (default NULL, will pull all taxa)
#' @param user_release String specifying the data release to download (default NULL, will pull latest release)
#' @param user_antibiotic_name String specifying an antibiotic to download data for (default NULL, will pull all antibiotics).
#' @importFrom arrow read_parquet
#' @importFrom RCurl getBinaryURL getURL
#' @return A data frame containing EBI genotype data
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' amp_sal_ebi <- import_ebi(
#'     user_genus="Salmonella",
#'     user_release="2025-12",
#'     user_antibiotic_name="ampicillin"
#' )
#' }
import_ebi <- function(user_genus=NULL, user_release=NULL, user_antibiotic_name=NULL) {
  
  # EBI source url
  ebi_url <- "ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/"

  if (!is.null(user_release)){
    ebi_geno <- RCurl::getBinaryURL(print(gsub("ftp\\:", "https\\:", paste0(ebi_url,user_release,"/","genotype.parquet"))))
    ebi_geno <- arrow::read_parquet(ebi_geno)
    
  }else{
  
    # get list of releases
    folders <- str_split(RCurl::getURL(ebi_url, dirlistonly = TRUE), "\n")[[1]]
    # get latest
    latest_release <- folders[!folders %in% c("releases.yml", "")] %>% 
      max()
    ebi_geno <- RCurl::getBinaryURL(gsub("ftp\\:", "https\\:", paste0(ebi_url,latest_release,"/","genotype.parquet")))
    ebi_geno <- arrow::read_parquet(ebi_geno)
    
  }
  
  # filter by genus if required
  if (!is.null(user_genus)){
    ebi_geno <- ebi_geno %>% 
      filter(genus==user_genus)
  }
  
  # filter by genus if required
  if (!is.null(user_antibiotic_name)){
    ebi_geno <- ebi_geno %>% 
      filter(antibiotic_name==user_antibiotic_name)
  }
  
  # return data frame
  return(ebi_geno)
}
