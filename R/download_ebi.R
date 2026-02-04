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
#' Download antimicrobial genotype or phenotype data from the EBI AMR Portal
#'
#' This function will retrieve genotype or phenotype data from the EBI AMR Portal, via FTP. The portal uses AMRfinderplus to identify AMR-associated genotypes, but the results are processed and not all fields returned by AMRfinderplus are included. See https://www.ebi.ac.uk/amr/about/#AMR-Genotypes for more information, and https://github.com/ncbi/amr/wiki/class-subclass for valid class and subclass terms.
#' Optionally, the function can also reformat the phenotype data for easy use with AMRgen functions (using `import_ebi_ast_ftp()`) and re-interpret assay measures using the latest breakpoints/ECOFF.
#' @param data String specifying the type of data to download, either "phenotype" or "genotype" (default "phenotype").
#' @param genus (Optional) String specifying a bacterial genus to filter on (default NULL, will pull all taxa).
#' @param species (Optional) String specifying a bacterial species to filter on (default NULL, will pull all taxa). Not used if genus is specified.
#' @param antibiotic (Optional) String specifying an antibiotic to to filter on (default NULL). Not used if class or subclass is specified.
#' @param geno_subclass (Optional) String specifying an antibiotic subclass to filter genotype data on (default NULL). Filter is based on string match, not identity, so e.g. subclass="TRIMETHOPRIM" will return all rows where the string "TRIMETHOPRIM" is included in the subclass field. Only used if `data`="genotype". Check NCBI Subclass for valid terms.
#' @param geno_class (Optional) String specifying an antibiotic subclass to filter on (default NULL). Filter is based on string match, not identity, so e.g. class="TRIMETHOPRIM" will return all rows where the string "TRIMETHOPRIM" is included in the class field. Only used if `data`="genotype" and subclass is not specified. Check NCBI Class for valid terms.
#' @param remove_dup (Optional) Logical specifying whether to clean up genotype data by removing duplicates for the same hit (default FALSE). Where a detected gene is associated with a subclass value that is actually a list of multiple drugs/classes, e.g. subclass="GENTAMICIN/TOBRAMYCIN", the EBI data table will have duplicate rows for the same gene hit, but with different values for the `antibiotic_name` and associated `antibiotic_ontology` and `antibiotic_ontology_link` annotation fields (e.g. one row each for gentamicin and tobramycin). To remove these duplicate rows (and the drug-specific annotation fields) and return only one row per hit (i.e. restoring AMRfinderplus output format), set this to TRUE.
#' @param reformat (Optional) Logical specifying whether to reformat the downloaded phenotype data for easy use with downstream AMRgen function, using `import_ebi_ast_ftp()`. Default `FALSE`. This does things like format the antibiotic, measurement, and phenotype columns to AMR package classes. If set to `TRUE` you can also turn on re-interpreting MIC/disk data using latest EUCAST/CLSI breakpoints. No columns are removed from the downloaded data frame, but key fields are renamed, see documentation for `format_ast()`.
#' @param interpret_eucast (Optional) Logical specifying whether to re-interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against EUCAST human breakpoints. These will be reported in a new column `pheno_eucast`, of class 'sir'. Only used when downloading phenotype data, with reformat set to `TRUE`.
#' @param interpret_clsi (Optional) Logical specifying whether to re-interpret the susceptibility phenotype (SIR) for each row based on the MIC or disk diffusion values, against CLSI human breakpoints. These will be reported in a new column `pheno_clsi`, of class 'sir'. Only used when downloading phenotype data, with reformat set to `TRUE`.
#' @param interpret_ecoff (Optional) Logical specifying whether to re-interpret the wildtype vs nonwildtype status for each row based on the MIC or disk diffusion values, against epidemiological cut-off (ECOFF) values. These will be reported in a new column `ecoff`, of class 'sir' and coded as 'R' (nonwildtype) or 'S' (wildtype). Only used when downloading phenotype data, with reformat set to `TRUE`.
#' @param release (Optional) String specifying the data release to download (default NULL, will pull latest release).
#' @importFrom arrow read_parquet
#' @importFrom RCurl getBinaryURL getURL
#' @return A data frame containing EBI genotype data
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' # download all phenotype data from EBI
#' pheno_ebi <- download_ebi()
#'
#' # download phenotype data from Dec 2025 release, and filter to Salmonella
#' pheno_salmonella <- download_ebi(
#'     genus="Salmonella",
#'     user_release="2025-12"
#' )
#' 
#' # reformat downloaded phenotype data to simplify use with AMRgen functions
#' pheno_salmonella <- import_ebi_ast_ftp(pheno_salmonella)
#' 
#' # download phenotype data for Staphylococcus aureus and reformat 
#' # for use with AMRgen functions
#' pheno_staph <- download_ebi(
#'     species="Staphylococcus aureus",
#'     reformat=T
#' )
#' 
#' # download phenotype data for Klebsiella quasipneumoniae, reformat 
#' # for use with AMRgen functions, and re-interpret phenotypes
#' pheno_kquasi_reinterpreted <- download_ebi(
#'     species="Klebsiella quasipneumoniae",
#'     reformat=T,
#'     interpret_eucast = TRUE, 
#'     interpret_clsi = TRUE, 
#'     interpret_ecoff = TRUE
#' )
#' 
#' ebi_geno <- download_ebi(data="genotype")
#'     
#' geno_kpn_tmp <- download_ebi(
#'     data="genotype",
#'     species="Klebsiella pneumoniae",
#'     geno_subclass="TRIMETHOPRIM"
#' )
#' }
download_ebi <- function(data="phenotype", antibiotic=NULL, 
                         genus=NULL, species=NULL,
                         geno_subclass=NULL, geno_class=NULL, remove_dup=FALSE,
                         release=NULL, reformat=FALSE,
                         interpret_eucast = FALSE, 
                         interpret_clsi = FALSE, 
                         interpret_ecoff = FALSE) {
  
  # EBI source url
  ebi_url <- "ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/"

  filename <- paste0(data,".parquet")
  
  cat(paste("Downloading", data, "data from EBI AMR portal (ftp://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/)\n"))
  
  if (!is.null(release)){
    user_release <- release
    cat(paste(" Requesting data from user-specified release", user_release, "\n"))
    
  } else {
    # get list of releases
    folders <- str_split(RCurl::getURL(ebi_url, dirlistonly = TRUE), "\n")[[1]]
    # get latest
    user_release <- folders[!folders %in% c("releases.yml", "")] %>% max()
    cat(paste("...Requesting data from latest release", user_release, "\n"))
  }
  
  ebi_dat <- RCurl::getBinaryURL(print(gsub("ftp\\:", "https\\:", paste0(ebi_url,user_release,"/",paste0(data,".parquet")))))
  
  cat(paste("...Reading parquet file\n"))
  ebi_dat <- arrow::read_parquet(ebi_dat)
  
  # filter by genus or species as required
  if (!is.null(genus)){
    user_genus <- genus
    cat(paste("...Filtering by genus", genus, "\n"))
    ebi_dat <- ebi_dat %>% 
      filter(genus==user_genus)
  } else if (!is.null(species)){
    cat(paste("...Filtering by species", species, "\n"))
    user_species <- species
    ebi_dat <- ebi_dat %>% 
      filter(species==user_species)
  }
  
  # filter by subclass/class/antibiotic as required
  if (!is.null(geno_subclass) & data=="genotype"){
    cat(paste("...Filtering by subclass", geno_subclass, "\n"))
    ebi_dat <- ebi_dat %>% 
      filter(grepl(geno_subclass, subclass))
  } else if (!is.null(geno_class) & data=="genotype"){
    cat(paste("...Filtering by class", geno_class, "\n"))
    ebi_dat <- ebi_dat %>% 
    filter(grepl(geno_class, class))
  } else if (!is.null(antibiotic)){
    cat(paste("...Filtering by antibiotic", antibiotic, "\n"))
    ebi_dat <- ebi_dat %>% 
      filter(antibiotic_name==antibiotic)
  }
  
  if (remove_dup) {
    cat("...Checking for duplicate rows ")
    before <- nrow(ebi_dat)
    ebi_dat <- ebi_dat %>% select(BioSample_ID:subclass) %>% distinct()
    after <- nrow(ebi_dat)
    if (after < before) {
      n <- before-after
      cat(paste("(removed", n, "rows)\n"))
    } else {
      cat(paste("(none detected)\n"))
    }
  }
  
  # reformat as per AMRgen import functions
  if (reformat & data=="phenotype") {
    cat(paste("Reformatting phenotype data for easy use with AMRgen functions\n"))
    ebi_dat <- import_ebi_ast_ftp(ebi_dat,
                                  interpret_eucast = interpret_eucast, 
                                  interpret_clsi = interpret_clsi, 
                                  interpret_ecoff = interpret_ecoff)
  }
  
  # return data frame
  return(ebi_dat)
}
