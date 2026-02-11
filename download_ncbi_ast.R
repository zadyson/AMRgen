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
#' Download NCBI antimicrobial susceptibility testing (AST) data
#'
#' This function downloads antimicrobial susceptibility testing (AST) data
#' from the NCBI Pathogen Detection database via the BioSample API. Data are
#' retrieved in batches, parsed from XML, and returned as a tidy tibble with
#' metadata including BioSample ID, Bioproject ID, and organism name.
#'
#' @param user_organism Character. Organism name for the search query 
#'   (e.g., `"Salmonella enterica"`). Required.
#' @param user_antibiotic Character. Optional antibiotic name to filter the
#'   returned data. Matching is case-insensitive.
#' @param user_max_records Integer. Maximum number of BioSample records to
#'   retrieve. Default is `10000`.
#' @param user_batch_size Integer. Number of records fetched per API request.
#'   Default is `200` which is recommended by NCBI.
#' @param user_sleep_time Numeric. Seconds to pause between batch requests to
#'   avoid overloading NCBI servers. Default is `0.34`.
#' @param user_force_antibiotic Logical. If `TRUE`, standardizes the antibiotic
#'   name using the \pkg{AMR} package before filtering. Default is `FALSE`.
#' @param user_reformat Logical. If `TRUE`, reformats the output using
#'   `AMRgen::import_ncbi_ast()` for compatibility with AMR analysis workflows.
#'   Default is `FALSE`.
#' @param user_interpret_eucast Logical. Passed to `import_ncbi_ast()`.
#'   If `TRUE`, interprets MIC values using EUCAST breakpoints. Default is
#'   `FALSE`.
#' @param user_interpret_clsi Logical. Passed to `import_ncbi_ast()`.
#'   If `TRUE`, interprets MIC values using CLSI breakpoints. Default is
#'   `FALSE`.
#' @param user_interpret_ecoff Logical. Passed to `import_ncbi_ast()`.
#'   If `TRUE`, interprets MIC values using ECOFF cutoffs. Default is `FALSE`.
#'
#' @details
#' The function constructs an Entrez query of the form:
#' `"<organism> AND antibiogram[filter]"`. XML records are downloaded
#' in batches, parsed, and combined into a single table. The resulting tibble
#' contains AST test results and associated metadata including:
#' 
#' - `id`: BioSample identifier
#' - `bioproject`: Bioproject accession ID
#' - `organism`: Organism name
#' - `Antibiotic`, `Phenotype`, `Measurement`, `Units`, `Method`, `System`,
#'   `Manufacturer`, `Panel`, `Standard`: AST data columns
#'
#' The function can optionally filter by a single antibiotic, and reformat
#' data for compatibility with AMRgen functions for antimicrobial resistance
#' analysis.
#'
#' @return
#' A tibble containing AST results and metadata.
#'
#' @examples
#' \dontrun{
#' # Download Salmonella AST data
#' ast <- download_ncbi_ast("Salmonella enterica")
#'
#' # Filter to ciprofloxacin
#' ast <- download_ncbi_ast(
#'   "Salmonella enterica",
#'   user_antibiotic = "ciprofloxacin",
#'   user_force_antibiotic = TRUE
#' )
#'
#' # Reformat for AMRgen workflow with EUCAST interpretation
#' ast <- download_ncbi_ast(
#'   "Escherichia coli",
#'   user_reformat = TRUE,
#'   user_interpret_eucast = TRUE
#' )
#' }
#'
#' @import rentrez
#' @import dplyr
#' @import tibble
#' @import XML
#' @importFrom AMR ab_name as.ab
#' @export
#' 
download_ncbi_ast <- function(user_organism = NULL, 
                              user_antibiotic = NULL,
                              user_max_records = 10000,
                              user_batch_size = 200,
                              user_sleep_time = 0.34,
                              user_force_antibiotic = FALSE,
                              user_reformat = FALSE,
                              user_interpret_eucast = FALSE,
                              user_interpret_clsi = FALSE,
                              user_interpret_ecoff = FALSE){
  
  # Build query term for entrez
  enterz_term <- paste0(tolower(user_organism), " AND antibiogram[filter]")
  
  
  # Search for ATS entries by pathogen 
  search <- rentrez::entrez_search(
    db = "biosample",
    term = enterz_term,
    retmax = user_max_records
  )
  
  # list samples by internal id
  ids <- search$ids
  
  cat(paste("Idenfified", length(ids), " ", tolower(user_organism), 
            " records from NCBI (https://www.ncbi.nlm.nih.gov/pathogens/ast/) \n"))
  
  # create empty data structures to store records
  all_ast_data <- NULL
  
  # Pull records in batches
  batch_size <- user_batch_size
  
  # Pause time between data pulls - treat server nicely
  sleep_time <- user_sleep_time
  
  # iterate though samples and pull records
  for (sample in seq(1, length(ids), by = batch_size)) {
    
    # account for remainder when batching
    if ((sample-1+batch_size > length(ids))){
      
      data_xml <- rentrez::entrez_fetch(db="biosample", 
                                        id=ids[sample:(length(ids))], 
                                        rettype="xml", 
                                        parsed=TRUE)   
      
      cat(paste("Downloading and processing records ", sample, " to ", 
                length(ids), "... \n"))
      
    }else{
      
      data_xml <- rentrez::entrez_fetch(db="biosample", 
                                        id=ids[sample:(sample-1+batch_size)], 
                                        rettype="xml", 
                                        parsed=TRUE)
      
      cat(paste("Downloading and processing records ", sample, " to ", 
                (sample-1+batch_size), "... \n"))
      
    }
    
    # Flatten XML
    data_list <- XML::xmlToList(data_xml)

    # extract records
    for (record in 1:length(data_list)){
      
      data_names <- as.character(unlist(
        data_list[record]$BioSample$Description$Comment$Table$Header))
      
      data_entries <- data_list[record]$BioSample$Description$Comment$Table$Body
      
      # retrieve AST data and reformat
      for (entry in 1:length(data_entries)){
        
        temp_entry <- data_list[record]$BioSample$Description$Comment$Table$Body[entry]
        
        # Preserve null data points by converting to NA (for headers)
        temp_entry <- temp_entry %>% 
          purrr::modify_tree(
            leaf = \(x) if(is.null(x)) NA else x,
            post = unlist
          ) %>% 
          t() %>%
          as.data.frame()
        
        # add col names here to prevent data mismatches when binding rows
        names(temp_entry) <- c(as.character(data_names))
        
        # Extract BioProject accession while accounting for two varying flattened 
        # xml node structures with duplicated identifiers
        bioproj <- pluck(data_list[record], "BioSample", "Links", 2, "Link", ".attrs", 3, .default = NA)
        
        if (is.na(bioproj)){
          bioproj <- pluck(data_list[record], "BioSample", "Links", "Link", ".attrs", 3, .default = NA)
        }else{
          bioproj <- NA # when no data listed - very rare
        }
        
        # Add accessions and organism cols to AST data
        temp_entry <- temp_entry %>%
          mutate(id = data_list[record]$BioSample$Ids$Id$text) %>%
          mutate(BioProject = bioproj) %>%
          mutate(organism = data_list[record]$BioSample$Description$Organism$OrganismName)
        
        # combine records
        all_ast_data <- bind_rows(all_ast_data, temp_entry)
        
      }  
    }
    
    # pause between record pulls
    Sys.sleep(sleep_time)
    cat(paste("Pausing download for ", sleep_time, "seconds to avoid overburdening the NCBI server... \n"))
    
  }
  
  # name and rearrange cols
  cat(paste("Ordering data... \n"))
  
  all_ast_data <- all_ast_data %>%
    select(id, BioProject, organism, everything()) %>%
    as_tibble()
  
  # filter data by drug where specified by user
  if (!is.null(user_antibiotic)){
    
    if (user_force_antibiotic){
      
      user_antibiotic <- na.omit(tolower(AMR::ab_name(AMR::as.ab(user_antibiotic))))
    }
    
    cat(paste("Filtering data by antibiotic", user_antibiotic, "... \n"))
    
    all_ast_data <- all_ast_data %>% 
      filter(Antibiotic == user_antibiotic)
  }
  
  # reformat as per AMRgen import functions
  if (user_reformat) {
    
    cat(paste("Reformatting phenotype data for easy use with AMRgen functions \n"))
    
    all_ast_data <- AMRgen::import_ncbi_ast(all_ast_data,
                                            sample_col = "id",
                                            interpret_eucast = user_interpret_eucast,
                                            interpret_clsi = user_interpret_clsi,
                                            interpret_ecoff = user_interpret_ecoff)
  }
  
  # return data 
  return(all_ast_data)
  
}
