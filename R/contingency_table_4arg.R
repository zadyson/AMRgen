#' Generate Contingency Tables for Genotype and Phenotype Data by Antibiotic
#'
#' This function creates contingency tables comparing genotype and phenotype resistance data 
#' for a given antibiotic or list of antibiotics. It uses the AMR package to standardize antibiotic 
#' names and processes data accordingly. It also requires an antibiotic subclass or list of subclasses 
#' to be specified by the user, which will hopefully be replaced once a mapping function is created (such as as.AFP_subclass()).
#'
#' @param genotype_data A dataframe containing genotype information, harmonised by HAMRonize, with the isolate identified stored in the "input_sequence_id", and the column "antimicrobial_agent" corresponding to AMRFinder/HAMRonize Subclasses.
#' @param phenotype_data A dataframe containing phenotype information, including a `#BioSample` column to identify isolates, a column of class sir, and a column of class ab. 
#' @param antibiotic_list A single antibiotic name or a vector of antibiotic names, in any format
#' @param ab_subclass_list A single antibiotic subclass name or a vector of the same length as antibiotic_list of antibiotic subclasses, in the same order as antibiotic_list. This is a temporary workaround until there is a function that can do this mapping automatically.
#' NOTE: The isolate names in the HAMRonize output may be lost, so the first part of the input_sequence_id column (which identified the contig) is used. 
#' NOTE: I have combined R and I for the final tables to obtain a 2x2 table, and had to convert from class SIR to class character to do this, otherwise I column persisted. 
#' NOTE: Output tables are currently not very pretty
#' 
#' @return A list of contingency tables. for each antibiotic.
#' @examples
#' 
#' genotype_data <- data.frame(input_sequence_id = c("Sample1.contig1", "Sample2.contig1", "Sample3.contig1", "Sample5.contig1"),
#'                             antimicrobial_agent = c("AMIKACIN", "QUINOLONE", "AMIKACIN", "BETA-LACTAM"),
#'                             gene_symbol = c("aph(2'')-IIa", "gyrA", "aph(2'')-IIa", "CTX-M-15"))
#' phenotype_data <- data.frame(`#BioSample` = c("Sample1", "Sample2", "Sample4", "Sample5"),
#'                              Antibiotic = as.ab(c("AMK", "CIP", "AMK", "AMK")),
#'                              Resistance_phenotype = as.sir(c("R", "I", "S", "R")))
#' colnames(phenotype_data)[1] <- "#BioSample"
#' antibiotic_list <- c("cipro", "amickacin")
#' ab_subclass_list <- c("QUINOLONE", "AMIKACIN")
#' contingency_table_by_subclass(genotype_data, phenotype_data, antibiotic_list, ab_subclass_list )
#' 
#' @import dplyr
#' @import tidyverse
#' @import AMR
#' @importFrom glue glue
#' @importFrom purrr purrr
#' 
#' @export
contingency_table_by_subclass <- function(genotype_data, phenotype_data, antibiotic_list, ab_subclass_list) {
  
  #' Convert a single string to a list if necessary
  if (is.character(antibiotic_list) && length(antibiotic_list) == 1) {
    antibiotic_list <- list(antibiotic_list)
  }
  
  #' Normalize antibiotic names using `as.ab`
  antibiotic_list <- purrr::map_chr(antibiotic_list, ~ as.ab(.))
  
  # Select the first "antibiotic" class column from the phenotype data
  antibiotic_column <- names(phenotype_data)[sapply(phenotype_data, inherits, "ab")][1]
  if (is.null(antibiotic_column)) {
    stop("No column of class 'ab' found in phenotype_data.")
  }
  
  # Select the first "SIR" class column from the phenotype data
  sir_column <- names(phenotype_data)[sapply(phenotype_data, inherits, "sir")][1]
  if (is.null(sir_column)) {
    stop("No column of class 'sir' found in phenotype_data.")
  }
  
  # Initialize an empty list to store results
  contingency_tables <- list()
  
  
  #' Process each antibiotic separately
  for (i in seq_along(antibiotic_list)) {
    antibiotic <- antibiotic_list[i]
    antibiotic_subclass <- ab_subclass_list[i]
    
    if (!antibiotic %in% unique(phenotype_data[[antibiotic_column]])) {
      stop(glue::glue("Antibiotic '{antibiotic}' not found in '{antibiotic_column}'"))
    }
    
    
    # Create a new dataframe with isolates from phenotype_data (#BioSample)
    output_df <- data.frame(
      input_isolate = unique(phenotype_data[["#BioSample"]]),
      genotype = 0,
      phenotype = NA,
      stringsAsFactors = FALSE
    )
    
    # Iterate over each isolate in the output dataframe
    for (j in seq_len(nrow(output_df))) {
      isolate <- output_df$input_isolate[j]
      
      # Check if any rows in genotype_data for this isolate contain antibiotic_subclass in the subclass column
      genotype_match <- genotype_data %>%
        filter(input_isolate == isolate) %>%
        filter(grepl(antibiotic_subclass, antimicrobial_agent, ignore.case = TRUE))  
      
      # If a match is found, set genotype to 1, otherwise 0
      if (nrow(genotype_match) > 0) {output_df$genotype[j] <- 1}
      
      
      # Check if phenotype_data has rows for this isolate and the specified antibiotic
      phenotype_match <- phenotype_data %>%
        filter(`#BioSample` == isolate & !!sym(antibiotic_column) == antibiotic) 
      
      if (nrow(phenotype_match) > 0) {
        output_df$phenotype[j] <- phenotype_match[[sir_column]][1]
        #' Ensure sir class variable is handled correctly, otherwise inputed as 1, 2, 3 ofr S, I, R
        output_df$phenotype <- as.sir(output_df$phenotype)
      }
    }
      
    #' Delete following 2 lines to get 3x2 table showinf R and I separately.  
    output_df$phenotype <- as.character(output_df$phenotype)  
    output_df$phenotype <- ifelse(output_df$phenotype == "I", "R", output_df$phenotype)
     
      
      # Generate contingency table
      table_basic <- table(output_df$genotype, output_df$phenotype)
      table_combined <- addmargins(table_basic)
      
      # Add percentages
      row_percent <- round(100 * prop.table(table_basic, margin = 1), 1)
      col_percent <- round(100 * prop.table(table_basic, margin = 2), 1)
      
      # Format percentages with "%" symbol
      row_percent <- apply(row_percent, c(1, 2), function(x) paste0(x, "%"))
      col_percent <- apply(col_percent, c(1, 2), function(x) paste0(x, "%"))
      
      # Store results in a list
      contingency_tables[[antibiotic]] <- list(
        counts_table = table_combined,
        row_percent = row_percent,
        col_percent = col_percent
      )
    }
    
    return(contingency_tables)
  }

