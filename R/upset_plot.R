require(ComplexUpset)

upset_plot <- function(binary_matrix, min_set_size = 10, mic_disk = 'mic', 
                      remove_NAs = TRUE, gene_determinants = NULL, colour_by='pheno',
                      plot_breakpoints=FALSE, organism=NULL, guideline = "EUCAST 2024", 
                      type="ECOFF", drug=NULL) {
  
  # Convert organism and antibiotic to mo and ab if they are provided
  if (!is.null(organism)) {
    organism <- as.mo(organism)
  }
  if (!is.null(drug)) {
    drug <- as.ab(drug)
  }
  
  # select gene determinants to plot
  if(is.null(gene_determinants)){
    col <- colnames(binary_matrix)
    cols_to_remove <- c("mic", "disk", "R", "NWT", "pheno")
    genes <- col[-1]
    genes <- setdiff(genes, cols_to_remove)
  }
  # otherwise set to the vector provided
  else{
    genes <- gene_determinants
  }

  # convert gene columns to logical
  upset_data <- binary_matrix %>% mutate(across(genes, as.logical))

  # rename first column to 'sample' for the first column, to avoid issues with duplicate column names, if col1 is "id"
  if(colnames(upset_data)[1] == "id"){
    upset_data <- upset_data %>% rename(sample = id)
  }

  # if plot breakpoints is set, then get the relevant info (mo, ab, type and guideline MUST be set)
  if (plot_breakpoints) {
    # Check if all required arguments have values
    if (is.null(organism) || is.null(guideline) || is.null(type) || is.null(drug)) {
      stop("When plot_breakpoints is TRUE, all of the following arguments must have values: organism, guideline, type, drug.")
    }

  }
  

  # grab the column to plot on the y axis - should be either mic or disk, so we can set the scale correctly
  if(mic_disk == 'mic'){
    if(remove_NAs){
      upset_data <- upset_data %>% filter(!is.na(mic))
    }
    plot <- upset(upset_data, genes, name="genetic determinant", min_size=min_set_size, width_ratio = 0.1,
  annotations = list(
    'MIC (mg/L)'=(ggplot(mapping=aes(y=mic)) + 
      geom_jitter(aes(colour=!!sym(colour_by))) + 
      #geom_hline(yintercept = as.mic(breakpoint_value)) +
      scale_y_mic()
  )))
  }
  if(mic_disk == 'disk'){
    if(remove_NAs){
      upset_data <- upset_data %>% filter(!is.na(disk))
    }
    plot <- upset(upset_data, genes, name="genetic determinant", min_size=min_set_size, width_ratio = 0.1,
  annotations = list(
    'disk (mm)'=(ggplot(mapping=aes(y=disk)) + 
      geom_jitter(aes(colour=!!sym(colour_by)))
  )))
  } 

  return(plot)
}