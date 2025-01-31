require(ComplexUpset)

amr_complexUpset <- function(binary_matrix, min_set_size = 10, mic_disk = 'mic', 
                      remove_NAs = TRUE, gene_determinants = NULL, colour_by='pheno',
                      plot_breakpoints=FALSE, organism=NULL, break_guide = "EUCAST 2024", 
                      break_type="ECOFF", drug=NULL, colour_values = c("#66c2a5", "#fdae61", "#d53e4f")) {
  
  # mic_disk must be either 'mic' or 'disk'
  if (! mic_disk %in% c('mic', 'disk')){
    stop("mic_disk must be either 'mic' or 'disk'. Please select one of these values to continue.")
  }
  
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
    if (is.null(organism) || is.null(break_guide) || is.null(break_type) || is.null(drug)) {
      stop("When plot_breakpoints is TRUE, all of the following arguments must have values: organism, guideline, type, drug.")
    }
    # grab the method to filter by depending on whether we're using mic or disk
    if (mic_disk == 'mic'){
      method_value <- "MIC"
    }
    else{ method_value <- "DISK" }
    # extract the breakpoints
    breakpoint_table <- AMR::clinical_breakpoints %>% 
      filter(ab == drug, mo == organism, guideline == break_guide, type == break_type, method == method_value)
    # if the table is empty, then return an error
    if(nrow(breakpoint_table) == 0){
      warning("No breakpoints available for bug/drug combination. Will not plot breakpoint lines on upset, continuing...")
    }
    # otherwise grab the values, and create the hlines
    else {
      break_r <- breakpoint_table %>% pull(breakpoint_R)
      break_s <- breakpoint_table %>% pull(breakpoint_S)
      # convert to mic values if needed
      if (mic_disk == 'mic'){
        break_r <- as.mic(break_r)
        break_s <- as.mic(break_s)
      }
      if (break_r == break_s){
        warning("Breakpoint value is the same for R and S, only plotting one line")
        break_r_line <- geom_hline(yintercept = break_r, linetype = "dashed", alpha = 0.6, colour="black")
        break_s_line <- NULL
      }
      else{
        break_r_line <- geom_hline(yintercept = as.mic(break_r), linetype = "dashed", alpha=0.6, colour="black")
        break_s_line <- geom_hline(yintercept = as.mic(break_s), linetype = "dashed", alpha = 0.6, colour="black")
      }
    }
  }
  else {
    break_r_line <- NULL
    break_s_line <- NULL
  }
    
  # remove NAs if required
  if (remove_NAs){
    upset_data <- upset_data %>% filter(!is.na(!!sym(mic_disk)))
  }
    
  # get the name of the y axis
  if(mic_disk == 'mic'){
    y_axis_name <- "MIC (mg/L)"
    scale_y <- scale_y_mic()
  }
  else{
    y_axis_name <- "Disk (mm)"
    scale_y <- NULL
  }
  
  # now plot
  plot <- upset(upset_data, genes, name="genetic determinant", min_size=min_set_size, width_ratio = 0.1,
  annotations = list(
    y_axis_name=(ggplot(mapping=aes(y=!!sym(mic_disk), colour=as.factor(!!sym(colour_by)))) + 
      geom_count() +
      break_r_line + # these values will be NULL if plot breakpoints isn't set
      break_s_line +
      scale_colour_manual(values=colour_values) +
      scale_y +
      labs(colour=colour_by, size="Number of\nisolates") +
      theme(legend.title = element_text(face="bold"))
  )))

  print(break_r, break_s)
  return(plot)
}