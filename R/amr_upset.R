### Function to generate "amr_upset"

#' amr_upset: Generate a series of plots for AMR gene and combination analysis
#'
#' This function generates a set of visualizations to analyze AMR gene combinations, MIC values, and gene prevalence
#' from a given binary matrix. It creates several plots, including MIC distributions, a bar plot for
#' the percentage of strains per combination, a dot plot for gene combinations in strains, and a plot for gene prevalence.
#' It also outputs a table summarising the MIC distribution (median, lower, upper) and number resistant, for each marker combination.
#'
#' @param binary_matrix A data frame containing the original binary matrix output from the `get_binary_matrix` function.
#'        Expected columns are and identifier (column 1, any name); 'pheno' (class sir, with S/I/R categories to colour points),
#'        'mic' (class mic, with MIC values to plot), and other columns representing gene presence/absence (binary coded, ie 
#'        1=present, 0=absent).
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots.
#'        Default is 2. Only genes with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#'        - "": Default (decreasing frequency of combinations),
#'        - "genes": Order by the number of genes in each combination.
#'        - "value": Order by the median assay value (MIC or disk zone) for each combination.
#' @param plot_set_size Logical indicating whether to include a bar plot showing the set size (i.e. 
#'        number of times each combination of markers is observed). Default is FALSE.
#' @param print_set_size Logical indicating whether, if plot_set_size is set to true, to print the number of strains
#'        with marker combination on the plot (default FALSE).
#' @param plot_category Logical indicating whether to include a stacked bar plot showing, for each marker combination,
#'         the proportion of samples with each phenotype classification (specified by the 'pheno' column in the input file). 
#'         Default is TRUE.
#' @param print_category_counts Logical indicating whether, if plot_category is set to true, to print the number of strains
#'        in each resistance category, for each marker combination in the plot (default FALSE).
#' @param boxplot_colour Colour for lines of the box plots summarising the MIC distribution for each marker combination,
#'        (default "grey").
#' @param assay A character string indicating whether to plot MIC or disk diffusion data.
#'        - "mic": Plot MIC data, stored in column 'mic' of class 'mic'.
#'        - "disk": Plot disk diffusion data, stored in column 'disk' of class 'disk'.

#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{plot}{A grid of plots displaying: (i) grid showing the marker combinations observed, MIC distribution per marker combination, frequency per marker and (optionally) phenotype classification and/or number of samples for each marker combination.}
#'     \item{summary}{Summary of each marker combination observed, including median MIC (and interquartile range) and positive predictive value for resistance (R).}
#'   }
#'
#' @details This function processes the provided binary matrix (`binary_matrix`), which is expected to contain data on gene 
#'          presence/absence, MIC values, and phenotype calls (S/I/R) (can be generated using `get_binary_matrix`).
#'          The function also includes an analysis of gene prevalence and an ordering option for visualizing combinations
#'          in different ways.
#'
#' @examples
#' \dontrun{
#' # Example usage
#'
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#' 
#' binary_matrix<- get_binary_matrix(geno_table=ecoli_geno, 
#'               pheno_table=ecoli_ast, 
#'               antibiotic="Ciprofloxacin", 
#'               drug_class_list=c("Quinolones"), 
#'               sir_col="pheno", 
#'               keep_assay_values=TRUE, 
#'               keep_assay_values_from = "mic"
#'            )
#' 
#' amr_upset(binary_matrix, min_set_size = 3, order = "value", assay="mic")
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom forcats fct_rev
#' @importFrom AMR as.mic
#' @import patchwork
#' 
#' @export
amr_upset <- function(binary_matrix, min_set_size = 2, order = "", 
                      plot_set_size=FALSE, plot_category=TRUE, 
                      print_category_counts=FALSE, print_set_size=FALSE,
                      boxplot_colour="grey", assay="mic") {

  # tidy up binary_matrix
  col <- colnames(binary_matrix) # get column names 
  
  # extract only the gene column names - need to exclude mic, disk, R, NWT (standard col names)
  # and the id column which will be the first col, doesn't matter what it's called
  # remaining columns will be the genes
  cols_to_remove <- c("mic", "disk", "R", "NWT", "pheno")
  genes <- col[-1]

  # gene names
  genes <- setdiff(genes, cols_to_remove)
  
  # check with have the expected assay column, with data
  if (!(assay %in% colnames(binary_matrix))) {
    stop(paste("input", deparse(substitute(binary_matrix)), "must have a column labelled ", assay))
  }
  data_rows <- binary_matrix %>% filter(!is.na(get(assay))) %>% nrow()
  if (data_rows==0) {
    stop(paste("input", deparse(substitute(binary_matrix)), "has no non-NA values in column ", assay))
  }

  # Add in a combination column and filter to samples with the required assay data (MIC or disk) only
  binary_matrix_wide <- binary_matrix %>% 
    filter(!is.na(get(assay))) %>%
    unite("combination_id", genes[1]:genes[length(genes)], remove = FALSE)  # add in combinations 
  
  # Make matrix longer 
  binary_matrix <- binary_matrix_wide %>% pivot_longer(cols = genes[1]:genes[length(genes)], names_to = "genes")
  
  ### Counts per combination, for bar plot - X axis = combination. Y axis = number of strains ###
  # This first to filter on combinations with enough data 
  combination_freq <- binary_matrix_wide %>% group_by(combination_id) %>%
    summarise(n = n()) %>%
    mutate(perc = 100 * n / sum(n)) %>% # count number with each combination
    filter(n >= min_set_size) # count filter
  
  if (nrow(combination_freq)==0) {stop(paste("No marker combinations pass the minimum frequency:", min_set_size))}
  # which have enough strains / data
  comb_enough_strains <- combination_freq %>% pull(combination_id)

  ### data for assay plot - dot plot. X axis = combination. Y axis = MIC or disk zone #####
  assay_plot <- binary_matrix %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(combination_id, get(assay), pheno) %>%
    summarise(n = n())  # count how many at each assay value, keep pheno for colour
  
  ### Gene prevalence plot (amongst strains with marker combinations exceeding the minimum set size)
  gene.prev <- binary_matrix %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(genes) %>%
    summarise(gene.prev = sum(value))

  ### Set order of genes for y axis in combination dot plot 
  ## (amongst strains with marker combinations exceeding the minimum set size, that will be included in plot)
  gene.order.desc <- gene.prev %>%
    arrange(desc(gene.prev)) %>%
    filter(gene.prev > 0) %>%
    pull(genes)

  # For gene prev plot
  gene.prev <- gene.prev %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc))
  
  # For combination dot plot 
  binary_matrix <- binary_matrix %>%
    filter(combination_id %in% comb_enough_strains) %>% 
    filter(genes %in% gene.order.desc) %>% 
    mutate(genes = factor(genes, levels = gene.order.desc))
  
  ### Point plot - X axis = combination. Y axis = genes. Lines joining genes in same strain ### 
  binary_matrix$point_size = 2 * binary_matrix$value # want dot size to be larger than 1 => can make 2/3/4 etc
  
  # Only plot a line between points if more than one gene in a combination 
  # get how many in each combination 
  multi_genes_combination_id_all <- binary_matrix %>%
    select(combination_id, genes, value, point_size) %>%
    group_by(combination_id) %>%
    filter(value == 1) %>%
    mutate(u = length(unique(genes))) %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc, ordered=T))
  
  # get only those with > 1
  multi_genes_combination_ids <- multi_genes_combination_id_all %>%
    filter(u > 1) %>%
    group_by(combination_id) %>%
    mutate(
      min = min(genes),
      max = max(genes)
    ) %>%
    ungroup()

  ### Set order of combination_id <- x axis
  # Default = decreasing frequency
  ordered_comb_order <- combination_freq %>%
    arrange(desc(perc)) %>%
    pull(combination_id)
  assay_plot$combination_id <- factor(assay_plot$combination_id, levels = ordered_comb_order)
  combination_freq$combination_id <- factor(combination_freq$combination_id, levels = ordered_comb_order)
  binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  
  # Do by # genes in combination (only want each id once)
  if (order == "genes") {
    ordered_comb_order <- multi_genes_combination_id_all %>%
      arrange(u) %>%
      filter(row_number() == 1) %>%
      pull(combination_id)
    assay_plot$combination_id <- factor(assay_plot$combination_id, levels = ordered_comb_order)
    combination_freq$combination_id <- factor(combination_freq$combination_id, levels = ordered_comb_order)
    binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  }
  # Do by # median assay value in combination (only want each id once)
  if(order == "value"){
    if (assay=="mic") {
      mic_medians <- binary_matrix_wide %>% group_by(combination_id) %>% summarise(median = median(mic))
      ordered_comb_order <- mic_medians %>% arrange(median) %>% pull(combination_id)
    }
    else {
      disk_medians <- binary_matrix_wide %>% group_by(combination_id) %>% summarise(median = median(as.double(disk)))
      ordered_comb_order <- disk_medians %>% arrange(-median) %>% pull(combination_id)
    }
    assay_plot$combination_id <- factor(assay_plot$combination_id, levels = ordered_comb_order)
    combination_freq$combination_id <- factor(combination_freq$combination_id, levels = ordered_comb_order)
    binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  }
  
  ##### Plots ###

  ### assay data plot (MIC/disk)
  g1 <- ggplot(data = assay_plot, aes(x=combination_id, y=`get(assay)`)) +
    geom_boxplot(colour = boxplot_colour) +
    geom_point(aes(size = n, colour = pheno), show.legend = TRUE) +
    theme_bw() +
    scale_size_continuous("Number of\nisolates") +
    scale_color_sir(name="Resistance\ncategory") +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  if (assay=="mic") {
    g1 <- g1 + 
      scale_y_mic() +
      ylab("MIC (mg/L)")
  }
  else {
    g1 <- g1 + 
      ylab("Disk zone (mm)")
  }
  
  ### Bar plot - set size 
  g2 <- ggplot(combination_freq, aes(x=combination_id, y = n)) + 
    geom_bar(stat = "identity") + 
    theme_bw() +
    scale_y_reverse("Set size") + 
    scale_x_discrete("group") +
    theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
  if (print_set_size) {
    g2 <- g2 + geom_text(aes(label = n), nudge_y = -.5)
  }
  
  # pheno category stacked barplot
  category_plot <- binary_matrix %>% select(id, pheno, combination_id) %>% distinct() %>%
    ggplot(aes(x=combination_id, fill=pheno)) +
    geom_bar(stat = "count", position = "fill") +
    scale_fill_sir() +
    theme_light() +
    labs(x = "", y = "Category") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none")
  if (print_category_counts) {
    category_plot <- category_plot + 
      geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(vjust = .5), size = 3)
  }
    
  ### Dot plot of combinations 
  g3 <- binary_matrix %>% 
    mutate(binary_comb=if_else(value>0, 1, 0)) %>%
    ggplot(aes(x=combination_id, y=fct_rev(genes))) + 
    geom_point(aes(size = binary_comb), show.legend = FALSE) + 
    theme_bw() + 
    scale_size_continuous(range = c(-1,2)) + 
    scale_y_discrete("Marker") + 
    geom_segment(data = multi_genes_combination_ids,
                 aes(x = combination_id, xend = combination_id, 
                     y = min, yend = max, group = combination_id),
                 color = "black") + # add lines
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

  ### Plot gene prev / set size
  g4 <- ggplot(gene.prev, aes(x = fct_rev(genes), y = gene.prev)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    scale_y_reverse("") +
    theme(axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank())
  
  # assemble plot
  final_plot <- patchwork::plot_spacer() + g1 + patchwork::plot_layout(ncol = 2, widths = c(1, 4), guides="collect")
  if (plot_category) { final_plot <- final_plot + patchwork::plot_spacer() + category_plot }
  final_plot <- final_plot + g4 + g3
  if (plot_set_size) { final_plot <- final_plot + patchwork::plot_spacer() + g2 }

  # set relative plotting heights
  if (plot_category & plot_set_size) {final_plot <- final_plot + patchwork::plot_layout(heights=c(2,1,2,1))}
  if (plot_category & !plot_set_size) {final_plot <- final_plot + patchwork::plot_layout(heights=c(2,1,2))}
  if (!plot_category & plot_set_size) {final_plot <- final_plot + patchwork::plot_layout(heights=c(2,2,1))}
  
  print(final_plot)
  
  # summary table
  summary <- binary_matrix_wide %>% 
    group_by(combination_id) %>%
    summarise(median = median(as.double(get(assay))), 
              q25=stats::quantile(as.double(get(assay)),0.25),
              q75=stats::quantile(as.double(get(assay)),0.75),
              ppv=mean(R, na.rm=T),
              R=sum(R, na.rm=T),
              n=n())
  
  # get names for summary
  combination_names <- binary_matrix_wide %>% 
    select(combination_id, any_of(genes)) %>% 
    distinct()
  
  combination_names <- combination_names %>% 
    mutate(marker_list = apply(., 1, function(row) {
      paste(names(combination_names)[-1][row[-1] == 1], collapse = ", ")
    }), .after=combination_id) %>%
    mutate(marker_count = rowSums(. == 1)) %>%
    select(combination_id, marker_list, marker_count)
  
  summary <- summary %>% 
    left_join(combination_names, by="combination_id") %>%
    select(-combination_id) %>%
    mutate(marker_list=if_else(is.na(marker_list), "-", marker_list)) %>%
    mutate(marker_count=if_else(is.na(marker_count), 0, marker_count)) %>% 
    relocate(marker_list, .before=median) %>% 
    relocate(marker_count, .before=median)
  
  return(list(plot=final_plot, summary=summary))
}
