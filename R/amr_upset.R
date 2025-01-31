### Function to generate "amr_upset"

#' amr_upset: Generate a series of plots for AMR gene and combination analysis
#'
#' This function generates a set of visualizations to analyze AMR gene combinations, MIC values, and gene prevalence
#' from a given binary matrix. It creates several plots, including MIC distributions, a bar plot for
#' the percentage of strains per combination, a dot plot for gene combinations in strains, and a plot for gene prevalence.
#'
#' @param binmat_orig A data frame containing the original binary matrix output from the `get_binary_matrix` function,
#'        with columns representing genes, resistance, MIC values, and metadata such as microorganism and antibiotic
#'        information. This needs to be updated / standardised in future versions of AMRgen. 
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots.
#'        Default is 2. Only genes with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#'        - "": Default (decreasing frequency of combinations),
#'        - "genes": Order by the number of genes in each combination,
#'        - "mic": Order by the median MIC of each combination. Default is decreasing frequency.
#'
#' @return A grid of plots displaying:
#'         - A plot for MIC values across combinations (g1),
#'         - A bar plot showing the percentage of strains in each combination (g2),
#'         - A dot plot of gene combinations in strains (g3),
#'         - A gene prevalence plot displaying the set size for each gene (g4).
#'
#' @details This function processes the provided binary matrix (`binmat_orig`), which is expected to contain data on gene 
#'          combinations found in the strain, MIC for that strain (e.g., resistance or susceptibility) and 
#'          corresponding MIC values for different genes.
#'          The function also includes an analysis of gene prevalence and an ordering option for visualizing combinations 
#'          in different ways.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- get_binary_matrix(geno_table, pheno_table, antibiotic="Ciprofloxacin", drug_class_list=c("Quinolones"), sir_col="Resistance phenotype", keep_assay_values=T, keep_assay_values_from = "mic") # Example input data from `get_binary_matrix` function, has mic, R/NWT, and gene lists
#' amr_upset(data, min_set_size = 3, order = "mic")
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom forcats fct_rev
#' @importFrom AMR as.mic
#' @import patchwork

amr_upset <- function(binmat_orig, min_set_size = 2, order = ""){
  ## Inputs
  # takes in binmat = output from get_binary_matrix function (in function.R on datacuration git)
  # takes in order <- single value. Default = decreasing frequency. 
  #         "genes" = # genes. "mic" = median mic 
  # default min set size is 2 (greater than or equal to this)
  
  # tidy up binmat
  col <- colnames(binmat_orig) # get column names 
  # Extract gene names only - not "resistance", "MIC" or "DD" values, or microorganism or antibiotic names
  # might need to change this if want to be more flexible on input? or input this list? 
  
  # extract only the gene column names - need to exclude mic, disk, R, NWT (standard col names)
  # and the id column which will be the first col, doesn't matter what it's called
  # remaining columns will be the genes
  cols_to_remove <- c("mic", "disk", "R", "NWT", "pheno")
  genes <- col[-1]

  # gene names
  genes <- setdiff(genes, cols_to_remove)

  # Add in a combination column 
  binmat_wide <- binmat_orig %>% 
    filter(!is.na(mic)) %>% # make this optionally disk also
    unite("combination_id", genes[1]:genes[length(genes)], remove = FALSE)  # add in combinations 
  # Make matrix longer 
  binmat <- binmat_wide %>% pivot_longer(cols = genes[1]:genes[length(genes)], names_to = "genes")
  
  ##### Data wrangling for plots ###
  ### Bar plot - X axis = combination. Y axis = number of strains ###
  # This first to filter on combinations with enough data 
  bar_plot <- binmat_wide %>% group_by(combination_id) %>%
    summarise(n = n()) %>%
    mutate(perc = 100 * n / sum(n)) %>% # count number with each combination 
    filter(n > 1) # count filter
  # which have enough strains / data 
  comb_enough_strains <- bar_plot %>% pull(combination_id)
  
  ### MIC plot - dot plot. X axis = combination. Y axis = MIC #####
  mic_plot <- binmat %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(combination_id, mic, R) %>%
    summarise(n = n()) # count how many at each MIC, keep resistant for colour
  
  ## TODO: currently broken as ab and mo not in the binary matrix
  ## Extract cutoff from AMR package to plot hline
  #cut_dat <- binmat_orig %>%
  #  dplyr::select(ab, mo) %>%
  #  dplyr::distinct() %>%
  #  dplyr::left_join(AMR::clinical_breakpoints) %>%
    # Magic values that ideally would not be hard coded
  #  dplyr::filter(method == "MIC", host == "human", guideline == "EUCAST 2024") 
  
  ### Gene prevalence plot 
  gene.prev <- binmat %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(genes) %>%
    summarise(gene.prev = sum(value))
  
  ################### ORDER y axis ########################
  ### Set order of genes <- y axis in dot plot 
  ## And filter out on set size mininum
  gene.order.desc <- gene.prev %>%
    arrange(desc(gene.prev)) %>%
    filter(gene.prev >= min_set_size) %>%
    pull(genes)
  
  # For gene prev plot 
  gene.prev <- gene.prev %>%
    filter(genes %in% gene.order.desc) %>% 
    mutate(genes = factor(genes, levels = gene.order.desc))
  
  # For combination dot plot 
  binmat <- binmat %>%
    filter(combination_id %in% comb_enough_strains) %>% 
    filter(genes %in% gene.order.desc) %>% 
    mutate(genes = factor(genes, levels = gene.order.desc))
  
  ############# Which have lines between in dot plot? 
  ### Point plot - X axis = combination. Y axis = genes. Lines joining genes in same strain ### 
  binmat$point_size = 2 * binmat$value # want dot size to be larger than 1 => can make 2/3/4 etc
  # Only plot a line between points if more than one gene in a strain 
  # get how many in each strain 
  multi_genes_combination_id_all <- binmat %>%
    group_by(combination_id) %>%
    filter(value == 1) %>%
    mutate(u = length(unique(genes))) %>% 
    filter(genes %in% gene.order.desc) %>% 
    mutate(genes = factor(genes, levels = gene.order.desc))
  # get only those with > 1
  multi_genes_combination_ids <- multi_genes_combination_id_all %>%
    filter(u > 1) %>%
    mutate(min = first(genes),
           max = last(genes)) %>%
    ungroup() 
  
  ################### ORDER x axis ########################
  
  
  ### Set order of combination_id <- x axis 
  # Default = decreasing frequency 
  ordered_comb_order <- bar_plot %>% arrange(desc(perc)) %>% pull(combination_id)
  mic_plot$combination_id <- factor(mic_plot$combination_id, levels = ordered_comb_order)
  bar_plot$combination_id <- factor(bar_plot$combination_id, levels = ordered_comb_order)
  binmat$combination_id <- factor(binmat$combination_id, levels = ordered_comb_order)
  # Do by # genes in combination (only want each id once)
  if(order == "genes"){
    ordered_comb_order <- multi_genes_combination_id_all %>% arrange(u) %>% filter(row_number()==1) %>% pull(combination_id)
    mic_plot$combination_id <- factor(mic_plot$combination_id, levels = ordered_comb_order)
    bar_plot$combination_id <- factor(bar_plot$combination_id, levels = ordered_comb_order)
    binmat$combination_id <- factor(binmat$combination_id, levels = ordered_comb_order)}
  # Do by # median mic in combination (only want each id once)
  if(order == "mic"){
    mic_medians <- binmat_wide %>% group_by(combination_id) %>% summarise(median = median(mic))
    ordered_comb_order <- mic_medians %>% arrange(median) %>% pull(combination_id)
    mic_plot$combination_id <- factor(mic_plot$combination_id, levels = ordered_comb_order)
    bar_plot$combination_id <- factor(bar_plot$combination_id, levels = ordered_comb_order)
    binmat$combination_id <- factor(binmat$combination_id, levels = ordered_comb_order)}
  
  
  
  ##### Plots ###
  ### AMR package colours
  colours_SIR <- c("#3CAEA3", "#F6D55C", "#ED553B")
  # currently only 0/1
  #names(colours_SIR) <- c("S", "I", "R")
  
  ### MIC plot  
  g1 <- ggplot(data = mic_plot, aes(combination_id, mic)) +
    geom_point(aes(size = n, colour = as.factor(R)), show.legend = TRUE) +
    #geom_hline(data = cut_dat, aes(yintercept = AMR::as.mic(breakpoint_S)), colour = colours_SIR["S"]) +
    #geom_hline(data = cut_dat, aes(yintercept = AMR::as.mic(breakpoint_R)), colour = colours_SIR["R"]) +
    theme_bw() +
    scale_y_mic() +
    ylab("Phenotype (MIC, mg/L)") +
    scale_x_discrete("group") +
    scale_size_continuous("Number of\nisolates") +  
    scale_color_manual("Resistance\nclass", values = colours_SIR) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ### Bar plot  
  g2 <- ggplot(bar_plot, aes(x=combination_id, y = perc)) + 
    geom_bar(stat = "identity") + 
    theme_bw() +
    #geom_text(aes(label = n), nudge_y = -.5) + 
    scale_y_continuous("Percentage") + 
    scale_x_discrete("group") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

  ### Dot plot of combinations 
  g3 <- ggplot(data=binmat, aes(combination_id, fct_rev(genes))) + 
    geom_point(aes(size = value), show.legend = FALSE) + 
    theme_bw() + 
    scale_size_continuous(range = c(-1,2)) + 
    scale_x_discrete("group") + 
    scale_y_discrete("gene") + 
    geom_segment(data = multi_genes_combination_ids,
                 aes(x = combination_id, xend = combination_id, 
                     y = min, yend = max, group = combination_id),
                 color = "black", linetype = "dashed") + # add lines => dashed ok? 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
  
  ### Plot gene prev / set size
  g4 <- ggplot(gene.prev, aes(x = fct_rev(genes), y = gene.prev)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    scale_y_reverse("Set size") +
    theme(#axis.text.y = element_blank(), # keep gene labels for now to make sure they align with other plot
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank())
  
  print(
    plot_spacer() + g1 + 
      plot_spacer() + g2 + 
      g4 + g3 + plot_layout(ncol = 2, widths = c(1, 4))
  )
  
}



