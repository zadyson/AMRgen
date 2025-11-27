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

#' Generate a Series of Plots for AMR Gene and Combination Analysis
#'
#' This function generates a set of visualizations to analyze AMR gene combinations, MIC values, and gene prevalence from a given binary matrix. It creates several plots, including MIC distributions, a bar plot for the percentage of strains per combination, a dot plot for gene combinations in strains, and a plot for gene prevalence. It also outputs a table summarizing the MIC distribution (median, lower, upper) and number resistant, for each marker combination.
#' @param binary_matrix A data frame containing the original binary matrix output from the `get_binary_matrix` function. Expected columns are an identifier (column 1, any name), `pheno` (class sir, with S/I/R categories to colour points), `mic` (class mic, with MIC values to plot), and other columns representing gene presence/absence (binary coded, i.e., 1 = present, 0 = absent).
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots. Default is 2. Only genes with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#' - "" (default): decreasing frequency of combinations
#' - "genes": order by the number of genes in each combination
#' - "value": order by the median assay value (MIC or disk zone) for each combination.
#' @param plot_set_size Logical indicating whether to include a bar plot showing the set size (i.e., number of times each combination of markers is observed). Default is FALSE.
#' @param print_set_size Logical indicating whether, if `plot_set_size` is TRUE, to print the number of strains with each marker combination on the plot. Default is FALSE.
#' @param plot_category Logical indicating whether to include a stacked bar plot showing, for each marker combination, the proportion of samples with each phenotype classification (specified by the `pheno` column in the input file). Default is TRUE.
#' @param print_category_counts Logical indicating whether, if `plot_category` is TRUE, to print the number of strains in each resistance category for each marker combination in the plot. Default is FALSE.
#' @param boxplot_colour Colour for lines of the box plots summarising the MIC distribution for each marker combination. Default is "grey".
#' @param assay A character string indicating whether to plot MIC or disk diffusion data. Must be one of:
#' - "mic": plot MIC data stored in column `mic`
#' - "disk": plot disk diffusion data stored in column `disk`
#' @importFrom AMR as.mic scale_color_sir scale_fill_sir scale_y_mic
#' @importFrom dplyr any_of arrange desc distinct filter group_by if_else left_join mutate n pull relocate row_number select summarise ungroup
#' @importFrom forcats fct_rev
#' @importFrom ggplot2 aes after_stat coord_flip element_blank geom_bar geom_boxplot geom_col geom_point geom_segment geom_text ggplot labs position_fill scale_size_continuous scale_x_discrete scale_y_discrete scale_y_reverse theme theme_bw theme_light ylab
#' @importFrom patchwork plot_layout plot_spacer
#' @importFrom stats median
#' @importFrom tidyr pivot_longer unite
#' @return A list containing the following elements:
#' - `plot`: A grid of plots displaying: (i) grid showing the marker combinations observed, MIC distribution per marker combination, frequency per marker and (optionally) phenotype classification and/or number of samples for each marker combination.
#' - `summary`: A data frame summarizing each marker combination observed, including median MIC (and interquartile range), number of resistant isolates, and positive predictive value for resistance.
#' @export
#' @examples
#' \dontrun{
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#' binary_matrix <- get_binary_matrix(
#'   geno_table = ecoli_geno,
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno",
#'   keep_assay_values = TRUE,
#'   keep_assay_values_from = "mic"
#' )
#' amr_upset(binary_matrix, min_set_size = 3, order = "value", assay = "mic")
#' }
amr_upset <- function(binary_matrix, min_set_size = 2, order = "",
                      plot_set_size = FALSE, plot_category = TRUE,
                      print_category_counts = FALSE, print_set_size = FALSE,
                      boxplot_colour = "grey", assay = "mic") {
  
  # tidy up binary_matrix
 # col <- colnames(binary_matrix) # get column names

  # extract only the gene column names - need to exclude mic, disk, R, NWT (standard col names)
  # and the id column which will be the first col, doesn't matter what it's called
  # remaining columns will be the genes
 # cols_to_remove <- c("mic", "disk", "R", "NWT", "pheno", "ecoff")
 # genes <- col[-1]

  # gene names
#  genes <- setdiff(genes, cols_to_remove)
  
  if (sum(!is.na(binary_matrix$pheno))==0) {
    if (sum(!is.na(binary_matrix$ecoff))>0) {
      binary_matrix <- binary_matrix %>% mutate(pheno=ecoff)
      cat(" Warning: no values in pheno column, colouring upset plot by ecoff column\n")
    }
    else {stop(" Failed to make upset plot as no values in field pheno or ecoff")}
    colour_label = "ECOFF\ncategory"
  } else {colour_label = "Resistance\ncategory"}
  
  # gene names
  genes <- binary_matrix %>% select(-any_of(c("id", "pheno", "ecoff", "R", "NWT", "mic", "disk"))) %>% colnames()

  # check with have the expected assay column, with data
  if (!(assay %in% colnames(binary_matrix))) {
    stop(paste("input", deparse(substitute(binary_matrix)), "must have a column labelled ", assay))
  }
  data_rows <- binary_matrix %>%
    filter(!is.na(get(assay))) %>%
    nrow()
  if (data_rows == 0) {
    stop(paste("input", deparse(substitute(binary_matrix)), "has no non-NA values in column ", assay))
  }

  # Add in a combination column and filter to samples with the required assay data (MIC or disk) only
  binary_matrix_wide <- binary_matrix %>%
    filter(!is.na(get(assay))) %>%
    unite("combination_id", genes[1]:genes[length(genes)], remove = FALSE) # add in combinations

  # Make matrix longer
  binary_matrix <- binary_matrix_wide %>%
    pivot_longer(cols = genes[1]:genes[length(genes)], names_to = "genes") %>%
    mutate(genes = gsub("\\.\\.", ":", genes)) %>%
    mutate(genes = gsub("`", "", genes))

  ### Counts per combination, for bar plot - X axis = combination. Y axis = number of strains ###
  # This first to filter on combinations with enough data
  combination_freq <- binary_matrix_wide %>%
    group_by(combination_id) %>%
    summarise(n = n()) %>%
    mutate(perc = 100 * n / sum(n)) %>% # count number with each combination
    filter(n >= min_set_size) # count filter

  if (nrow(combination_freq) == 0) {
    stop(paste("No marker combinations pass the minimum frequency:", min_set_size))
  }
  # which have enough strains / data
  comb_enough_strains <- combination_freq %>% pull(combination_id)

  ### data for assay plot - dot plot. X axis = combination. Y axis = MIC or disk zone #####
  assay_plot <- binary_matrix %>%
    filter(combination_id %in% comb_enough_strains) %>%
    group_by(combination_id, get(assay), pheno) %>%
    summarise(n = n()) # count how many at each assay value, keep pheno for colour

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
  binary_matrix$point_size <- 2 * binary_matrix$value # want dot size to be larger than 1 => can make 2/3/4 etc

  # Only plot a line between points if more than one gene in a combination
  # get how many in each combination
  multi_genes_combination_id_all <- binary_matrix %>%
    select(combination_id, genes, value, point_size) %>%
    group_by(combination_id) %>%
    filter(value == 1) %>%
    mutate(u = length(unique(genes))) %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc, ordered = TRUE))

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
  if (order == "value") {
    if (assay == "mic") {
      # use mic class median (ie ignoring range indicators <>=) for the purpose of ordering columns
      mic_medians <- binary_matrix_wide %>%
        group_by(combination_id) %>%
        summarise(median = median(mic))
      ordered_comb_order <- mic_medians %>%
        arrange(median) %>%
        pull(combination_id)
    } else {
      disk_medians <- binary_matrix_wide %>%
        group_by(combination_id) %>%
        summarise(median = median(as.double(disk)))
      ordered_comb_order <- disk_medians %>%
        arrange(-median) %>%
        pull(combination_id)
    }
    assay_plot$combination_id <- factor(assay_plot$combination_id, levels = ordered_comb_order)
    combination_freq$combination_id <- factor(combination_freq$combination_id, levels = ordered_comb_order)
    binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  }

  ##### Plots ###

  ### assay data plot (MIC/disk)
  g1 <- ggplot(data = assay_plot, aes(x = combination_id, y = `get(assay)`)) +
    geom_boxplot(colour = boxplot_colour) +
    geom_point(aes(size = n, colour = pheno), show.legend = TRUE) +
    theme_bw() +
    scale_size_continuous("Number of\nisolates") +
    scale_color_sir(name = colour_label) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
    
  if (assay == "mic") {
    g1 <- g1 +
      scale_y_mic() +
      ylab("MIC (mg/L)")
  } else {
    g1 <- g1 +
      ylab("Disk zone (mm)")
  }

  ### Bar plot - set size
  g2 <- ggplot(combination_freq, aes(x = combination_id, y = n)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_y_reverse("Set size") +
    scale_x_discrete("group") +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  if (print_set_size) {
    g2 <- g2 + geom_text(aes(label = n), nudge_y = -.5)
  }

  # pheno category stacked barplot
  category_plot <- binary_matrix %>%
    select(id, pheno, combination_id) %>%
    distinct() %>%
    ggplot(aes(x = combination_id, fill = pheno)) +
    geom_bar(stat = "count", position = "fill") +
    scale_fill_sir() +
    theme_light() +
    labs(x = "", y = "Category") +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )
  if (print_category_counts) {
    category_plot <- category_plot +
      geom_text(aes(label = after_stat(count)), stat = "count", position = position_fill(vjust = .5), size = 3)
  }

  ### Dot plot of combinations

  g3 <- binary_matrix <- binary_matrix %>%
    mutate(binary_comb = if_else(value > 0, 1, 0)) %>%
    ggplot(aes(x = combination_id, y = fct_rev(genes))) +
    geom_point(aes(size = binary_comb), show.legend = FALSE) +
    theme_bw() +
    scale_size_continuous(range = c(-1, 2)) +
    scale_y_discrete(name = "Marker") +
    geom_segment(
      data = multi_genes_combination_ids,
      aes(
        x = combination_id, xend = combination_id,
        y = min, yend = max, group = combination_id
      ),
      color = "black"
    ) + # add lines
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  ### Plot gene prev / set size
  g4 <- gene.prev %>%
    ggplot(aes(x = fct_rev(genes), y = gene.prev)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    scale_y_reverse("") +
    theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  # assemble plot
  final_plot <- plot_spacer() + g1 + plot_layout(ncol = 2, widths = c(1, 4), guides = "collect")
  if (plot_category) {
    final_plot <- final_plot + plot_spacer() + category_plot
  }
  final_plot <- final_plot + g4 + g3
  if (plot_set_size) {
    final_plot <- final_plot + plot_spacer() + g2
  }

  # set relative plotting heights
  if (plot_category & plot_set_size) {
    final_plot <- final_plot + plot_layout(heights = c(2, 1, 2, 1))
  }
  if (plot_category & !plot_set_size) {
    final_plot <- final_plot + plot_layout(heights = c(2, 1, 2))
  }
  if (!plot_category & plot_set_size) {
    final_plot <- final_plot + plot_layout(heights = c(2, 2, 1))
  }

  print(final_plot)

  # summary table (ignore MIC values expressed as ranges, when calculating median/IQR)
  if (assay == "mic") {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(
        median = median(as.double(as.mic(mic, keep_operators = FALSE)), na.rm = TRUE),
        q25 = stats::quantile(as.double(as.mic(mic, keep_operators = FALSE)), 0.25, na.rm = TRUE),
        q75 = stats::quantile(as.double(as.mic(mic, keep_operators = FALSE)), 0.75, na.rm = TRUE),
        ppv = mean(R, na.rm = TRUE),
        R = sum(R, na.rm = TRUE),
        n = n()
      )
  } else {
    summary <- binary_matrix_wide %>%
      group_by(combination_id) %>%
      summarise(
        median = median(as.double(as.disk(disk)), na.rm = TRUE),
        q25 = stats::quantile(as.double(as.disk(disk)), 0.25, na.rm = TRUE),
        q75 = stats::quantile(as.double(as.disk(disk)), 0.75, na.rm = TRUE),
        ppv = mean(R, na.rm = TRUE),
        R = sum(R, na.rm = TRUE),
        n = n()
      )
  }

  # get names for summary
  combination_names <- binary_matrix_wide %>%
    select(combination_id, any_of(genes)) %>%
    distinct()

  combination_names <- combination_names %>%
    mutate(marker_list = apply(., 1, function(row) {
      paste(names(combination_names)[-1][row[-1] == 1], collapse = ", ")
    }), .after = combination_id) %>%
    mutate(marker_count = rowSums(. == 1)) %>%
    select(combination_id, marker_list, marker_count)

  summary <- summary %>%
    left_join(combination_names, by = "combination_id") %>%
    select(-combination_id) %>%
    mutate(marker_list = if_else(is.na(marker_list), "-", marker_list)) %>%
    mutate(marker_count = if_else(is.na(marker_count), 0, marker_count)) %>%
    relocate(marker_list, .before = median) %>%
    relocate(marker_count, .before = median) %>%
    mutate(marker_list = gsub("\\.\\.", ":", marker_list)) %>%
    mutate(marker_list = gsub("`", "", marker_list))

  return(list(plot = final_plot, summary = summary))
}
