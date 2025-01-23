#' Generate a Stacked Bar Plot of MIC Values Colored by Gene Symbol for Each Antibiotic
#'
#' This function creates a stacked bar plot using `ggplot2`, where the x-axis represents MIC (Minimum Inhibitory Concentration) values, 
#' the y-axis represents their frequency, and the bars are colored by the gene symbols. The plot is faceted by antibiotic. 
#' @param pheno_data A data.frame containing columns "BioSample", "Antibiotic", 
#' "Measurement sign", "MIC (mg/L)", "Disk diffusion (mm)","Resistance phenotype".
#' BioSample is the ID used to match the entries in geno_data. 
#' MIC (Minimum Inhibitory Concentration) and Disk diffusion are numeric values, the other character strings.
#' @param geno_data A data.frame containing columns "Name", "Gene symbol", 
#' "Class" and "Subclass"  
#'
#' @param pathogen_mo A string in mo format (?????) specifying what pathogen 
#' is used to select the breakpoints, e.g. 'B_STPHY_AURS' 
#' or group the data by different pathogens.
#' 
#' @param measure_method A string specifying the method used to measure antimicrobial resistance.
#' The default is 'mic' for Minimum Inhibitory Concentration, the other option is 'disk' for Disk Diffusion'.
#' 
#' @param breakpoint_source A string representing the source for breakpoint data. By default, it is set 
#' to the EUCAST standard from two years ago.
#' 
#' @param breakpoint_type A string specifying the type of breakpoint to use, default is 'human'.
#' 
#' @param color_by A string specifying the variable used to color the plot elements. Options include 
#' 'gene' (default) to color by gene, or 'class' to color by gene class. Any other strings also default to gene.
#' 
#' @param plot_title A string specifying the title of the plot. The default title is 
#' "Frequency distribution of MIC Value for different genotypes".
#' 
#' @param leg_pos A string specifying the position of the plot's legend. The default is 'bottom', 
#' but other options include 'top', 'left', or 'right'.
#' 
#' @param ncol_leg An integer specifying the number of columns for the legend. By default, it is 
#' calculated based on the number of unique gene symbols in `geno_data`, rounded to 10 items per column, but minimally 2.
#' This ensures the legend is well-organized when there are many gene symbols.
#' 
#' @return A ggplot2 object representing the stacked bar plot.
#' 
#' @details
#' - The function automatically groups the data by `Antibiotic`, `MIC (mg/L)`, and `Gene symbol` and calculates the frequency of occurrences.
#' - The plot uses the `viridis` color scale with the "turbo" palette by default. You can replace it with your custom palette if desired.
#'
#' @examples
#' AMRgen_barplot_mic_gen(pheno_data %>%
#'                          filter(Antibiotic %in% c('tetracycline','ceftaroline', 'levofloxacin')), 
#'                        geno_data, 
#'                        measure_method = 'MIC',
#'                        pathogen_mo = 'B_STPHY_AURS')
#' 
#' AMRgen_barplot_mic_gen(pheno_data %>%
#'                          filter(Antibiotic %in% c('erythromycin','penicillin', 'levofloxacin')), 
#'                        geno_data, 
#'                        color_by = 'class',
#'                        measure_method = 'DISK',
#'                        pathogen_mo = 'B_STPHY_AURS')
#'
#' @import ggplot2
#' @import dplyr
#' @import viridis
#' @import AMR
#' @export
#' 
AMRgen_barplot_mic_gen <- function(pheno_data, geno_data,
                                   sample_col_name = NULL, pheno_sample_col = NULL, geno_sample_col = NULL, # default is null and assume first col in each, otherwise user 
                                   # provides specifics or a single name if it's the same in both
                                   pathogen_mo,
                                   abs_to_plot, # list of antibiotics you want to plot
                                   
                                   measure_method = 'mic',
                                   breakpoint_source = paste('EUCAST', as.numeric(format(Sys.Date(), "%Y"))-2),
                                   breakpoint_type = 'ECOFF',
                                   color_by = 'gene', #alt.: 'class'
                                   
                                   plot_title = "Frequency distribution of MIC Value for different genotypes",
                                   leg_pos = 'bottom',
                                   ncol_leg = max(2, round(length(unique(geno_data$marker))/10))
                                   ){
  # todo: 
  # what about measurement signs? (ignoring for now)
  # add a lot more parameters (ggplot2)

  #TODO (JH):
  # - fix sample ID column matching (provide as argument to function, or assume the first column and use that name)

  # - fix name of columns in pheno (should be lowercase mic or disk)
  # - references to "Antibiotic" should now be "drug_agent" as this is the new name of this col, which is automatically ab class

  source("helpers.R")

  # 1. get breakpoint data ----
  if (toupper(breakpoint_source) %in% clinical_breakpoints$guideline) {
    breakpoints <- toupper(breakpoint_source) # Replace with actual path
    message("Using ", breakpoint_source, " breakpoints.")
  } else {
    message("Invalid breakpoint source. Using default: " ,
         paste('EUCAST', as.numeric(format(Sys.Date(), "%Y"))-2))
    breakpoints <- paste('EUCAST', as.numeric(format(Sys.Date(), "%Y"))-2)
  }
  
  if (breakpoint_type %in% clinical_breakpoints$type) {
    breakpoint_type <- breakpoint_type # Replace with actual path
  } else {
    message("Invalid breakpoint type. Using default: human")
    breakpoint_type <- 'human'
  }
  
  if(!measure_method %in% c('mic', 'disk')){
    stop('Please choose one of the methods "mic" or "disk" for method.')
  }

  # get unique antibiotics in pheno file
  ab_bp <- data.frame('ab' = unique(pheno_data$drug_agent))

  # select only those we want to plot
  ab_bp <- ab_bp %>% filter(ab %in% abs_to_plot)
           
  clinical_breakpoints_here <- clinical_breakpoints %>%
    filter(type == breakpoint_type,
           guideline == breakpoint_source,
           method == measure_method,
           mo == pathogen_mo,
           ab %in% ab_bp$ab)
  
  ab_bp <- ab_bp %>%
    left_join(clinical_breakpoints_here %>%
                select(ab, breakpoint_S, breakpoint_R, disk_dose),
              by = 'ab')
  
  ## NEED TO RETURN A WARNING HERE IF THERE ARE NO GUIDELINES FOR THE CHOSEN DRUGS
  
  # 2. prep data ----
  # get the column names to merge on - if provided use those, otherwise use first col as default
  sample_id_results <- compare_geno_pheno_id(geno_data, pheno_data, geno_sample_col, pheno_sample_col)

  # filter data to get only samples in both
  pheno_data_overlap <- pheno_data %>% filter(!!sym(pheno_sample_col) %in% sample_id_results$overlap_ids)
  geno_data_overlap <- geno_data %>% filter(!!sym(geno_sample_col) %in% sample_id_results$overlap_ids)

  # we want to select only pheno data and geno data matching our chosen Abs
  # HOWEVER - I think we want to also include the drug_class in the pheno data, because otherwise we're going to miss entries where we have no agent, but we do have the class
  ab_groups_plot <- antibiotics %>% filter(ab %in% abs_to_plot) %>% pull(group)

  pheno_data_overlap <- pheno_data_overlap %>% filter(ab %in% abs_to_plot)
  geno_data_overlap <- geno_data_overlap %>% filter(drug_agent %in% abs_to_plot | drug_class %in% ab_groups_plot)

  # now we want to colour bars by combos of genes in each sample, so first need to
  # list the full set of markers for a drug in each genome
  geno_data_markers <- geno_data_overlap %>% 
    group_by(Name, drug_class) %>% 
    summarise(markers = paste(marker, collapse = ", "), .groups="drop")

  # merge with the pheno data (we can now as we have one row per sample!!)
  # only merge on rows where the drug class matches
  merged_data <- pheno_data_overlap %>%
    mutate(drug_class = ab_group(drug_agent), .after=drug_agent) %>%
    left_join(geno_data_markers, 
      by = c(setNames(geno_sample_col, pheno_sample_col), "drug_class"))
  
  # now we want to plot the frequency of a sample that is resistant to a drug
  # note that geom_bar will autocount
  ggplot(pheno_data_overlap, aes(x=mic)) + 
    geom_bar() + 
    facet_wrap(~drug_agent)
  
  # todo add disk breakpoint
  mic_summary <- merged_data %>%
    left_join(ab_bp, join_by("drug_agent" == "ab")) %>%
    mutate(fill_col = case_when(
      color_by == 'class' ~ ifelse(is.na(drug_class), 'NA', drug_class), 
      TRUE ~ markers))
  
  #3. plot ----
  ggplot(mic_summary, aes_string(x = measure_method, fill = "fill_col")) +
    geom_bar() +  # Use stack for stacking bars
    facet_wrap(~ drug_agent, scales = "free", ncol = NULL) +  # Create a panel for each antibiotic
    theme_minimal() +
    #exchange to our own color palette?
    #scale_fill_viridis_d(option = "turbo") +  # Use "turbo" or other options like "viridis", "plasma"
    labs(
      title = plot_title,
      x = "MIC (mg/L)",
      y = "Frequency",
      fill = "Gene Symbol"
    ) +
    theme(
      legend.position = leg_pos,  # Move the legend to the top
      legend.title = element_blank(),
      legend.text = element_text(size = 8),  # Reduce legend text size
      legend.key.size = unit(0.4, "cm"),  # Reduce the size of legend keys
      axis.text.x = element_text(hjust = 1),
      strip.text = element_text(size = 10, face = "bold") # Customize panel titles
    ) +
    guides(fill = guide_legend(ncol = ncol_leg)) + # Arrange legend items into 4 columns
    geom_vline(data = ab_bp, aes(xintercept = `breakpoint_S`), 
               color = "black", linetype = "dashed") +  # Add breakpoint_S lines
    geom_vline(data = ab_bp, aes(xintercept = `breakpoint_R`), 
               color = "black", linetype = "dashed")   # Add breakpoint_R lines
}

# questions
# can we really expect the input data to look like this?
# should the user be able to paas all the parameters you theoretically could to ggplot2?


# example ----
library(tidyverse)
library(ggplot2)
library(AMR)
library(viridis) #might replace by own colors
setwd('S:/OE/FG37/BotheC/AMRHackathon')
# what is add reference ditribution note? 
# each row is drug-bur combi i.e. one tested antibiotic for one sample
#afpSA:
# focus on name and gene symbol
pheno_data <- read_tsv(fs::path(getwd(), 'Data visualisation', 'astSA_matched.tsv'))
# columns: id, gene, ab it encodes res for and corresponding class
geno_data <- read_tsv(fs::path(getwd(), 'Data visualisation', 'afpSA_matched.tsv'))
AMRgen_barplot_mic_gen(pheno_data %>%
                         filter(Antibiotic %in% c('tetracycline','ceftaroline', 'levofloxacin')), 
                       geno_data, 
                       measure_method = 'MIC',
                       pathogen_mo = 'B_STPHY_AURS')

AMRgen_barplot_mic_gen(pheno_data %>%
                         filter(Antibiotic %in% c('erythromycin','penicillin', 'levofloxacin')), 
                       geno_data, 
                       color_by = 'class',
                       measure_method = 'DISK',
                       pathogen_mo = 'B_STPHY_AURS')

