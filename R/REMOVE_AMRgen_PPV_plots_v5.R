#' Identify and visualize solo markers for antibiotic resistance
#'
#' This function identifies solo markers (genes uniquely found in one strain) 
#' for a given antibiotic and visualizes their association with phenotypic resistance.
#'
#' @param ast_matched Dataset containing antimicrobial susceptibility test (AST) results.
#' @param antibiotic Antibiotic of interest.
#' @param afp_matched Dataset with gene information for strains.
#' @param refgene_class Optional, specifies resistance gene class to filter by.
#' @param refgene_subclass Optional, specifies resistance gene subclass to filter by.
#' @param plot_cols Specifies colors for resistance phenotypes in the plots.
#' 
#' @import ggplot2 
#' @import dplyr 
#' @import forcats
#' @import patchwork
#' 
#' @return A list containing solo stats, solo count, and plots.
#' 
#' @examples
#' # Example usage (assuming datasets `ast_matched` and `afp_matched` are loaded)
#' result <- solo_marker_atb(ast_matched, "Ciprofloxacin", afp_matched, refgene_class = "Class A")
#' 
#' @export
#' 
#' Original version by Kat Holt, slight modifications by Laura Phillips and Elisa Sosa 
#' Outlook: - make an overall plot with all antibiotics, sort by isolates (to highlight interesting cases)
#'          - Filter for PPV >= 0.8 (or any other value) and add antibiotic class to the graph



solo_marker_atb <- function(ast_data, antibiotic, afp_data, microorganism_name, species_name = NULL, hq_filter_col = NULL, refgene_class = NULL, refgene_subclass = NULL, 
                            plot_cols = c("resistant" = "indianred", 
                                          "intermediate" = "#FFBF47", 
                                          "susceptible" = "lightgrey", 
                                          "non-susceptible (resistant & intermediate)" = "#FF7F50")) {
  
  # This function first matches identifiers across phenotype and genotype resistance data followed by identification and visualisation
  # of solo markers (genes) for a given antibiotic
  # Solo markers are genes uniquely found in one strain and are analyzed for their association 
  # with phenotypic resistance (e.g., resistant, susceptible).
  
  # Arguments:
  # - ast_data: Dataset containing antimicrobial susceptibility test (AST) results.
  # - antibiotic: Antibiotic of interest.
  # - afp_data: Dataset of genotype results.
  # - microorganism_name: Specify a specific microorganism of interest.
  # - species_name: Optional, specifies the species name to filter by. Default is NULL.
  # - hq_filter_col: A character string specifying the name of a column in `afp_data` that contains logical (TRUE/FALSE) values to filter high-quality rows. Default is NULL.
  # - refgene_class: Optional, specifies resistance gene class to filter by. Default is NULL.
  # - refgene_subclass: Optional, specifies resistance gene subclass to filter by. Default is NULL.
  # - plot_cols: Specifies colors for resistance phenotypes in the plots.
  
  # Step 1: Match and filter data for plotting
  
  #Strip any non-alphanumeric values from the start of all column names
    cleaned_names <- gsub("^[^a-zA-Z0-9]+", "", names(afp_data)) 
    names(afp_data) <- cleaned_names
  
    cleaned_names2 <- gsub("^[^a-zA-Z0-9]+", "", names(ast_data)) 
    names(ast_data) <- cleaned_names2
  
  # Optionally filter afp data by high-quality rows only.  
    if (!is.null(hq_filter_col)) { 
      afp_data <- afp_data %>% dplyr::filter(!!sym(hq_filter_col) == TRUE)
    }
    
  # Download microorganism relevant complexes from AMR package and join to afp data
    afp_data <- dplyr::left_join(
      afp_data, 
      AMR::microorganisms.groups %>% 
        dplyr::filter(grepl(microorganism_name, mo_name)) %>% 
        dplyr::rename(complex = mo_group_name, Species = mo_name) %>% 
        dplyr::select(complex, Species), 
      by = "Species" 
    )

  # Check for matched identifiers across afp and ast data
    matched_count <- sum(unique(ast_data$BioSample) %in% afp_data$Name)
    
  # Filter matched data to create dataframe of AST data with matched AFP data and vice versa
    ast_matched <- ast_data %>% 
      dplyr::filter(BioSample %in% afp_data$Name)
    
    afp_matched <- afp_data %>% 
      dplyr::filter(Name %in% ast_data$BioSample)
    
    # Check count of unique sample names
    ast_unique_count <- ast_matched %>% 
      dplyr::pull(BioSample) %>% 
      unique() %>% 
      length()
    
    afp_unique_count <- afp_matched %>% 
      dplyr::pull(Name) %>% 
      unique() %>% 
      length()
    
    print(paste0("There are ", ast_unique_count, " unique sample names identified in the genotype data and ", afp_unique_count, " unique sample names identified in the phenotype data"))
    
    # Conditionally filter for specific species if species_name is provided
    if (!is.null(species_name)) {
      afp_matched <- afp_matched %>% 
        dplyr::filter(Species == species_name)
    } else {
      afp_matched <- afp_matched
    }
    
    ast_matched <- ast_matched %>% 
      dplyr::filter(BioSample %in% afp_matched$Name)

    # Validate the antibiotic
    if (!antibiotic %in% unique(ast_matched$Antibiotic)) {
      stop(paste("Error: The antibiotic", antibiotic, "is not found in the dataset."))
    }
    
    # Step 2: Extract relevant resistance genes based on class or subclass
    if (!is.null(refgene_class) && !is.na(refgene_class)) {
      genes <- afp_matched %>% 
        dplyr::filter(Class %in% refgene_class) %>% 
        dplyr::pull(`Gene symbol`)
    } else if (!is.null(refgene_subclass) && !is.na(refgene_subclass)) {
      genes <- afp_matched %>% 
        dplyr::filter(stringr::str_detect(Subclass, refgene_subclass[1])) %>% 
        dplyr::pull(`Gene symbol`)
      if (length(refgene_subclass) > 1) {
        for (i in 2:length(refgene_subclass)) {
          genes <- c(genes, afp_matched %>% 
                       dplyr::filter(stringr::str_detect(Subclass, refgene_subclass[i])) %>% 
                       dplyr::pull(`Gene symbol`))
        }
      }
    } else {
      stop("Error: You must specify either refgene_class or refgene_subclass.")
    }
    
    if (length(genes) == 0) {
      stop("Error: No genes found for the specified refgene_class or refgene_subclass.")
    }
    
    # Step 3: Identify strains with only one relevant gene (solo strains)
    afp_solo_strains <- afp_matched %>% 
      dplyr::filter(`Gene symbol` %in% genes) %>% 
      dplyr::group_by(Name) %>% 
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% 
      dplyr::filter(n == 1) %>% 
      dplyr::pull(Name)
    
    if (length(afp_solo_strains) == 0) {
      stop(paste("Error: No solo strains found for the antibiotic", antibiotic, 
                 "and the specified gene class or subclass."))
    }
    
    # Step 4: Match AST results with solo strains and genes
    afp_solo_ast <- afp_matched %>% 
      dplyr::filter(`Gene symbol` %in% genes) %>% 
      dplyr::filter(Name %in% afp_solo_strains) %>% 
      dplyr::right_join(ast_matched %>% 
                          dplyr::filter(Antibiotic == antibiotic & !is.na(`Resistance phenotype`)), 
                        by = c("Name" = "BioSample"))
    
    # Step 5: Count the number of solo markers per gene
    solo_count <- afp_solo_ast %>% 
      dplyr::group_by(`Gene symbol`) %>% 
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    
    # Initialize variables for later plots
    solo_stats <- NULL
    solo_count_plot <- NULL
    combined_plot <- NULL
    
    # Step 6: Generate plots if there are solo markers
    if (nrow(solo_count) > 0) {
      solo_count_plot <- solo_count %>% 
        dplyr::filter(!is.na(`Gene symbol`)) %>%
        ggplot2::ggplot(ggplot2::aes(x = `Gene symbol`, y = n)) + 
        ggplot2::geom_col() + 
        ggplot2::coord_flip() + 
        ggplot2::theme_light()
      
      # PPV calculations
      solo_stats_R <- afp_solo_ast %>% 
        dplyr::group_by(`Gene symbol`) %>% 
        dplyr::summarise(total = sum(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible")), 
                         n = sum(`Resistance phenotype` == "resistant"), 
                         p = n / total, 
                         se = sqrt(p * (1 - p) / total), 
                         ci.lower = max(0, p - 1.96 * se), 
                         ci.upper = min(1, p + 1.96 * se), 
                         .groups = "drop") %>% 
        dplyr::mutate(category = "resistant")
      
      solo_stats_NWT <- afp_solo_ast %>% 
        dplyr::group_by(`Gene symbol`) %>% 
        dplyr::summarise(total = sum(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible")), 
                         n = sum(`Resistance phenotype` %in% c("resistant", "intermediate")), 
                         p = n / total, 
                         se = sqrt(p * (1 - p) / total), 
                         ci.lower = max(0, p - 1.96 * se), 
                         ci.upper = min(1, p + 1.96 * se), 
                         .groups = "drop") %>% 
        dplyr::mutate(category = "non-susceptible (resistant & intermediate)")
      
      solo_stats <- dplyr::bind_rows(solo_stats_R, solo_stats_NWT) %>% 
        dplyr::relocate(category, .before = n) %>%
        dplyr::mutate(Antibiotic = antibiotic) %>%
        dplyr::filter(!is.na(`Gene symbol`))
      
      pd <- ggplot2::position_dodge(width=0.8) # position dodge for coefficient plots
      
      ppv_plot <- solo_stats %>% 
        dplyr::filter(!is.na(`Gene symbol`)) %>%
        ggplot2::ggplot(ggplot2::aes(y = `Gene symbol`, group = category, col = category)) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = 2) + 
        ggplot2::geom_linerange(ggplot2::aes(xmin = ci.lower, xmax = ci.upper), position = pd) + 
        ggplot2::geom_point(ggplot2::aes(x = p), position = pd) + 
        ggplot2::theme_bw() +
        ggplot2::scale_y_discrete(position = "right") + 
        ggplot2::labs(y = "", x = "PPV", col = "Category") + 
        ggplot2::scale_colour_manual(values = plot_cols) + 
        ggplot2::xlim(0, 1)
      
      solo_pheno_plot <- afp_solo_ast %>% 
        dplyr::filter(!is.na(`Gene symbol`)) %>%
        dplyr::mutate(`Resistance phenotype` = forcats::fct_relevel(`Resistance phenotype`, "susceptible", "intermediate", "resistant")) %>%
        ggplot2::ggplot(ggplot2::aes(x = `Gene symbol`, fill = `Resistance phenotype`)) + 
        ggplot2::geom_bar(stat = "count", position = "fill") + 
        ggplot2::scale_fill_manual(values = plot_cols) + 
        ggplot2::coord_flip() + 
        ggplot2::theme_light() + 
        ggplot2::labs(x = "", y = "Proportion", fill = "Phenotype")
      
      combined_plot <- solo_pheno_plot + 
        ppv_plot + 
        patchwork::plot_layout(axes = "collect", guides = "collect") + 
        patchwork::plot_annotation(title = "Solo markers for resistance genes", subtitle = paste("vs", antibiotic, "phenotype"))
    }
    
    return(list(solo_stats = solo_stats, 
                solo_count = solo_count, 
                solo_count_plot = solo_count_plot, 
                combined_plot = combined_plot))
}


afp_file <-read_tsv("../../ATB_Enterobacter_AFP.tsv.gz") 
ast_file <-read_delim("../../AST_Enterobacter.tsv.gz") 


test <- solo_marker_atb(ast_file, antibiotic = "ciprofloxacin", afp_file, microorganism_name = "Enterobacter", species_name = "Enterobacter hormaechei_A", hq_filter_col = "HQ", refgene_class = "QUINOLONE")

test$solo_stats

test$combined_plot
  