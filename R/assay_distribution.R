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
#
#' Generate a Stacked Bar Plot of Assay Values Colored by a Variable
#'
#' This function creates a stacked bar plot using `ggplot2`, where the x-axis represents MIC (Minimum Inhibitory Concentration) or disk values, the y-axis indicates their frequency, and the bars are colored by a variable (by default, colours indicate whether the assay value is expressed as a range or not). Plots can optionally be faceted on an additional categorical variable. If breakpoints are provided, or species and drug are provided so we can extract EUCAST breakpoints, vertical lines indicating the S/R breakpoints and ECOFF will be added to the plot. 
#' @param pheno_table Phenotype table in standard format as per import_ast().
#' @param antibiotic Name of the antibiotic.
#' @param measure Field name containing the assay measurements to plot (default "mic").
#' @param facet_var (optional) Field name containing a field to facet on (default NULL).
#' @param colour_by (optional) Field name containing a field to colour bars by (default NULL, which will colour each bar to indicate whether the value is expressed as a range or not)
#' @param species (optional) Name of species, so we can retrieve breakpoints to print at the top of the plot to help interpret it.
#' @param bp_site (optional) Breakpoint site to retrieve (only relevant if also supplying `species` to retrieve breakpoints, and not supplying breakpoints via `bp_S`, `bp_R`, `ecoff`).
#' @param species (optional) Name of species, so we can retrieve breakpoints to print at the top of the plot to help interpret it.
#' @param bp_S (optional) S breakpoint
#' @param bp_R (optional) R breakpoint
#' @param ecoff (optional) ECOFF breakpoint
#' @param marker_free_strains (optional) Vector of sample names to select to get their own plot. Most useful for defining the set of strains with no known markers associated with the given antibiotic, so you can view the distribution of assay values for strains expected to be wildtype, which can help to identify issues with the assay.
#' @param cols (optional) Manual colour scale to use for plot. If NULL, `colour_by` variable is of class 'sir', bars will by default be coloured using standard SIR colours.
#' @importFrom ggplot2 aes element_text facet_wrap geom_bar geom_vline ggplot labs theme scale_fill_manual sym
#' @return A list containing
#' \item{plot}{Main plot with all samples that have assay data for the given antibiotic}
#' \item{plot_nomarkers}{Additional plot showing only those samples listed in `marker_free_strains`}
#' @export
assay_by_var <- function(pheno_table, antibiotic, measure="mic", facet_var=NULL, 
                         species=NULL, bp_site=NULL, bp_S=NULL, bp_R=NULL, ecoff=NULL,
                         marker_free_strains=NULL, colour_by=NULL, cols=NULL) {
  if (measure %in% colnames(pheno_table)) {
    pheno_table <- pheno_table %>%
      filter(!is.na(get(measure))) %>%
      filter(drug_agent==as.ab(antibiotic)) %>%
      arrange(get(measure))
  } else {stop(paste0("No '", measure, "' column in input table"))}
  if (!is.null(facet_var)){
    if (!(facet_var %in% colnames(pheno_table))) {
      stop(paste0("Facet variable '", facet_var, "' not found in input table"))
    }
  }
  if (!is.null(colour_by)) {
    if(!(colour_by %in% colnames(pheno_table))) {
      cat(paste0("WARNING: Colour variable '", colour_by, "' not found in input table\n"))
      colour_by <- NULL
    }
  }
  
  if(is.null(colour_by)) {
    pheno_table <- pheno_table %>%
      mutate(range=if_else(grepl("<",get(measure)), "range", "value"))
    colour_by <- "range"
    cols <- c(range="maroon", value="navy", `NA`="grey")
  }
  
  # if species provided, but breakpoints aren't, check breakpoints to annotate plot
  if (!is.null(species) & is.null(bp_S) & is.null(bp_R) & is.null(ecoff)) {
    if (measure %in% c("mic", "disk")) {
      ecoff <- safe_execute(getBreakpoints(species=species, guide="EUCAST 2025", antibiotic=antibiotic, "ECOFF") %>% filter(method==toupper(measure)) %>% pull(breakpoint_S))
      bp_S <- safe_execute(unlist(checkBreakpoints(species=species, guide="EUCAST 2025", antibiotic=antibiotic, bp_site=bp_site, assay=toupper(measure))[1]))
      bp_R <- safe_execute(unlist(checkBreakpoints(species=species, guide="EUCAST 2025", antibiotic=antibiotic, bp_site=bp_site, assay=toupper(measure))[2]))
    } 
  } 
  if (!is.null(bp_S) | !is.null(bp_R) | !is.null(ecoff)) {
    if (measure=="mic") {subtitle=paste("ECOFF:", ecoff, "S <=", bp_S, "R>", bp_R)}
    else if (measure=="disk") {subtitle=paste("ECOFF:", ecoff, "S >=", bp_S, "R<", bp_R)}
  } else {subtitle=NULL}
  
  pheno_data <- pheno_table %>% count(factor(!!sym(measure)), !!sym(colour_by))
  colnames(pheno_data)[1] <- measure
  
  # get x coordinates for breakpoints after converting to factor
  if (!is.null(bp_S) | !is.null(bp_R) | !is.null(ecoff)) {
    assay_order <- pheno_table %>% count(factor(!!sym(measure)))
    colnames(assay_order)[1] <- measure
    bp_S <- c(1:nrow(assay_order))[assay_order[,1]==bp_S]
    bp_R <- c(1:nrow(assay_order))[assay_order[,1]==bp_R]
    ecoff <- c(1:nrow(assay_order))[assay_order[,1]==ecoff]
  } 
  
  # plot distribution per variable
  if (nrow(pheno_table)>0) {
    plot_all <- pheno_table %>%
      ggplot(aes(x=factor(!!sym(measure)), fill=!!sym(colour_by))) +
      geom_bar() +
      labs(x="Measurement", y="count", fill="Value", subtitle=subtitle,
           title=paste(antibiotic, "assay distributions for all samples")) +
      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
    if (!is.null(cols)) {
      plot_all <- plot_all + scale_fill_manual(values=cols)
    }
    if (!is.null(facet_var)) { 
      if (pheno_table %>% filter(!is.na(get(facet_var))) %>% nrow() > 0) {
        plot_all <- plot_all + facet_wrap(~get(facet_var), ncol=1, scales="free_y")
      }
    }
    # add breakpoints
    if (!is.null(bp_S)) {
      plot_all <- plot_all + geom_vline(xintercept=bp_S, color="#3CAEA3", linetype=2) }
    if (!is.null(bp_R)) {
      plot_all <- plot_all + geom_vline(xintercept=bp_R, color="#ED553B", linetype=2) }
    if (!is.null(ecoff)) {
      plot_all <- plot_all + geom_vline(xintercept=ecoff, color="grey", linetype=2) }
  } else {plot_all <- NULL}
  
  # samples with no markers
  if (!is.null(marker_free_strains)) {
    pheno_table_nomarkers <- pheno_table %>% filter(id %in% marker_free_strains)
    if (nrow(pheno_table_nomarkers)>0) {
      plot_nomarkers <- pheno_table_nomarkers %>%
        ggplot(aes(x=factor(!!sym(measure)), fill=!!sym(colour_by))) +
        geom_bar() +
        labs(x="Measurement", y="count", fill="Value", subtitle=subtitle,
             title=paste(antibiotic, "assay distributions for samples with no markers identified")) +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
      if (!is.null(cols)) {
        plot_nomarkers <- plot_nomarkers + scale_fill_manual(values=cols)
      }
      if (!is.null(facet_var)) { 
        if (pheno_table_nomarkers %>% filter(!is.na(get(facet_var))) %>% nrow()>0) {
          plot_nomarkers <- plot_nomarkers + facet_wrap(~get(facet_var), ncol=1, scales="free_y")
        }
      }
      # add breakpoints
      if (!is.null(bp_S)) {
        plot_nomarkers <- plot_nomarkers + geom_vline(xintercept=bp_S, color="#3CAEA3", linetype=2) }
      if (!is.null(bp_R)) {
        plot_nomarkers <- plot_nomarkers + geom_vline(xintercept=bp_R, color="#ED553B", linetype=2) }
      if (!is.null(ecoff)) {
        plot_nomarkers <- plot_nomarkers + geom_vline(xintercept=ecoff, color="grey", linetype=2) }
    } else {plot_nomarkers <- NULL}
  } else {plot_nomarkers <- NULL}
  
  return(list(plot_nomarkers=plot_nomarkers, plot=plot_all))
}

