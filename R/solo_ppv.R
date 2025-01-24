require(AMR)
require(patchwork)

#' Perform Solo PPV Analysis for AMR Markers
#'
#' This function performs a Positive Predictive Value (PPV) analysis for AMR markers
#' associated with a given antibiotic and drug class. The function calculates the PPV
#' for solo markers (those with only one genetic marker relevant to the drug class) 
#' and visualizes the results using various plots. It returns a list containing summary 
#' statistics for each solo marker, and associated plots showing the breakdown of 
#' resistance phenotypes, and PPV (with 95% confidence interval) for each marker.
#'
#' @param geno_table A data frame containing genotype data, including at least one column labeled
#'   `drug_class` for drug class information and one column for sample identifiers (specified via
#'   `geno_sample_col` otherwise it is assumed the first column contains identifiers).
#' 
#' @param pheno_table A data frame containing phenotype data, which must include a column `drug_agent`
#'   (with the antibiotic information) and a column with the resistance interpretation (S/I/R, colname 
#'   specified via `sir_col`).
#' 
#' @param antibiotic A character string specifying the antibiotic of interest to filter phenotype data.
#'   The value must match one of the entries in the `drug_agent` column of `pheno_table`.
#' 
#' @param drug_class_list A character vector of drug classes to filter genotype data for markers 
#'   related to the specified antibiotic. Markers in `geno_table` will be filtered based on whether
#'   their `drug_class` matches any value in this list.
#' 
#' @param geno_sample_col A character string (optional) specifying the column name in `geno_table`
#'   containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column 
#'   contains identifiers.
#' 
#' @param pheno_sample_col A character string (optional) specifying the column name in `pheno_table`
#'   containing sample identifiers. Defaults to `NULL`, in which case it is assumed the first column 
#'   contains identifiers.
#' 
#' @param sir_col A character string specifying the column name in `pheno_table` that contains
#'   the resistance interpretation (SIR) data. The values should be interpretable as "R" (resistant), 
#'   "I" (intermediate), or "S" (susceptible).
#'   
#' @param plot_cols A named vector of colors for the plot. The names should be the phenotype categories 
#'   (e.g., "R", "I", "S", "NWT"), and the values should be valid color names or hexadecimal color codes.
#'   Default colors are provided for resistant ("R"), intermediate ("I"), susceptible ("S"), and non-wild-type ("NWT").
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{solo_stats}{A dataframe summarizing the PPV for resistance (R vs S/I) and NWT (R/I vs S), 
#'     including the number of positive hits, sample size, PPV, and 95% confidence intervals for each marker.}
#'     \item{combined_plot}{A combined ggplot object showing the PPV plot for the solo markers, and a bar plot 
#'     for the phenotype distribution.}
#'     \item{solo_binary}{A dataframe with binary values indicating the presence or absence of the solo markers.}
#'     \item{amr_binary}{A dataframe with binary values for the AMR markers, based on the input genotype and phenotype data.}
#'   }
#'
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot aes geom_vline geom_linerange geom_point theme_bw labs scale_y_discrete
#'   scale_colour_manual coord_flip geom_bar geom_text theme_light position_fill
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map
#' @importFrom patchwork plot_layout plot_annotation
#'
#' @examples
#' geno_table <- parse_amrfp("testdata/Ecoli_AMRfinderplus.tsv", "Name")
#' pheno_table <- import_ncbi_ast("testdata/Ecoli_AST_NCBI.tsv")
#' 
#' soloPPV_cipro <- solo_ppv_analysis(geno_table = geno_table, 
#'                             pheno_table = pheno_table, 
#'                             antibiotic = "Ciprofloxacin", 
#'                             drug_class_list = c("Quinolones"),
#'                             sir_col="Resistance phenotype")
#'                             
#' soloPPV_mero <- solo_ppv_analysis(geno_table = geno_table, 
#'                             pheno_table = pheno_table, 
#'                             antibiotic="Meropenem",
#'                             drug_class_list=c("Carbapenems", "Cephalosporins"),
#'                             sir_col="Resistance phenotype")
#' 
#' # View the results
#' soloPPV_cipro$solo_stats
#' soloPPV_cipro$combined_plot
#' soloPPV_mero$solo_stats
#' soloPPV_mero$combined_plot
#'
#' @export
solo_ppv_analysis <- function(geno_table, pheno_table, antibiotic, drug_class_list,
                              geno_sample_col=NULL, pheno_sample_col=NULL, sir_col=NULL,
                              plot_cols=c("R"="IndianRed", "I"="orange", "S"="lightgrey",
                                          "NWT"="navy")) {
  
  # check there is a SIR column specified
  if (is.null(sir_col)) {stop("Please specify a column with S/I/R values, via the sir_col parameter.")}
  if (!(sir_col %in% colnames(pheno_table))) {
    stop(paste0("Column '", sir_col,"' not found in input phenotype data. Please specify a valid column with S/I/R values, via the sir_col parameter."))
  }
  
  # get binary matrix
  amr_binary <- getBinMat(geno_table, pheno_table, antibiotic=antibiotic, 
                          drug_class_list=drug_class_list, 
                          geno_sample_col=geno_sample_col, pheno_sample_col=pheno_sample_col, 
                          sir_col=sir_col)
  
  # get solo markers
  marker_counts <- amr_binary %>% select(-id, -pheno, -R, -NWT) %>% rowSums()
  solo_binary <- amr_binary %>% mutate(solo=marker_counts==1) %>% 
    filter(solo) %>% 
    pivot_longer(!c(id, pheno, R, NWT, solo)) %>% 
    filter(value==1) %>% 
    filter(!is.na(pheno))
  
  if (nrow(solo_binary)==0) { stop("No solo markers found") }
  
  # summarise numerator, denominator, proportion, 95% CI - for R and NWT
  solo_stats_R <- solo_binary %>% 
    group_by(name) %>% 
    summarise(x=sum(R, na.rm=T), 
              n=n(), 
              p=x/n,
              se=sqrt(p*(1-p)/n), 
              ci.lower=max(0,p-1.96*se), 
              ci.upper=min(1,p+1.96*se)) %>%
    mutate(category="R")

  solo_stats_NWT <- solo_binary %>% 
    group_by(name) %>% 
    summarise(x=sum(NWT), 
              n=n(), 
              p=x/n,
              se=sqrt(p*(1-p)/n), 
              ci.lower=max(0,p-1.96*se), 
              ci.upper=min(1,p+1.96*se)) %>%
    mutate(category="NWT")
  
  solo_stats <- bind_rows(solo_stats_R, solo_stats_NWT) %>% 
    relocate(category, .before=x)
  
  # plots
  
  pd <- position_dodge(width=0.8) # position dodge for coefficient plots
  
  ppv_plot <- solo_stats %>% 
    ggplot(aes(y=name, group=category, col=category)) +
    geom_vline(xintercept=0.5, linetype=2) +
    geom_linerange(aes(xmin=ci.lower, xmax=ci.upper), position=pd) +
    geom_point(aes(x=p), position=pd) + 
    theme_bw() +
    scale_y_discrete(labels=paste0("(n=",solo_stats$n,")"), position="right") + 
    labs(y="", x="PPV", col="Category") + 
    scale_colour_manual(values=plot_cols) + xlim(0,1)
  
  solo_pheno_plot <- solo_binary %>%
    ggplot(aes(x=name, fill=pheno)) + 
    geom_bar(stat="count", position="fill") + 
    scale_fill_manual(values=plot_cols) + 
    coord_flip() + 
    geom_text(aes(label=..count..), stat="count", position=position_fill(vjust = .5), size=3) + 
    scale_x_discrete(solo_stats$name) +
    theme_light() + labs(x="", y="Proportion", fill="Phenotype")
  
  header <- paste("Solo markers for class:", paste0(drug_class_list, collapse=", ")) 
  
  combined_plot <- solo_pheno_plot + 
    ppv_plot + 
    patchwork::plot_layout(axes="collect", guides="collect") + 
    patchwork::plot_annotation(title=header, subtitle=paste("vs phenotype for drug:", antibiotic))
  
  return(list(solo_stats=solo_stats, combined_plot=combined_plot, solo_binary=solo_binary, amr_binary=amr_binary))
  
}