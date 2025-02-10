#' Plot Estimates from a Table of Results
#'
#' This function creates a ggplot object visualizing logistic regression coefficients
#' with their 95% confidence intervals. Significant markers are highlighted based
#' on a specified p-value threshold.
#'
#' @param tbl A data frame or tibble containing the logistic regression results.
#'   Expected columns are:
#'   - `marker`: The name of the marker (e.g., variable name).
#'   - `pval`: The p-value for each marker.
#'   - `ci.lower`: The lower bound of the confidence interval.
#'   - `ci.upper`: The upper bound of the confidence interval.
#'   - `est`: The estimated coefficient.
#' @param sig (optional) The significance threshold for p-values. Defaults to 0.05.
#' @param sig_colors (optional) A vector of two colors to represent significant and non-significant estimates.
#' @param x_title (optional) The title for the x-axis. Defaults to "Coefficient (95% CI)".
#' @param y_title (optional) The title for the y-axis. Defaults to "Variant".
#' @param title (optional) The main title of the plot. If NULL, no title is added.
#' @param axis_label_size (optional) The font size of the axis labels. Defaults to 9.
#'
#' @return A ggplot object showing the logistic regression coefficients with confidence
#'   intervals. Significant markers (p-value < \code{sig}) are colored differently.
#'
#' @examples
#' # Example dataset
#' tbl <- tibble::tibble(
#'   marker = c("(Intercept)", "var1", "var2", "var3"),
#'   pval = c(0.1, 0.03, 0.2, 0.04),
#'   ci.lower = c(-0.2, 0.1, -0.3, 0.2),
#'   ci.upper = c(0.5, 0.8, 0.4, 1.1),
#'   est = c(0.2, 0.5, 0.1, 0.7)
#' )
#'
#' # Plot
#' plot_estimates(tbl)
#'
#' @export
plot_estimates <- function(tbl, sig = 0.05, 
                           sig_colors=c("grey", "blue4"),
                           x_title="Coefficient (95% CI)",
                           y_title="Variant", 
                           title=NULL,
                           axis_label_size=9) {
  plot <- tbl %>%
    mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>%
    filter(marker != "(Intercept)") %>%
    mutate(marker = gsub("`", "", marker)) %>%
    ggplot(aes(y = marker, col = sig_binary)) +
    geom_vline(xintercept = 0) +
    geom_linerange(aes(xmin = ci.lower, xmax = ci.upper)) +
    geom_point(aes(x = est)) +
    scale_color_manual(values = sig_colors) +
    theme_light() +
    theme(axis.text.x=element_text(size=axis_label_size),
          axis.text.y=element_text(size=axis_label_size)) +
    labs(title = title,
         x = x_title,
         y = y_title,
         col = paste0("p<", sig))

  return(plot)
}

#' Plot to Compare Two Sets of Estimates
#' 
#' This function compares two sets of estimates by creating a plot that overlays the estimates and confidence intervals for both sets. It can also display the estimates in two separate plots.
#' 
#' @param tbl1 A tibble containing the first set of summary statistics (e.g., coefficients, p-values, CI) for each variant.
#'   Expected columns are:
#'   - `marker`: The name of the marker (e.g., variable name).
#'   - `pval`: The p-value for each marker.
#'   - `ci.lower`: The lower bound of the confidence interval.
#'   - `ci.upper`: The upper bound of the confidence interval.
#'   - `est`: The estimated coefficient.
#' @param tbl2 A tibble containing the second set of summary statistics for each variant (same format as tbl1).
#' @param single_plot A boolean indicating whether to make a single combined plot (TRUE), or plot each dataset side-by-side (FALSE).
#' @param title1 Title for tbl1 data. If single_plot, this will be the legend label for tbl1 data; otherwise it will be the title for the tbl1 plot.
#' @param title2 Title for tbl2 data. If single_plot, this will be the legend label for tbl2 data; otherwise it will be the title for the tbl2 plot.
#' @param title The main title of the combined plot, if single_plot is TRUE.
#' @param sig For individual plots, the p-value threshold for data points to highlight as significant. Defaults to 0.05.
#' @param colors For combined plot, a vector of two colors to represent the two input datasets.
#' @param x_title The title for the x-axis. Defaults to "Coefficient (95% CI)".
#' @param y_title The title for the y-axis. Defaults to "Variant".
#' @param axis_label_size The font size of the axis labels. Defaults to 9.
#' 
#' @return A ggplot object displaying the comparison of the two sets of estimates.
#' @export
compare_estimates <- function(tbl1, 
                              tbl2, 
                              title1=NULL, title2=NULL, title=NULL,
                              sig = 0.05, 
                              colors=c("maroon", "blue4"),
                              x_title="Coefficient (95% CI)",
                              y_title="Variant",
                              axis_label_size=9,
                              single_plot=TRUE,
                              pd=position_dodge(width = 0.8)) {
  
  if (!single_plot) {
    plot1 <- plot_estimates(tbl1, 
                            sig = sig, 
                            sig_colors=c("grey", colors[1]),
                            x_title=x_title,
                            y_title=y_title, 
                            title=title1,
                            axis_label_size=axis_label_size)
    
    
    plot2 <- plot_estimates(tbl2, 
                            sig = sig, 
                            sig_colors=c("grey", colors[2]),
                            x_title=x_title,
                            y_title=y_title, 
                            title=title2,
                            axis_label_size=axis_label_size)
    
    plot <- plot1 + plot2
  }
  
  else {
    tbl1 <- tbl1 %>% mutate(group=title1)
    
    plot <- tbl2 %>% mutate(group=title2) %>% full_join(tbl1) %>%
      mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>%
      filter(marker != "(Intercept)") %>%
      ggplot(aes(y = marker)) +
      geom_vline(xintercept = 0) +
      geom_linerange(aes(xmin = ci.lower, xmax = ci.upper, col=as.factor(group)), position=pd) +
      geom_point(aes(x = est,  col=as.factor(group)), position=pd) +
      scale_color_manual(values = colors) +
      theme_light() +
      theme(axis.text.x=element_text(size=axis_label_size),
            axis.text.y=element_text(size=axis_label_size)) +
      labs(title=title,
           x = x_title,
           y = y_title,
           col = "Group")
  }
  
  return(plot)
}

#' @noRd
#' @method autoplot model_summary
#' @export
autoplot.model_summary <- function(object, sig = 0.05, 
                                   sig_colors=c("grey", "blue4"),
                                   x_title="Coefficient (95% CI)",
                                   y_title="Variant", 
                                   title=NULL,
                                   axis_label_size=9) {
  plot_estimates(model_summary,
                 sig = sig, 
                 sig_colors=sig_colors,
                 x_title=x_title,
                 y_title=y_title, 
                 title=title,
                 axis_label_size=axis_label_size)
}
  
  

#' Extract Details from a logistf Model
#' 
#' This function extracts and formats the estimates, confidence intervals, and p-values from a fitted logistf model.
#' 
#' @param model A fitted logistf model object.
#' 
#' @return A tibble containing the estimates, confidence intervals, and p-values for each predictor in the model.
#' 
#' Example
#' library(logistf)
#' model <- logistf(R ~ ., data=dat)
#' model_details <- logistf_details(model)
#' autoplot(model_details)
#'
#' @export
logistf_details <- function(model) {
  
  # Create summary tibble
  model_summary <- cbind(est=model$coefficients, 
                         ci.lower=model$ci.lower, 
                         ci.upper=model$ci.upper, 
                         pval=model$prob) %>%
    as_tibble(rownames="marker")

  structure(model_summary, class = c("model_summary", class(model_summary)))
}

#' @noRd
#' @export
print.model_summary <- function(x, ...) {
  class(x) <- class(x)[!class(x) == "model_summary"]
  print(x, ...)
  message("Use ggplot2::autoplot() on this output to visualise")
}

#' Extract Details from a Generalized Linear Model
#' 
#' This function extracts and formats the estimates, confidence intervals, and p-values from a fitted glm model.
#' 
#' @param model A fitted glm model object.
#' 
#' @return A tibble containing the estimates, confidence intervals, and p-values for each predictor in the model.
#' 
#' Example
#' model <- glm(R ~ ., data=dat, family = binomial(link = "logit"))
#' model_details <- glm_details(model)
#' autoplot(model_details)
#' @importFrom stats confint
#' @export
glm_details <- function(model) {
  
  # get CI data
  ci <- stats::confint(model) %>% 
    as_tibble(rownames="marker") %>%
    rename(ci.lower=`2.5 %`, ci.upper=`97.5 %`)
  
  # Create summary tibble
  model_summary <- summary(model)$coef %>% 
    as_tibble(rownames="marker") %>% 
    rename(est=Estimate, pval=`Pr(>|z|)`) %>% 
    select(marker, est, pval) %>% 
    left_join(ci)
  
  structure(model_summary, class = c("model_summary", class(model_summary)))
}


#' AMR Logistic Regression Analysis
#'
#' Performs logistic regression to analyze the relationship between genetic markers and phenotype (R, and NWT) for a specified antibiotic.
#'
#' @param geno_table A data frame containing the genotype data.
#' @param pheno_table A data frame containing the phenotypic data.
#' @param antibiotic A character string specifying the antibiotic to model using logistic regression.
#' @param drug_class_list A vector of drug class names. Used to subset the relevant markers for analysis.
#' @param geno_sample_col (Optional) A character string specifying the column in `geno_table` that identifies the sample IDs. Defaults to `NULL`.
#' @param pheno_sample_col (Optional) A character string specifying the column in `pheno_table` that identifies the sample IDs. Defaults to `NULL`.
#' @param sir_col (Optional) A character string specifying the column in `pheno_table` that contains the phenotype values (e.g., resistance/susceptibility). Defaults to `"pheno"`.
#' @param ecoff_col (Optional) A character string specifying the column in `pheno_table` containing the ECOFF (epidemiological cutoff) values. Defaults to `"ecoff"`.
#' @param maf (Optional) An integer specifying the minimum allele frequency (MAF) threshold. Markers with a MAF lower than this value will be excluded. Defaults to 10.
#' @param glm (Optional) Change to TRUE to fit model with glm. Otherwise fit model with logistf (default).
#' @param single_plot (Optional) A logical value. If `TRUE`, a single plot is produced comparing the estimates for resistance (`R`) and non-resistance (`NWT`). Otherwise, two plots are printed side-by-side. Defaults to `TRUE`.
#' @param colors (Optional) A vector of two colors, to use for R and NWT models in the plots. Defaults to `c("maroon", "blue4")`.
#' @param axis_label_size (Optional) A numeric value controlling the size of axis labels in the plot. Defaults to 9.
#'
#' @return A list with three components:
#' \item{bin_mat}{The binary matrix of genetic data and phenotypic resistance information.}
#' \item{modelR}{The fitted logistic regression model for resistance (`R`).}
#' \item{modelNWT}{The fitted logistic regression model for non-resistance (`NWT`).}
#' \item{plot}{A ggplot object comparing the estimates for resistance and non-resistance with corresponding statistical significance indicators.}
#'
#' @examples
#' # Example usage of the amr_logistic function
#' result <- amr_logistic(geno_table = import_amrfp(ecoli_geno_raw, "Name"), 
#'                       pheno_table = ecoli_ast, 
#'                       antibiotic = "Ciprofloxacin", 
#'                       drug_class_list = c("Quinolones"), 
#'                       maf = 10)
#'
#' # To access the plot:
#' print(result$plot)
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import logistf
#' @export
amr_logistic <- function(geno_table, pheno_table, antibiotic, drug_class_list,
                         geno_sample_col = NULL, pheno_sample_col = NULL,
                         sir_col = "pheno", ecoff_col = "ecoff",
                         maf=10, glm=FALSE, single_plot=TRUE, 
                         colors=c("maroon", "blue4"),
                         axis_label_size=9) {
  
  bin_mat <- get_binary_matrix(geno_table = geno_table, 
                               pheno_table = pheno_table, 
                               antibiotic = antibiotic, 
                               drug_class_list = drug_class_list, 
                               geno_sample_col = geno_sample_col,
                               pheno_sample_col = pheno_sample_col,
                               sir_col = sir_col, 
                               ecoff_col = ecoff_col) 
  
  if (glm) {
    print ("Fitting logistic regression models using glm")
    modelR <- glm(R ~ ., data=bin_mat %>% select(-c(id,pheno,NWT)) %>% select(where(~ sum(., na.rm=T) >= maf)), family=binomial(link="logit"))
    modelR <- glm_details(modelR)
    modelNWT <- glm(NWT ~ ., data=bin_mat %>% select(-c(id,pheno,R)) %>% select(where(~ sum(., na.rm=T) >= maf)), family=binomial(link="logit"))
    modelNWT <- glm_details(modelNWT)
  }
  else {
    print ("Fitting logistic regression models using logistf")
    modelR <- logistf::logistf(R ~ ., data=bin_mat %>% select(-c(id,pheno,NWT)) %>% select(where(~ sum(., na.rm=T) >= maf)))
    modelR <- logistf_details(modelR)
    modelNWT <- logistf::logistf(NWT ~ ., data=bin_mat %>% select(-c(id,pheno,R)) %>% select(where(~ sum(., na.rm=T) >= maf)))
    modelNWT <- logistf_details(modelNWT)
  }
  
  plot <- compare_estimates(modelR, modelNWT, 
                            single_plot=single_plot, 
                            title1="R", title2="NWT", 
                            colors=colors, axis_label_size=axis_label_size)
  if (single_plot) {
          ggtitle(label=paste("R and NWT for",antibiotic), 
                  subtitle=paste("for", paste(drug_class_list, collapse=","), "markers present in at least", maf, "samples"))
  }
  
  print(plot)
  
  return(list(bin_mat=bin_mat,
              modelR=modelR,
              modelNWT=modelNWT,
              plot=plot))
}


merge_logreg_soloppv <- function(model, solo_stats, title=NULL) {
  
  combined <- model %>% 
    full_join(solo_stats, by="marker", suffix=c(".est",".ppv")) %>%
    filter(marker != "(Intercept)")
  
  plot <- plot_combined_stats(combined, title=paste(title))
  
  print(plot)
  
  return(list(combined=combined, plot=plot))
}

plot_combined_stats <- function(combined_stats, sig=0.05, title=NULL) {
  combined_stats %>% 
    mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>% 
    ggplot(aes(y=est, x=ppv, col=as.factor(sig_binary))) + 
    geom_point() + 
    geom_linerange(aes(xmin = ci.lower.ppv, xmax = ci.upper.ppv)) +
    geom_linerange(aes(ymin = ci.lower.est, ymax = ci.upper.est)) +
    labs(y="Logistic regression coefficient", x="PPV", 
         col = paste0("logreg p<", sig),
         title=title) + 
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=0.5) + 
    theme_bw()
}
