#' Plot Estimates from a Table of Results
#'
#' This function creates a ggplot object visualizing logistic regression coefficients with their 95% confidence intervals. Significant markers are highlighted based on a specified p-value threshold.
#' @param tbl A data frame or tibble containing the logistic regression results. Expected columns are:
#' - `marker`: The name of the marker (e.g., variable name).
#' - `pval`: The p-value for each marker.
#' - `ci.lower`: The lower bound of the confidence interval.
#' - `ci.upper`: The upper bound of the confidence interval.
#' - `est`: The estimated coefficient.
#' @param sig (optional) The significance threshold for p-values. Defaults to 0.05.
#' @param sig_colors (optional) A vector of two colors to represent significant and non-significant estimates.
#' @param x_title (optional) The title for the x-axis. Defaults to "Coefficient (95% CI)".
#' @param y_title (optional) The title for the y-axis. Defaults to "Variant".
#' @param title (optional) The main title of the plot. If NULL, no title is added.
#' @param axis_label_size (optional) The font size of the axis labels. Defaults to 9.
#' @param marker_order (optional) Vector indicating the order of the markers to be plotted on the y-axis.
#' @importFrom dplyr filter if_else mutate
#' @importFrom ggplot2 aes element_text geom_linerange geom_point geom_vline ggplot labs scale_color_manual scale_y_discrete theme theme_light
#' @return A ggplot object showing the logistic regression coefficients with confidence intervals. Significant markers (p-value < `sig`) are colored differently.
#' @export
#' @examples
#' tbl <- tibble::tibble(
#'   marker = c("(Intercept)", "var1", "var2", "var3"),
#'   pval = c(0.1, 0.03, 0.2, 0.04),
#'   ci.lower = c(-0.2, 0.1, -0.3, 0.2),
#'   ci.upper = c(0.5, 0.8, 0.4, 1.1),
#'   est = c(0.2, 0.5, 0.1, 0.7)
#' )
#'
#' plot_estimates(tbl)
plot_estimates <- function(tbl, sig = 0.05,
                           sig_colors = c(`FALSE` = "grey", `TRUE` = "blue4"),
                           x_title = "Coefficient (95% CI)",
                           y_title = "Variant",
                           title = NULL,
                           axis_label_size = 9,
                           marker_order = NULL) {
  if (!is.null(sig)) {
    legend_title <- paste0("p<", sig)
    legend_position <- "right"
    tbl <- tbl %>% mutate(sig_binary = if_else(pval < sig, TRUE, FALSE))
  } else {
    legend_title <- ""
    legend_position <- "none"
    tbl <- tbl %>% mutate(sig_binary = TRUE)
  }

  plot <- tbl %>%
    filter(marker != "(Intercept)") %>%
    mutate(marker = gsub("`", "", marker)) %>%
    ggplot(aes(y = marker, col = sig_binary)) +
    geom_vline(xintercept = 0) +
    geom_linerange(aes(xmin = ci.lower, xmax = ci.upper)) +
    geom_point(aes(x = est)) +
    scale_color_manual(values = sig_colors) +
    theme_light() +
    theme(
      axis.text.x = element_text(size = axis_label_size),
      axis.text.y = element_text(size = axis_label_size),
      legend.position = legend_position
    ) +
    labs(
      title = title,
      x = x_title,
      y = y_title,
      col = legend_title
    )

  if (!is.null(marker_order)) {
    plot <- plot + scale_y_discrete(limits = marker_order)
  }

  return(plot)
}

#' Plot to Compare Two Sets of Estimates
#'
#' This function compares two sets of estimates by creating a plot that overlays the estimates and confidence intervals for both sets. It can also display the estimates in two separate plots.
#'
#' @param tbl1 A tibble containing the first set of summary statistics (e.g., coefficients, p-values, CI) for each variant. Expected columns are:
#'  - `marker`: The name of the marker (e.g., variable name).
#'  - `pval`: The p-value for each marker.
#'  - `ci.lower`: The lower bound of the confidence interval.
#'  - `ci.upper`: The upper bound of the confidence interval.
#'  - `est`: The estimated coefficient.
#' @param tbl2 A tibble containing the second set of summary statistics for each variant (same format as tbl1).
#' @param single_plot A boolean indicating whether to make a single combined plot (TRUE), or plot each dataset side-by-side (FALSE).
#' @param title1 Title for tbl1 data. If single_plot, this will be the legend label for tbl1 data; otherwise it will be the title for the tbl1 plot.
#' @param title2 Title for tbl2 data. If single_plot, this will be the legend label for tbl2 data; otherwise it will be the title for the tbl2 plot.
#' @param title (optional) The main title of the combined plot, if single_plot is TRUE.
#' @param sig (optional) For individual plots, the p-value threshold for data points to highlight as significant. Defaults to 0.05.
#' @param colors (optional) For combined plot, a vector of two colors to represent the two input datasets.
#' @param x_title (optional) The title for the x-axis. Defaults to "Coefficient (95% CI)".
#' @param y_title (optional) The title for the y-axis. Defaults to "Variant".
#' @param axis_label_size (optional) The font size of the axis labels. Defaults to 9.
#' @param pd (optional) Position dodge, i.e. spacing for the 2 estimates to be positioned above/below the line. Default 'position_dodge(width = 0.8)'
#' @param marker_order (optional) Vector indicating the order of the markers to be plotted on the y-axis.
#' @importFrom dplyr bind_rows filter if_else mutate
#' @importFrom ggplot2 aes element_text geom_linerange geom_point geom_vline ggplot labs position_dodge scale_color_manual scale_y_discrete theme theme_light
#' @return A ggplot object displaying the comparison of the two sets of estimates.
#' @export
compare_estimates <- function(tbl1,
                              tbl2,
                              title1 = NULL, title2 = NULL, title = NULL,
                              sig = 0.05,
                              colors = c("maroon", "blue4"),
                              x_title = "Coefficient (95% CI)",
                              y_title = "Variant",
                              axis_label_size = 9,
                              single_plot = TRUE,
                              pd = position_dodge(width = 0.8),
                              marker_order = NULL) {
  if (!single_plot) {
    plot1 <- plot_estimates(tbl1,
      sig = sig,
      sig_colors = c("grey", colors[1]),
      x_title = x_title,
      y_title = y_title,
      title = title1,
      axis_label_size = axis_label_size
    )

    plot2 <- plot_estimates(tbl2,
      sig = sig,
      sig_colors = c("grey", colors[2]),
      x_title = x_title,
      y_title = y_title,
      title = title2,
      axis_label_size = axis_label_size
    )

    if (!is.null(marker_order)) {
      plot1 <- plot1 + scale_y_discrete(limits = marker_order)
      plot2 <- plot2 + scale_y_discrete(limits = marker_order)
    }

    plot <- plot1 + plot2
  } else {
    tbl1 <- tbl1 %>% mutate(group = title1)

    plot <- tbl2 %>%
      mutate(group = title2) %>%
      bind_rows(tbl1) %>%
      mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>%
      filter(marker != "(Intercept)") %>%
      ggplot(aes(y = marker)) +
      geom_vline(xintercept = 0) +
      geom_linerange(aes(xmin = ci.lower, xmax = ci.upper, col = as.factor(group)), position = pd) +
      geom_point(aes(x = est, col = as.factor(group)), position = pd) +
      scale_color_manual(values = colors) +
      theme_light() +
      theme(
        axis.text.x = element_text(size = axis_label_size),
        axis.text.y = element_text(size = axis_label_size)
      ) +
      labs(
        title = title,
        x = x_title,
        y = y_title,
        col = "Group"
      )

    if (!is.null(marker_order)) {
      plot <- plot + scale_y_discrete(limits = marker_order)
    }
  }

  return(plot)
}

#' @noRd
#' @method autoplot model_summary
#' @importFrom ggplot2 autoplot
#' @export
autoplot.model_summary <- function(object,
                                   sig = 0.05,
                                   sig_colors = c("grey", "blue4"),
                                   x_title = "Coefficient (95% CI)",
                                   y_title = "Variant",
                                   title = NULL,
                                   axis_label_size = 9,
                                   ...) {
  plot_estimates(
    object,
    sig = sig,
    sig_colors = sig_colors,
    x_title = x_title,
    y_title = y_title,
    title = title,
    axis_label_size = axis_label_size
  )
}



#' Extract Details from a logistf Model
#'
#' This function extracts and formats the estimates, confidence intervals, and p-values from a fitted logistf model.
#' @param model A fitted logistf model object.
#' @importFrom dplyr as_tibble
#' @return A tibble containing the estimates, confidence intervals, and p-values for each predictor in the model.
#' @export
#' @examples
#' library(logistf)
#' library(ggplot2)
#' model <- logistf(case ~ age + oc + vic + vicl + vis + dia, data = sex2)
#' model_details <- logistf_details(model)
#' autoplot(model_details)
logistf_details <- function(model) {
  # Create summary tibble
  model_summary <- cbind(
    est = model$coefficients,
    ci.lower = model$ci.lower,
    ci.upper = model$ci.upper,
    pval = model$prob
  ) %>%
    as_tibble(rownames = "marker")

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
#' @param model A fitted glm model object.
#' @importFrom dplyr as_tibble left_join rename select
#' @return A tibble containing the estimates, confidence intervals, and p-values for each predictor in the model.
#' @export
#' @examples
#' # Generate example data
#' set.seed(1)
#' dat <- data.frame(
#'   R = rbinom(100, 1, 0.3),
#'   geneA = rbinom(100, 1, 0.2),
#'   geneB = rbinom(100, 1, 0.5),
#'   geneC = rbinom(100, 1, 0.3)
#' )
#'
#' # Fit logistic regression model and extract details
#' model <- glm(R ~ ., data = dat, family = binomial(link = "logit"))
#' model_details <- glm_details(model)
#'
#' # Plot model summary
#' ggplot2::autoplot(model_details)
glm_details <- function(model) {
  # get CI data
  ci <- stats::confint(model) %>%
    as_tibble(rownames = "marker") %>%
    rename(ci.lower = `2.5 %`, ci.upper = `97.5 %`)

  # Create summary tibble
  model_summary <- summary(model)$coef %>%
    as_tibble(rownames = "marker") %>%
    rename(est = Estimate, pval = `Pr(>|z|)`) %>%
    select(marker, est, pval) %>%
    left_join(ci, by = "marker")

  structure(model_summary, class = c("model_summary", class(model_summary)))
}


#' AMR Logistic Regression Analysis
#'
#' Performs logistic regression to analyze the relationship between genetic markers and phenotype (R, and NWT) for a specified antibiotic.
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
#' @param marker_col (Optional) Name of the column containing the marker identifiers, whose unique values will be treate as predictors in the regression. Defaults to `"marker"`.
#' @importFrom dplyr any_of select where
#' @importFrom ggplot2 ggtitle
#' @importFrom stats glm
#' @return A list with the following components:
#' - `bin_mat`: The binary matrix of genetic data and phenotypic resistance information.
#' - `modelR`: The fitted logistic regression model for resistance (`R`).
#' - `modelNWT`: The fitted logistic regression model for non-resistance (`NWT`).
#' - `plot`: A ggplot object comparing the estimates for resistance and non-resistance with corresponding statistical significance indicators.
#' @export
#' @examples
#' # Example usage of the amr_logistic function
#' result <- amr_logistic(
#'   geno_table = import_amrfp(ecoli_geno_raw, "Name"),
#'   pheno_table = ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   maf = 10
#' )
#' # To access the plot:
#' print(result$plot)
amr_logistic <- function(geno_table, pheno_table, antibiotic, drug_class_list,
                         geno_sample_col = NULL, pheno_sample_col = NULL,
                         sir_col = "pheno", ecoff_col = "ecoff",
                         maf = 10, glm = FALSE, single_plot = TRUE,
                         colors = c("maroon", "blue4"),
                         axis_label_size = 9, marker_col="marker") {

  bin_mat <- get_binary_matrix(
    geno_table = geno_table,
    pheno_table = pheno_table,
    antibiotic = antibiotic,
    drug_class_list = drug_class_list,
    geno_sample_col = geno_sample_col,
    pheno_sample_col = pheno_sample_col,
    sir_col = sir_col,
    ecoff_col = ecoff_col,
    marker_col = marker_col
  )

  if (glm) {
    print("Fitting logistic regression models using glm")
    modelR <- glm(R ~ ., data = bin_mat %>% select(-any_of(c("id", "pheno", "mic", "disk", "NWT"))) %>% select(where(~ sum(., na.rm = T) >= maf)), family = stats::binomial(link = "logit"))
    modelR <- glm_details(modelR) %>% mutate(marker=gsub("\\.\\.", ":", marker)) %>% mutate(marker=gsub("`", "", marker))
    modelNWT <- glm(NWT ~ ., data = bin_mat %>% select(-any_of(c("id", "pheno", "mic", "disk", "R"))) %>% select(where(~ sum(., na.rm = T) >= maf)), family = stats::binomial(link = "logit"))
    modelNWT <- glm_details(modelNWT) %>% mutate(marker=gsub("\\.\\.", ":", marker)) %>% mutate(marker=gsub("`", "", marker))
  } else {
    print("Fitting logistic regression models using logistf")
    modelR <- logistf::logistf(R ~ ., data = bin_mat %>% select(-any_of(c("id", "pheno", "mic", "disk", "NWT"))) %>% select(where(~ sum(., na.rm = T) >= maf)), pl = FALSE)
    modelR <- logistf_details(modelR) %>% mutate(marker=gsub("\\.\\.", ":", marker)) %>% mutate(marker=gsub("`", "", marker))
    modelNWT <- logistf::logistf(NWT ~ ., data = bin_mat %>% select(-any_of(c("id", "pheno", "mic", "disk", "R"))) %>% select(where(~ sum(., na.rm = T) >= maf)), pl = FALSE)
    modelNWT <- logistf_details(modelNWT) %>% mutate(marker=gsub("\\.\\.", ":", marker)) %>% mutate(marker=gsub("`", "", marker))
  }

  plot <- compare_estimates(modelR, modelNWT,
    single_plot = single_plot,
    title1 = "R", title2 = "NWT",
    colors = colors, axis_label_size = axis_label_size
  )
  if (single_plot) {
    ggtitle(
      label = paste("R and NWT for", antibiotic),
      subtitle = paste("for", paste(drug_class_list, collapse = ","), "markers present in at least", maf, "samples")
    )
  }

  print(plot)

  return(list(
    bin_mat = bin_mat,
    modelR = modelR,
    modelNWT = modelNWT,
    plot = plot
  ))
}


#' Merge Logistic Regression and Solo PPV Statistics
#'
#' This function merges logistic regression model statistics with solo PPV statistics and creates a combined plot.
#' @param model A data frame containing logistic regression model statistics.
#' @param solo_stats A data frame containing solo PPV statistics.
#' @param title An optional title for the plot.
#' @importFrom dplyr filter full_join
#' @return A list containing:
#' - `combined`: A merged data frame of logistic regression and PPV statistics.
#' - `plot`: A ggplot object showing the relationship between model estimates and PPV.
#' @export
#' @examples
#' \dontrun{
#' soloPPV_cipro <- solo_ppv_analysis(ecoli_geno, ecoli_ast,
#'   antibiotic = "Ciprofloxacin",
#'   drug_class_list = c("Quinolones"),
#'   sir_col = "pheno"
#' )
#' logistic_cipro <- amr_logistic(ecoli_geno, ecoli_ast,
#'   "Ciprofloxacin", c("Quinolones"),
#'   maf = 5
#' )
#' allstatsR <- merge_logreg_soloppv(logistic_cipro$modelR,
#'   soloPPV_cipro$solo_stats %>% filter(category == "R"),
#'   title = "Quinolone markers vs Cip R"
#' )
#' }
merge_logreg_soloppv <- function(model, solo_stats, title = NULL) {
  combined <- model %>%
    full_join(solo_stats, by = "marker", suffix = c(".est", ".ppv")) %>%
    filter(marker != "(Intercept)")

  plot <- plot_combined_stats(combined, title = paste(title))

  print(plot)

  return(list(combined = combined, plot = plot))
}

#' Plot Combined Statistics
#'
#' This function creates a plot of combined logistic regression and solo PPV statistics.
#' @param combined_stats A data frame containing combined statistics from logistic regression and solo PPV.
#' @param sig A significance level for the logistic regression p-values. Default is 0.05.
#' @param title An optional title for the plot.
#' @importFrom dplyr if_else mutate
#' @importFrom ggplot2 aes geom_hline geom_linerange geom_point geom_vline ggplot labs theme_bw
#' @return A ggplot2 object representing the combined plot.
#' @export
plot_combined_stats <- function(combined_stats, sig = 0.05, title = NULL) {
  combined_stats %>%
    mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>%
    ggplot(aes(y = est, x = ppv, col = as.factor(sig_binary))) +
    geom_point() +
    geom_linerange(aes(xmin = ci.lower.ppv, xmax = ci.upper.ppv)) +
    geom_linerange(aes(ymin = ci.lower.est, ymax = ci.upper.est)) +
    labs(
      y = "Logistic regression coefficient", x = "PPV",
      col = paste0("logreg p<", sig),
      title = title
    ) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.5) +
    theme_bw()
}


#' Plot Combined Statistics of Logistic Regression and Solo PPV
#'
#' This function creates a plot comparing logistic regression coefficients with PPV values from solo marker analysis. It highlights markers based on statistical significance of the logistic regression.
#' @param combined_stats A data frame containing combined statistics from logistic regression and solo PPV analysis.
#' @param sig A numeric value specifying the significance threshold for p-values. Default is 0.05.
#' @param title An optional character string specifying the plot title.
#' @importFrom dplyr if_else mutate
#' @importFrom ggplot2 aes geom_hline geom_linerange geom_point geom_vline ggplot labs theme_bw
#' @return A ggplot2 object visualizing the relationship between PPV and logistic regression estimates, with confidence intervals and significance annotation.
#' @export
plot_solo_logReg <- function(combined_stats, sig = 0.05, title = NULL) {
  combined_stats %>%
    mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>%
    ggplot(aes(y = est, x = ppv, col = as.factor(sig_binary))) +
    geom_point() +
    geom_linerange(aes(xmin = ci.lower.ppv, xmax = ci.upper.ppv)) +
    geom_linerange(aes(ymin = ci.lower.est, ymax = ci.upper.est)) +
    labs(
      y = "Logistic regression coefficient", x = "PPV",
      col = paste0("logreg p<", sig),
      title = title
    ) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0.5) +
    theme_bw()
}
