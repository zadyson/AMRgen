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
#' @param sig_colors For combined plot, a vector of two colors to represent the two input datasets.
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
                              sig_colors=c("maroon", "blue4"),
                              x_title="Coefficient (95% CI)",
                              y_title="Variant",
                              axis_label_size=9,
                              single_plot=TRUE) {
  
  if (!single_plot) {
    plot1 <- plot_estimates(tbl1, 
                            sig = sig, 
                            sig_colors=c("grey", sig_colors[1]),
                            x_title=x_title,
                            y_title=y_title, 
                            title=title1,
                            axis_label_size=axis_label_size)
    
    
    plot2 <- plot_estimates(tbl2, 
                            sig = sig, 
                            sig_colors=c("grey", sig_colors[2]),
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
      geom_linerange(aes(xmin = ci.lower, xmax = ci.upper, col=as.factor(group)), position=position_dodge(width=0.8)) +
      geom_point(aes(x = est,  col=as.factor(group)), position=position_dodge(width=0.8)) +
      scale_color_manual(values = sig_colors) +
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
autoplot.model_summary <- function(object, ...) {
  plot_estimates(model_summary)
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
#' logistf_details(model)
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
#' glm_details(model)
#' 
#' @export
glm_details <- function(model) {
  
  # get CI data
  ci <- confint(model) %>% 
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
