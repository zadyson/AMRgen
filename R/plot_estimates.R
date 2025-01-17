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
#' @param sig A numeric value specifying the significance threshold for p-values. 
#'   Defaults to 0.05.
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
plot_estimates <- function(tbl, sig=0.05) {
  plot <- tbl %>% 
    mutate(sig_binary = if_else(pval < sig, TRUE, FALSE)) %>%
    filter(marker != "(Intercept)") %>%
    mutate(marker = gsub("`", "", marker)) %>%
    ggplot(aes(y = marker, col = sig_binary)) + 
    geom_vline(xintercept = 0) + 
    geom_linerange(aes(xmin = ci.lower, xmax = ci.upper)) + 
    geom_point(aes(x = est)) + 
    scale_color_manual(values = c("grey", "blue4")) + 
    theme_light() + 
    labs(x = "Logistic regression coefficient (95% CI)", col = paste0("p<", sig)) + 
    labs(y = "Variant")
  
  return(plot)
}
