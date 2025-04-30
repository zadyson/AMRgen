#' Retrieve Available Antimicrobial Wild Type Distributions from EUCAST
#'
#' Run this function to get an updated list of antimicrobial distributions currently supported by EUCAST. This retrieves live info from <https://mic.eucast.org>.
#' @param ... arguments passed on to function, currently unused
#' @importFrom rvest read_html html_element html_children html_text2 html_attrs
#' @importFrom AMR ab_name
#' @importFrom dplyr %>%
#' @importFrom tidyr pivot_longer
#' @examples
#' eucast_supported_ab_distributions()
#' @export
eucast_supported_ab_distributions <- function(...) {
  font_url <- get("font_url", envir = asNamespace("AMR"))

  if (is.null(AMRgen_env$eucast_ab_select_list)) {
    if (interactive()) message("Retrieving list of antimicrobials from ", font_url("https://mic.eucast.org", "mic.eucast.org"), "...", appendLF = FALSE)
    url <- "https://mic.eucast.org/search/?search[method]=mic"
    page <- read_html(url)
    select_list <- page %>%
      html_element("#search_antibiotic") %>%
      html_children()
    select_values <- select_list %>%
      html_attrs() %>%
      unlist()
    select_names <- select_list %>% html_text2()
    select_names <- select_names[!grepl("...", names(select_values), fixed = TRUE)]
    select_names_AMR <- as.character(as.ab(select_names, flag_multiple_results = FALSE, info = FALSE, fast_mode = TRUE))
    select_values <- select_values[!is.na(select_names_AMR)]
    select_names_AMR <- select_names_AMR[!is.na(select_names_AMR)]
    names(select_names_AMR) <- select_values
    AMRgen_env$eucast_ab_select_list <- select_names_AMR
    if (interactive()) message("OK")
  }

  if (isTRUE(list(...)$invisible)) {
    return(invisible())
  }
  out <- ab_name(AMRgen_env$eucast_ab_select_list)
  names(out) <- AMRgen_env$eucast_ab_select_list
  sort(out)
}
