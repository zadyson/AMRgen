library(dplyr)
library(tidyr)
library(AMR)
library(rvest)

# set up a local environment to store the website results in:
AMRgen_env <- new.env()

# gets all antbiotic options from mic.eucast.org:
eucast_supported_ab_distributions <- function(...) {
  font_url <- get("font_url", envir = asNamespace("AMR"))

  if (is.null(AMRgen_env$eucast_ab_select_list)) {
    if (interactive()) message("Retrieving list of antimicrobials from ", font_url("https://mic.eucast.org", "mic.eucast.org"), "...", appendLF = FALSE)
    url <- "https://mic.eucast.org/search/?search[method]=mic"
    page <- read_html(url)
    select_list <- page %>% html_element("#search_antibiotic") %>% html_children()
    select_values <- select_list %>% html_attrs() %>% unlist()
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

# get the distributions:
get_eucast_amr_distribution <- function(ab, mo = NULL, method = "MIC", as_freq_table = TRUE) {
  font_url <- get("font_url", envir = asNamespace("AMR"))
  font_italic <- get("font_italic", envir = asNamespace("AMR"))

  # retrieve available antimicrobials online
  eucast_supported_ab_distributions(invisible = TRUE)

  method <- trimws(tolower(method)[1])
  method <- ifelse(method == "mic", "mic", "diff")
  ab <- trimws(toupper(ab)[1])
  ab_coerced <- as.ab(ab)
  mo_coerced <- NULL
  if (!is.null(mo)) {
    mo <- trimws(toupper(mo)[1])
    mo_coerced <- as.mo(mo)
  }
  if (ab != ab_coerced) {
    if (interactive()) message("Returning antimicrobial wild type distributions for ",
                               ab_name(ab_coerced, language = NULL, tolower = TRUE), " (", ab_coerced, ", ", ab_atc(ab_coerced, only_first = TRUE), ")",
                               ifelse(!is.null(mo_coerced), paste0(" in ", font_italic(suppressWarnings(mo_name(mo_coerced, language = NULL, keep_synonyms = TRUE)))), ""))
  }
  ab <- ab_coerced
  ab_index <- names(AMRgen_env$eucast_ab_select_list)[AMRgen_env$eucast_ab_select_list == ab]
  url <- paste0("https://mic.eucast.org/search/?search[method]=", method, "&search[antibiotic]=", ab_index, "&search[limit]=999")
  if (interactive()) message("From: ", font_url(url))
  tbl <- read_html(url) %>%
    html_element("#search-results-table") %>%
    html_table()
  colnames(tbl)[1] <- "microorganism"
  tbl <- tbl %>%
    mutate(microorganism_code = suppressWarnings(as.mo(microorganism, language = NULL, keep_synonyms = TRUE, info = FALSE)),
                  .after = 1)

  if (!is.null(mo)) {
    mo <- as.mo(mo)
    tbl <- tbl %>%
      filter(microorganism_code == mo) %>%
      select(-1, -2)
  }

  col_nms <- colnames(tbl)
  col_nms <- gsub("[^a-zA-Z0-9.]+", "_", col_nms)
  col_nms <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", col_nms)
  col_nms <- gsub("([a-zA-Z])([0-9])", "\\1_\\2", col_nms)
  col_nms <- tolower(col_nms)
  colnames(tbl) <- col_nms
  colnames(tbl)[colnames(tbl) == "_t_ecoff"] <- "ecoff"

  if (isTRUE(as_freq_table)) {
    tbl <- tbl %>%
      select(matches("^([0-9.]+|microorganism|microorganism_code)$")) %>%
      pivot_longer(-matches("^(microorganism|microorganism_code)$"),
                          names_to = ifelse(method == "mic", "mic", "disk_diffusion"),
                          values_to = "count")
    if (method == "mic") {
      tbl$mic <- as.mic(tbl$mic)
    } else {
      tbl$disk_diffusion <- as.disk(tbl$disk_diffusion)
    }
  }

  if ("count" %in% colnames(tbl)) {
    tbl$count[tbl$count == "No data available"] <- 0
    tbl$count <- as.integer(tbl$count)
  }

  tbl
}

# retrieve ECOFF
get_ecoff <- function(ab, mo = NULL) {

}

