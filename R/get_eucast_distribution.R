#' Retrieve Antimicrobial Wild Type Distributions from EUCAST
#'
#' These functions allow retrieval of antimicrobial wild type distributions, live from [eucast.org](https://mic.eucast.org).
#' @param ab antimicrobial, can be anything understood by [ab_name()]
#' @param mo microorganism, can be anything understood by [mo_name()] (can be left blank)
#' @param method either `"MIC"` or `"disk"`/`"diff"`
#' @param as_freq_table either `TRUE` or `FALSE`, to return result as frequency table
#' @details
#' In December 2024, EUCAST had 176 distributions available, namely for these antimicrobials:
#'
#' Amikacin, amoxicillin, amoxicillin/clavulanic acid, amphotericin B, ampicillin, ampicillin/sulbactam, anidulafungin, apramycin, aspoxicillin, avilamycin, azithromycin, aztreonam, aztreonam/avibactam, bacitracin, bedaquiline, benzylpenicillin, capreomycin, cefaclor, cefadroxil, cefalexin, cefaloridine, cefalotin, cefapirin, cefazolin, cefdinir, cefepime, cefepime/tazobactam, cefepime/zidebactam, cefiderocol, cefixime, cefoperazone, cefoperazone/sulbactam, cefoselis, cefotaxime, cefotetan, cefovecin, cefoxitin, cefpirome, cefpodoxime, cefpodoxime/clavulanic acid, cefquinome, ceftaroline, ceftazidime, ceftazidime/avibactam, ceftibuten, ceftiofur, ceftobiprole, ceftolozane/tazobactam, ceftriaxone, cefuroxime, cephradine, chloramphenicol, chlortetracycline, ciprofloxacin, clarithromycin, clavulanic acid, clinafloxacin, clindamycin, clofazimine, cloxacillin, colistin, cycloserine, dalbavancin, danofloxacin, daptomycin, delafloxacin, delamanid, dicloxacillin, difloxacin, doripenem, doxycycline, enrofloxacin, eravacycline, ertapenem, erythromycin, ethambutol, ethionamide, faropenem, fidaxomicin, florfenicol, flucloxacillin, fluconazole, flucytosine, flumequine, fosfomycin, fusidic acid, gamithromycin, gatifloxacin, gemifloxacin, gentamicin, imipenem, imipenem/relebactam, isavuconazole, isoniazid, itraconazole, kanamycin, ketoconazole, lefamulin, levofloxacin, lincomycin, linezolid, loracarbef, marbofloxacin, mecillinam, meropenem, meropenem/vaborbactam, metronidazole, micafungin, minocycline, moxifloxacin, mupirocin, nalidixic acid, narasin, neomycin, netilmicin, nitrofurantoin, nitroxoline, norfloxacin, norvancomycin, ofloxacin, omadacycline, orbifloxacin, oritavancin, oxacillin, oxolinic acid, oxytetracycline, pefloxacin, phenoxymethylpenicillin, piperacillin, piperacillin/tazobactam, pirlimycin, posaconazole, pradofloxacin, pristinamycin, pyrazinamide, quinupristin/dalfopristin, retapamulin, rezafungin, rifabutin, rifampicin, roxithromycin, secnidazole, sitafloxacin, spectinomycin, spiramycin, streptomycin, sulbactam, sulfadiazine, sulfamethoxazole, sulfisoxazole, tedizolid, teicoplanin, telavancin, telithromycin, temocillin, terbinafine, tetracycline, thiamphenicol, tiamulin, ticarcillin, ticarcillin/clavulanic acid, tigecycline, tildipirosin, tilmicosin, tobramycin, trimethoprim, trimethoprim/sulfamethoxazole, tulathromycin, tylosin, tylvalosin, vancomycin, viomycin, and voriconazole.
#' @importFrom rvest read_html html_element html_children html_text2 html_table html_attrs
#' @importFrom AMR as.ab ab_name as.mo as.mic as.disk
#' @importFrom dplyr mutate filter select matches %>%
#' @importFrom tidyr pivot_longer
#' @rdname get_eucast_amr_distribution
#' @export
#' @examples
#' get_eucast_mic_distribution("cipro")
#' get_eucast_mic_distribution("cipro", as_freq_table = FALSE)
#'
#' get_eucast_mic_distribution("cipro", "K. pneumoniae")
#' get_eucast_disk_distribution("cipro", "K. pneumoniae")
#'
#'
#' # Plotting ----------------------------------------------------------------
#'
#' library(ggplot2)
#'
#' mic_data <- get_eucast_mic_distribution("cipro", "K. pneumoniae")
#' mics <- rep(mic_data$mic, mic_data$count)
#' autoplot(mics, ab = "cipro", mo = "K. pneumoniae", title = "Look at my MICs!")
#'
#' disk_data <- get_eucast_disk_distribution("cipro", "K. pneumoniae")
#' disks <- rep(disk_data$disk_diffusion, disk_data$count)
#' autoplot(disks, ab = "cipro", mo = "K. pneumoniae", title = "Look at my diffusion zones!")
#'
#'
#' # Comparing With User Values ----------------------------------------------
#'
#' my_mic_values <- AMR::random_mic(500)
#' comparison <- compare_mic_with_eucast(my_mic_values, ab = "cipro", mo = "K. pneumoniae")
#' comparison
#' autoplot(comparison)
get_eucast_amr_distribution <- function(ab, mo = NULL, method = "MIC", as_freq_table = TRUE) {
  if (is.null(AMRgen_env$eucast_ab_select_list)) {
    message("Retrieving list of antimicrobials from eucast.org...", appendLF = FALSE)
    url <- "https://mic.eucast.org/search/?search[method]=mic"
    page <- read_html(url)
    select_list <- page %>% html_element("#search_antibiotic") %>%  html_children()
    select_values <- select_list %>% html_attrs() %>% unlist()
    select_names <- select_list %>% html_text2()
    select_names <- select_names[!grepl("...", names(select_values), fixed = TRUE)]
    select_names_AMR <- as.character(as.ab(select_names, flag_multiple_results = FALSE, info = FALSE, fast_mode = TRUE))
    select_values <- select_values[!is.na(select_names_AMR)]
    select_names_AMR <- select_names_AMR[!is.na(select_names_AMR)]
    names(select_names_AMR) <- select_values
    AMRgen_env$eucast_ab_select_list <- select_names_AMR
    message("OK\n")
  }

  method <- trimws(tolower(method)[1])
  method <- ifelse(method == "mic", "mic", "diff")
  ab <- trimws(toupper(ab)[1])
  ab_coerced <- as.ab(ab)
  if (ab != ab_coerced) {
    message("Returning antimicrobial wild type distributions for ", ab_name(ab_coerced, language = NULL, tolower = TRUE))
  }
  ab <- ab_coerced
  ab_index <- names(AMRgen_env$eucast_ab_select_list)[AMRgen_env$eucast_ab_select_list == ab]
  url <- paste0("https://mic.eucast.org/search/?search[method]=", method, "&search[antibiotic]=", ab_index, "&search[limit]=999")
  font_url <- get("font_url", envir = asNamespace("AMR"))
  message("From: ", font_url(url))
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

  tbl
}

#' @rdname get_eucast_amr_distribution
#' @export
get_eucast_mic_distribution <- function(ab, mo = NULL, as_freq_table = TRUE) {
  get_eucast_amr_distribution(ab = ab, mo = mo, method = "mic", as_freq_table = as_freq_table)
}

#' @rdname get_eucast_amr_distribution
#' @export
get_eucast_disk_distribution <- function(ab, mo = NULL, as_freq_table = TRUE) {
  get_eucast_amr_distribution(ab = ab, mo = mo, method = "diff", as_freq_table = as_freq_table)
}

#' @rdname get_eucast_amr_distribution
#' @param mics MIC values of class [`mic`][AMR::as.mic()]
#' @importFrom AMR is.mic
#' @importFrom dplyr select filter as_tibble
#' @details
#' The `compare_*_with_eucast()` functions allow to compare a user range with EUCAST distributions. Use [autoplot()] on the output to visualise the results.
#' @export
compare_mic_with_eucast <- function(mics, ab, mo = NULL) {
  if (!is.mic(mics)) {
    stop("`mics` must be of class 'mic', please see ?AMR::as.mic", call. = FALSE)
  }
  if (is.null(mo)) {
    stop("`mo` must be filled in for comparison", call. = FALSE)
  }
  distr <- get_eucast_mic_distribution(ab = ab, mo = mo, as_freq_table = TRUE)
  vals <- rep(distr[[1]], distr[[2]])
  user_mic <- mics |>
    table() |>
    as.data.frame() |>
    select(value = mics, user = Freq)
  eucast_mic <- vals |>
    table() |>
    as.data.frame() |>
    select(eucast = Freq)
  total <- user_mic |>
    cbind(eucast_mic) |>
    filter(user + eucast > 0) |>
    as_tibble()
  structure(total, class = c("compare_eucast", class(total)))
}

#' @rdname get_eucast_amr_distribution
#' @param disks Disk diffusion values of class [`disk`][AMR::as.disk()]
#' @importFrom AMR is.disk
#' @importFrom dplyr select filter as_tibble
#' @export
compare_disk_with_eucast <- function(disks, ab, mo = NULL) {
  if (!is.disk(mics)) {
    stop("`disks` must be of class 'disk', please see ?AMR::as.disk", call. = FALSE)
  }
  if (is.null(mo)) {
    stop("`mo` must be filled in for comparison", call. = FALSE)
  }
  distr <- get_eucast_disk_distribution(ab = ab, mo = mo, as_freq_table = TRUE)
  vals <- rep(distr[[1]], distr[[2]])
  user_disk <- disks |>
    table() |>
    as.data.frame() |>
    select(value = disks, user = Freq)
  eucast_disk <- vals |>
    table() |>
    as.data.frame() |>
    select(eucast = Freq)
  total <- user_disk |>
    cbind(eucast_disk) |>
    filter(user + eucast > 0) |>
    as_tibble()
  structure(total, class = c("compare_eucast", class(total)))
}

#' @noRd
#' @export
print.compare_eucast <- function(x, ...) {
  class(x) <- class(x)[!class(x) == "compare_eucast"]
  print(x, ...)
  message("Use ggplot2::autoplot() on this output to visualise")
}

#' @noRd
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 autoplot ggplot aes labs geom_col
#' @importFrom AMR scale_x_mic
#' @method autoplot compare_eucast
#' @export
autoplot.compare_eucast <- function(object, ...) {
  long <- object |>
    mutate(User = user / sum(user),
           EUCAST = eucast / sum(eucast)) |>
    select(-user, -eucast) |>
    pivot_longer(-value, names_to = "Source", values_to = "count")
  ggplot(long, aes(x = value, y = count, fill = Source)) +
    geom_col(position = "dodge") +
    labs(x = "Measurement Value", y = "Density") +
    scale_x_mic(keep_operators = "none")
}
