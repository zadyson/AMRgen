#' Get and Compare Antimicrobial Wild Type Distributions from EUCAST
#'
#' These functions allow retrieval of antimicrobial wild type distributions, live from [eucast.org](https://mic.eucast.org).
#' @param ab antimicrobial, can be anything understood by [`ab_name()`][AMR::ab_name()]
#' @param mo microorganism, can be anything understood by [`mo_name()`][AMR::mo_name()] (can be left blank)
#' @param method either `"MIC"` or `"disk"`/`"diff"`
#' @param as_freq_table either `TRUE` (default) or `FALSE`, to return result as frequency table
#' @param mics MIC values, will be coerced with [`as.mic()`][AMR::as.mic()]
#' @param disks Disk diffusion values, will be coerced with [`as.disk()`][AMR::as.disk()]
#' @details
#' The `compare_*_with_eucast()` functions allow to compare a user range with EUCAST distributions. Use [ggplot2::autoplot()] on the output to visualise the results.
#'
#' ### Supported Antimicrobials
#'
#' In December 2024, EUCAST had 176 distributions available, namely for these antimicrobials:
#'
#' Amikacin, amoxicillin, amoxicillin/clavulanic acid, amphotericin B, ampicillin, ampicillin/sulbactam, anidulafungin, apramycin, aspoxicillin, avilamycin, azithromycin, aztreonam, aztreonam/avibactam, bacitracin, bedaquiline, benzylpenicillin, capreomycin, cefaclor, cefadroxil, cefalexin, cefaloridine, cefalotin, cefapirin, cefazolin, cefdinir, cefepime, cefepime/tazobactam, cefepime/zidebactam, cefiderocol, cefixime, cefoperazone, cefoperazone/sulbactam, cefoselis, cefotaxime, cefotetan, cefovecin, cefoxitin, cefpirome, cefpodoxime, cefpodoxime/clavulanic acid, cefquinome, ceftaroline, ceftazidime, ceftazidime/avibactam, ceftibuten, ceftiofur, ceftobiprole, ceftolozane/tazobactam, ceftriaxone, cefuroxime, cephradine, chloramphenicol, chlortetracycline, ciprofloxacin, clarithromycin, clavulanic acid, clinafloxacin, clindamycin, clofazimine, cloxacillin, colistin, cycloserine, dalbavancin, danofloxacin, daptomycin, delafloxacin, delamanid, dicloxacillin, difloxacin, doripenem, doxycycline, enrofloxacin, eravacycline, ertapenem, erythromycin, ethambutol, ethionamide, faropenem, fidaxomicin, florfenicol, flucloxacillin, fluconazole, flucytosine, flumequine, fosfomycin, fusidic acid, gamithromycin, gatifloxacin, gemifloxacin, gentamicin, imipenem, imipenem/relebactam, isavuconazole, isoniazid, itraconazole, kanamycin, ketoconazole, lefamulin, levofloxacin, lincomycin, linezolid, loracarbef, marbofloxacin, mecillinam, meropenem, meropenem/vaborbactam, metronidazole, micafungin, minocycline, moxifloxacin, mupirocin, nalidixic acid, narasin, neomycin, netilmicin, nitrofurantoin, nitroxoline, norfloxacin, norvancomycin, ofloxacin, omadacycline, orbifloxacin, oritavancin, oxacillin, oxolinic acid, oxytetracycline, pefloxacin, phenoxymethylpenicillin, piperacillin, piperacillin/tazobactam, pirlimycin, posaconazole, pradofloxacin, pristinamycin, pyrazinamide, quinupristin/dalfopristin, retapamulin, rezafungin, rifabutin, rifampicin, roxithromycin, secnidazole, sitafloxacin, spectinomycin, spiramycin, streptomycin, sulbactam, sulfadiazine, sulfamethoxazole, sulfisoxazole, tedizolid, teicoplanin, telavancin, telithromycin, temocillin, terbinafine, tetracycline, thiamphenicol, tiamulin, ticarcillin, ticarcillin/clavulanic acid, tigecycline, tildipirosin, tilmicosin, tobramycin, trimethoprim, trimethoprim/sulfamethoxazole, tulathromycin, tylosin, tylvalosin, vancomycin, viomycin, and voriconazole.
#'
#' For the current list, run [eucast_supported_ab_distributions()].
#' @importFrom rvest read_html html_element html_table
#' @importFrom AMR as.ab ab_name ab_atc as.mo as.mic as.disk mo_name
#' @importFrom dplyr mutate filter select matches %>%
#' @importFrom tidyr pivot_longer
#' @rdname get_eucast_amr_distribution
#' @export
#' @examples
#' get_eucast_mic_distribution("cipro")
#'
#' # not returning as frequency table
#' get_eucast_mic_distribution("cipro", as_freq_table = FALSE)
#'
#' # specify microorganism to only get results for that pathogen
#' get_eucast_mic_distribution("cipro", "K. pneumoniae")
#'
#' get_eucast_disk_distribution("cipro", "K. pneumoniae")
#'
#'
#' # Plotting ----------------------------------------------------------------
#'
#' mic_data <- get_eucast_mic_distribution("cipro", "K. pneumoniae")
#' mics <- rep(mic_data$mic, mic_data$count)
#' ggplot2::autoplot(mics, ab = "cipro", mo = "K. pneumoniae", title = "Look at my MICs!")
#'
#' disk_data <- get_eucast_disk_distribution("cipro", "K. pneumoniae")
#' disks <- rep(disk_data$disk_diffusion, disk_data$count)
#' ggplot2::autoplot(disks, ab = "cipro", mo = "K. pneumoniae", title = "Look at my diffusion zones!")
#'
#'
#' # Comparing With User Values ----------------------------------------------
#'
#' my_mic_values <- AMR::random_mic(500)
#' comparison <- compare_mic_with_eucast(my_mic_values, ab = "cipro", mo = "K. pneumoniae")
#' comparison
#' ggplot2::autoplot(comparison)
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
    if (interactive()) {
      message(
        "Returning antimicrobial wild type distributions for ",
        ab_name(ab_coerced, language = NULL, tolower = TRUE), " (", ab_coerced, ", ", ab_atc(ab_coerced, only_first = TRUE), ")",
        ifelse(!is.null(mo_coerced), paste0(" in ", font_italic(suppressWarnings(mo_name(mo_coerced, language = NULL, keep_synonyms = TRUE)))), "")
      )
    }
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
    mutate(
      microorganism_code = suppressWarnings(as.mo(microorganism, language = NULL, keep_synonyms = TRUE, info = FALSE)),
      .after = 1
    )

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
        values_to = "count"
      )
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
#' @importFrom AMR is.mic as.mic
#' @importFrom dplyr select filter as_tibble
#' @export
compare_mic_with_eucast <- function(mics, ab, mo = NULL) {
  if (!is.mic(mics)) {
    mics <- as.mic(mics)
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
#' @importFrom AMR is.disk as.disk
#' @importFrom dplyr select filter as_tibble
#' @export
compare_disk_with_eucast <- function(disks, ab, mo = NULL) {
  if (!is.disk(disks)) {
    disks <- as.disk(disks)
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
# @importFrom AMR scale_x_mic
#' @method autoplot compare_eucast
#' @export
autoplot.compare_eucast <- function(object, ...) {
  long <- object |>
    mutate(
      User = user / sum(user),
      EUCAST = eucast / sum(eucast)
    ) |>
    select(-user, -eucast) |>
    pivot_longer(-value, names_to = "Source", values_to = "count")
  ggplot(long, aes(x = value, y = count, fill = Source)) +
    geom_col(position = "dodge") +
    labs(x = "Measurement Value", y = "Density") # +
  # scale_x_mic(keep_operators = "none")
}
