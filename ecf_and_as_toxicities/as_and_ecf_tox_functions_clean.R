# =============================================================================
# as_and_ecf_tox_functions_clean.R
# Functions for analysing sigma and anti-sigma factor toxicity
# Author: Rebecca Wolters
#
# Experimental context:
#   This experiment measures the effect of expressing sigma (ECF) and
#   anti-sigma (AS) factors on bacterial growth and luminescence.
#   Key metadata columns:
#     antisigma   — anti-sigma factor name (e.g. "AS14_1324", "AS26_4464")
#     tox         — toxicity level / AS version denoting by numeric code (numeric code, e.g. 0, 1, 3.1, 3.2)
#     copies      — plasmid copy number
#     aTc         — aTc inducer concentration (ng/ml), drives AS expression
#     arabinose   — arabinose concentration (%), drives ECF/sigma expression
#     as_origin   — simplified AS name extracted from antisigma (e.g. "AS14")
# =============================================================================

# --- Shared plot theme -------------------------------------------------------

theme_tox <- function(base_size = 10) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(-0.15, "cm")
    )
}

# Colour scale shared across growth curve and bar plots
# Order: no-aTc/no-arab, aTc-only, arab-only, both inducers
tox_colors <- c("grey10", "cornflowerblue", "mediumorchid", "orchid4")


# =============================================================================
# DATA I/O
# =============================================================================

#' Read raw plate reader data files from a folder
#'
#' @param folder Path to folder containing raw .csv files
#' @param pattern File name pattern to match (default "tox.csv")
read_in_plates <- function(folder, pattern = "tox.csv") {
  list.files(folder, pattern = pattern, full.names = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(file = value) |>
    dplyr::mutate(
      data = purrr::map(file, wellr::plate_read_biotek),
      exp  = dplyr::row_number()
    ) |>
    dplyr::select(-file) |>
    tidyr::unnest(data) |>
    dplyr::mutate(time = time / 60 / 60)   # seconds -> hours
}

#' Read metainfo files from a folder (plater format)
#'
#' @param folder Path to folder containing metainfo .csv files
#' @param pattern File name pattern to match
read_in_meta <- function(folder, pattern = "metainfo.csv") {
  list.files(folder, pattern = pattern, full.names = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(file = value) |>
    dplyr::mutate(
      data = purrr::map(file, plater::read_plate),
      exp  = dplyr::row_number()
    ) |>
    dplyr::select(-file) |>
    tidyr::unnest(data)
}


# =============================================================================
# BLANKING & QC
# =============================================================================

#' Plot blank OD readings before induction for QC
#'
#' @param .data Raw joined dataframe
#' @param blank_variable Column name to plot on y-axis (string)
#' @param preinduction_time Time of induction (hours)
#' @param variable_to_color Column name to colour by (string)
check_blank <- function(.data, blank_variable, preinduction_time, variable_to_color) {
  .data |>
    dplyr::filter(strain == "blank", time < preinduction_time) |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = .data[[blank_variable]],
      color = as.factor(.data[[variable_to_color]])
    )) +
    ggplot2::geom_point() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(color = variable_to_color) +
    theme_tox()
}

#' Get average blank OD value before induction
#'
#' @param .data Raw dataframe
#' @param preinduction_time Time of induction (hours)
get_blank_average <- function(.data, preinduction_time) {
  .data |>
    dplyr::filter(strain == "blank", time < preinduction_time) |>
    dplyr::summarise(blank_od600_mean = mean(od600, na.rm = TRUE))
}

#' Background correct OD600 and calculate lum/OD
#'
#' @param .data Raw joined dataframe
#' @param average_blank_od Average blank OD value (from get_blank_average)
#' @param pc_factor Pathlength correction factor (96-well: 4.898, 384-well: 3.09)
blank_plate <- function(.data, average_blank_od, pc_factor) {
  .data |>
    tidyr::drop_na(strain) |>
    dplyr::mutate(
      od600_blank = od600 * pc_factor - average_blank_od,
      lum_blank   = lum / od600_blank
    )
}


# =============================================================================
# AVERAGING
# =============================================================================

#' Average blanked data across replicates
#'
#' @param .data Blanked dataframe (df_blank)
#' @param group_vars Character vector of grouping columns (must include "time")
average_plate <- function(.data, group_vars) {
  .data |>
    dplyr::summarise(
      dplyr::across(
        dplyr::matches("blank"),
        list(mean = mean, sd = sd),
        .names = "{.col}_{.fn}"
      ),
      replicate = paste(well, collapse = " "),
      .by       = all_of(group_vars)
    )
}

#' Average growth rate and doubling time data across replicates
#'
#' @param .data Growth rate dataframe (df_gr)
#' @param group_vars Character vector of grouping columns
average_gr_data <- function(.data, group_vars) {
  .data |>
    dplyr::summarise(
      dplyr::across(
        dplyr::matches("gr|dt"),
        list(mean = mean, sd = sd),
        .names = "{.col}_{.fn}"
      ),
      .by = all_of(group_vars)
    )
}

#' Average toxicity (% OD difference from WT) across replicates
#'
#' @param .data Dataframe with toxicity_blank column
#' @param group_vars Character vector of grouping columns
average_tox_data <- function(.data, group_vars) {
  .data |>
    dplyr::summarise(
      dplyr::across(
        dplyr::matches("toxicity_blank"),
        list(mean = mean, sd = sd),
        .names = "{.col}_{.fn}"
      ),
      .by = all_of(group_vars)
    )
}


# =============================================================================
# GROWTH RATE
# =============================================================================

#' Fit linear model to log(OD600) to extract growth rate and doubling time
#'
#' Fits a linear model to log-transformed OD600 within a defined time window.
#' Growth rate (gr) = slope of log(OD) vs time.
#' Doubling time (dt) = log(2) / gr.
#'
#' @param .data Blanked dataframe (df_blank)
#' @param group_vars Character vector of grouping columns (must include rep, exp)
#' @param starttime Start of exponential growth window (hours)
#' @param endtime End of exponential growth window (hours)
get_growthrate <- function(.data, group_vars, starttime, endtime) {
  .data |>
    dplyr::filter(
      dplyr::between(time, starttime, endtime),
      strain != "blank"
    ) |>
    dplyr::group_by(dplyr::across(all_of(group_vars))) |>
    dplyr::mutate(od = log(od600_blank)) |>
    dplyr::group_modify(~ broom::tidy(lm(od ~ time, data = .x))) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      gr = estimate,
      dt = log(2) / gr
    ) |>
    dplyr::ungroup()
}


# =============================================================================
# TOXICITY CALCULATION
# =============================================================================

#' Calculate relative toxicity as % OD difference compared to wildtype
#'
#' For each strain, calculates how much its OD600 deviates from the wildtype
#' strain at the same condition. Positive = lower OD than WT (toxic);
#' negative = higher OD than WT.
#'
#' Formula: (1 - gmo_od / wt_od) * 100
#'
#' @param .data Averaged dataframe (df_av) filtered to time == max(time)
#' @param wt_strain Strain ID of the wildtype reference (string, e.g. "SV01")
#' @param aTc_concentration Numeric vector of aTc concentrations to include
#' @param basedon_gr Logical — if TRUE use growth rate; if FALSE use endpoint OD
get_relative_tox <- function(.data, wt_strain, aTc_concentration, basedon_gr = FALSE) {

  percent_diff <- function(wt_od, gmo_od) {
    (1 - gmo_od / wt_od) * 100
  }

  .data |>
    dplyr::filter(aTc %in% aTc_concentration) |>
    dplyr::mutate(
      as_origin      = stringr::str_extract(antisigma, "([A-Z]{2}\\d{2})"),
      toxicity_blank = percent_diff(
        wt_od  = od600_blank_mean[strain == wt_strain],
        gmo_od = od600_blank_mean
      )
    ) |>
    dplyr::filter(strain != wt_strain)
}

#' Calculate luminescence fold-repression between aTc-induced and uninduced
#'
#' Pivots data wide on aTc concentration and computes fold-repression
#' (lum without aTc / lum with aTc). Uses averaged data (df_av).
#' Filters to max arabinose and final time point by default.
#'
#' @param .data Averaged dataframe (df_av)
#' @param atc_uninduced aTc concentration for uninduced condition (default 0)
#' @param atc_induced aTc concentration for induced condition (default 100)
get_fold_repression <- function(.data, atc_uninduced = 0, atc_induced = 100) {
  .data |>
    dplyr::filter(
      time     == max(time),
      arabinose == max(arabinose),
      !antisigma %in% c("blank", "none", "Ptet-lux", "wildtype")
    ) |>
    dplyr::mutate(
      as_origin = stringr::str_extract(antisigma, "([A-Z]{2}\\d{2})")
    ) |>
    dplyr::select(aTc, antisigma, lum_blank_mean, lum_blank_sd, tox, as_origin) |>
    tidyr::pivot_wider(
      names_from  = "aTc",
      values_from = dplyr::matches("lum")
    ) |>
    dplyr::mutate(
      fc_mean = .data[[paste0("lum_blank_mean_", atc_uninduced)]] /
                .data[[paste0("lum_blank_mean_", atc_induced)]],
      fc_sd   = sqrt(
        (.data[[paste0("lum_blank_sd_", atc_uninduced)]]  /
         .data[[paste0("lum_blank_mean_", atc_uninduced)]])^2 +
        (.data[[paste0("lum_blank_sd_", atc_induced)]]    /
         .data[[paste0("lum_blank_mean_", atc_induced)]])^2
      ) * fc_mean
    )
}


# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

#' Run ANOVA + Tukey HSD on growth rate across inducer combinations per strain
#'
#' Tests whether growth rate differs significantly across the four inducer
#' combinations (aTc × arabinose) within each strain.
#' Returns a tidy table of pairwise comparisons with adjusted p-values.
#'
#' @param .data Raw growth rate dataframe (df_gr), unaveraged
#' @param antisigma_exclude Character vector of antisigma values to exclude
#' @param tox_levels Numeric vector of tox levels to include
run_anova_tukey <- function(.data,
                             antisigma_exclude = c("none", "blank", "Ptet-lux"),
                             tox_levels        = c(0, 3.1, 3.2)) {
  df_for_anova <- .data |>
    dplyr::filter(
      !antisigma %in% antisigma_exclude,
      tox        %in% tox_levels
    ) |>
    dplyr::select(strain, antisigma, aTc, arabinose, copies, tox, gr) |>
    dplyr::mutate(
      aTc            = as.factor(aTc),
      arabinose      = as.factor(arabinose),
      interaction_var = interaction(aTc, arabinose)
    )

  aov_model <- aov(gr ~ interaction_var * strain, data = df_for_anova)
  tukey     <- rstatix::tukey_hsd(aov_model)

  tukey |>
    dplyr::select(term, group1, group2, estimate, p.adj, p.adj.signif)
}

#' Run paired t-tests comparing induced vs uninduced growth rate per antisigma
#'
#' For each antisigma factor, tests whether aTc induction significantly
#' changes growth rate compared to the uninduced condition.
#' Returns a tidy dataframe with t-statistics, p-values, and significance labels.
#'
#' @param .data Averaged growth rate dataframe (df_gr_av)
#' @param antisigma_exclude Character vector of antisigma values to exclude
#' @param tox_levels Numeric vector of tox levels to include
#' @param atc_low aTc concentration for uninduced condition (default 0)
#' @param atc_high aTc concentration for induced condition (default 100)
run_t_tests <- function(.data,
                         antisigma_exclude = c("none", "blank", "Ptet-lux"),
                         tox_levels        = c(0, 3.1, 3.2),
                         atc_low           = 0,
                         atc_high          = 100) {
  .data |>
    dplyr::filter(
      !antisigma %in% antisigma_exclude,
      tox        %in% tox_levels,
      aTc        %in% c(atc_low, atc_high)
    ) |>
    dplyr::reframe(
      {
        low_vals  <- gr_mean[aTc == atc_low]
        high_vals <- gr_mean[aTc == atc_high]

        if (length(low_vals) < 2 || length(high_vals) < 2) {
          data.frame(
            t_statistic = NA_real_, df = NA_real_,
            p_value = NA_real_, mean_diff = NA_real_
          )
        } else {
          result <- t.test(low_vals, high_vals,
                           paired = length(low_vals) == length(high_vals))
          data.frame(
            t_statistic = result$statistic,
            df          = result$parameter,
            p_value     = result$p.value,
            mean_diff   = result$estimate[[1]]
          )
        }
      },
      .by = c(antisigma, tox, copies)
    ) |>
    dplyr::mutate(
      significance = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        p_value >= 0.05 ~ "ns",
        is.na(p_value)  ~ NA_character_
      )
    )
}


# =============================================================================
# PLOTTING — QC
# =============================================================================

#' Plot data in plate layout format for QC
#'
#' @param .data Blanked dataframe (df_blank)
#' @param variable_to_plot Bare column name to plot on y-axis
#' @param variable_to_color Bare column name to colour by
#' @param number_of_experiment Experiment number to display (integer)
#' @param number_of_colors Number of colours for palette
plot_plate <- function(.data,
                       variable_to_plot,
                       variable_to_color,
                       number_of_experiment,
                       number_of_colors) {
  .data |>
    dplyr::filter(exp == number_of_experiment) |>
    dplyr::mutate(
      col = stringr::str_extract(well, "\\d+$"),
      row = factor(stringr::str_extract(well, "^\\w"), levels = LETTERS)
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = {{ variable_to_plot }},
      color = as.factor({{ variable_to_color }})
    )) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(col), rows = ggplot2::vars(row)) +
    ggplot2::geom_vline(xintercept = 2, linetype = "solid", colour = "grey") +
    ggplot2::scale_color_manual(
      values = MetBrewer::met.brewer("Renoir", number_of_colors, type = "continuous")
    ) +
    ggplot2::scale_y_log10(minor_breaks = minor_breaks) +
    theme_tox()
}


# =============================================================================
# PLOTTING — GROWTH CURVES
# =============================================================================

#' Plot OD600 kinetics faceted by antisigma factor
#'
#' Coloured by the interaction of arabinose × aTc concentrations, so all four
#' inducer combinations are visually distinguishable.
#'
#' @param .data Averaged dataframe (df_av)
#' @param tox_level Numeric vector of tox levels to include
plot_growthcurves <- function(.data, tox_level) {
  .data |>
    dplyr::filter(tox %in% tox_level) |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = as.numeric(od600_blank_mean),
      color = factor(interaction(arabinose, aTc))
    )) +
    ggplot2::geom_point(size = 0.25) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin   = pmax(od600_blank_mean - od600_blank_sd, 0.01),
        ymax   = od600_blank_mean + od600_blank_sd,
        colour = factor(interaction(arabinose, aTc))
      ),
      alpha = 0.6,
      linewidth = 0.2
    ) +
    ggplot2::facet_wrap(facets = "antisigma") +
    ggplot2::geom_vline(xintercept = 2, linetype = "solid", colour = "grey") +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0, 10)
    ) +
    ggplot2::scale_y_log10(labels = scales::label_number()) +
    ggplot2::scale_color_manual(values = tox_colors) +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("OD"[600]),
      color = "arabinose (%) + aTc (ng/ml)"
    ) +
    theme_tox() +
    ggplot2::annotation_logticks(sides = "l")
}

#' Plot OD600 kinetics for all strains overlaid in a single panel
#'
#' Useful for comparing growth across all antisigma factors simultaneously.
#' Add facet_grid() after the call to split by inducer condition.
#'
#' @param .data Averaged dataframe (df_av), pre-filtered to strains of interest
plot_growthcurves_allinone <- function(.data) {
  .data |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = as.numeric(od600_blank_mean),
      color = factor(interaction(arabinose, aTc))
    )) +
    ggplot2::geom_point(size = 0.25) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin   = pmax(od600_blank_mean - od600_blank_sd, 0.01),
        ymax   = od600_blank_mean + od600_blank_sd,
        colour = factor(interaction(arabinose, aTc))
      ),
      alpha = 0.6,
      linewidth = 0.2
    ) +
    ggplot2::geom_vline(xintercept = 2, linetype = "solid", colour = "grey") +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0, 10)
    ) +
    ggplot2::scale_y_log10(labels = scales::label_number()) +
    ggplot2::scale_color_manual(values = tox_colors) +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("OD"[600]),
      color = "arabinose (%) + aTc (ng/ml)"
    ) +
    theme_tox() +
    ggplot2::annotation_logticks(sides = "l")
}


# =============================================================================
# PLOTTING — GROWTH RATE BAR PLOTS
# =============================================================================

#' Plot growth rate or doubling time as faceted bar plot (faceted by antisigma × arabinose)
#'
#' @param .data Averaged growth rate dataframe (df_gr_av)
#' @param tox_level Numeric vector of tox levels to include
#' @param variable_to_plot Bare column name of the y-axis variable (e.g. dt_mean)
#' @param variable_to_plot_sd Bare column name of the SD (e.g. dt_sd)
#' @param y_axis_label String for the y-axis label
plot_bar_faceted <- function(.data,
                              tox_level,
                              variable_to_plot,
                              variable_to_plot_sd,
                              y_axis_label) {
  .data |>
    dplyr::filter(tox %in% tox_level) |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(aTc),
      y    = {{ variable_to_plot }},
      fill = as.factor(aTc)
    )) +
    ggplot2::geom_col(colour = "black", position = "dodge") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = {{ variable_to_plot }} - {{ variable_to_plot_sd }},
        ymax = {{ variable_to_plot }} + {{ variable_to_plot_sd }}
      ),
      width    = 0.5,
      colour   = "black",
      alpha    = 0.8
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(antisigma),
      rows = ggplot2::vars(arabinose)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +
    ggplot2::scale_fill_manual(values = c("white", "grey26")) +
    ggplot2::labs(
      x    = "aTc (ng/ml)",
      y    = y_axis_label,
      fill = "aTc (ng/ml)"
    ) +
    theme_tox() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#' Plot growth rate or doubling time as grouped bar plot (grouped by antisigma)
#'
#' Bars are grouped by the interaction of aTc × arabinose concentrations.
#' A dashed reference line marks the wildtype growth rate (default 0.25 doublings/h).
#'
#' @param .data Averaged growth rate dataframe (df_gr_av) with interaction_var column
#' @param tox_level Numeric vector of tox levels to include
#' @param variable_to_plot Bare column name of y-axis variable (e.g. gr_mean)
#' @param variable_to_plot_sd Bare column name of SD (e.g. gr_sd)
#' @param y_axis_label String for the y-axis label
#' @param wt_gr Reference growth rate for dashed line (default 0.25)
plot_bar_grouped <- function(.data,
                              tox_level,
                              variable_to_plot,
                              variable_to_plot_sd,
                              y_axis_label,
                              wt_gr = 0.25) {
  .data |>
    dplyr::filter(tox %in% tox_level) |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(antisigma),
      y    = {{ variable_to_plot }},
      fill = interaction_var
    )) +
    ggplot2::geom_col(
      stat     = "identity",
      width    = 0.5,
      position = ggplot2::position_dodge(width = 0.5),
      colour   = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = {{ variable_to_plot }} - {{ variable_to_plot_sd }},
        ymax = {{ variable_to_plot }} + {{ variable_to_plot_sd }}
      ),
      width    = 0.2,
      position = ggplot2::position_dodge(width = 0.5),
      colour   = "black"
    ) +
    # Dashed reference line at wildtype growth rate
    ggplot2::geom_hline(
      yintercept = wt_gr,
      linetype   = "dashed",
      colour     = "brown"
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(-0.05, 0.45)
    ) +
    ggplot2::scale_fill_manual(values = tox_colors) +
    ggplot2::labs(
      x    = expression("anti-" * sigma),
      y    = y_axis_label,
      fill = "arabinose (%) + aTc (ng/ml)"
    ) +
    theme_tox() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )
}


# =============================================================================
# PLOTTING — ENDPOINT & KINETIC PLOTS
# =============================================================================

#' Plot endpoint OD600 (final time point) as bar plot
#'
#' @param .data Averaged dataframe (df_av)
#' @param tox_level Numeric vector of tox levels to include
plot_endpoint_od <- function(.data, tox_level) {
  .data |>
    dplyr::filter(
      tox      %in% tox_level,
      time     == max(time),
      !antisigma %in% c("none", "blank")
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(aTc),
      y    = od600_blank_mean,
      fill = as.factor(aTc)
    )) +
    ggplot2::geom_col(
      width    = 0.5,
      position = ggplot2::position_dodge(width = 0.5),
      colour   = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = od600_blank_mean - od600_blank_sd,
        ymax = od600_blank_mean + od600_blank_sd
      ),
      width    = 0.2,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(antisigma),
      rows = ggplot2::vars(arabinose)
    ) +
    ggplot2::scale_fill_manual(values = c("white", "grey26")) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(0.1, 1.2)) +
    ggplot2::labs(
      x    = expression("anti-" * sigma),
      y    = expression("OD"[600]),
      fill = "aTc (ng/ml)"
    ) +
    theme_tox() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#' Plot luminescence at induction time point as bar plot
#'
#' Filters to the time window just after induction (between 6.002 and 6.09 h)
#' to capture the luminescence level at the point of AS induction.
#'
#' @param .data Averaged dataframe (df_av)
#' @param tox_level Numeric vector of tox levels to include
plot_endpoint_lum <- function(.data, tox_level) {
  .data |>
    dplyr::filter(
      tox %in% tox_level,
      dplyr::between(time, 6.002, 6.09)
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(antisigma),
      y    = lum_blank_mean,
      fill = as.factor(aTc)
    )) +
    ggplot2::geom_col(
      width    = 0.5,
      position = ggplot2::position_dodge(width = 0.5),
      colour   = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(lum_blank_mean - lum_blank_sd, 1),
        ymax = lum_blank_mean + lum_blank_sd
      ),
      width    = 0.2,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(arabinose)) +
    ggplot2::scale_fill_manual(values = c("goldenrod1", "white")) +
    ggplot2::scale_y_log10(
      labels       = scales::label_log(),
      minor_breaks = minor_breaks,
      expand       = c(0, 0)
    ) +
    ggplot2::labs(
      x    = expression("anti-" * sigma),
      y    = expression("RLU/OD"[600]),
      fill = "aTc (ng/ml)"
    ) +
    theme_tox() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}

#' Plot luminescence kinetics over time faceted by antisigma factor
#'
#' @param .data Averaged dataframe (df_av)
#' @param tox_level Numeric vector of tox levels to include
plot_kinetic_lum <- function(.data, tox_level) {
  .data |>
    dplyr::filter(tox %in% tox_level) |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = as.numeric(lum_blank_mean),
      color = as.factor(interaction(aTc, arabinose))
    )) +
    ggplot2::geom_point(size = 0.25) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin  = pmax(lum_blank_mean - lum_blank_sd, 10),
        ymax  = lum_blank_mean + lum_blank_sd,
        alpha = 0.5
      ),
      linewidth = 0.1
    ) +
    ggplot2::facet_wrap(facets = ggplot2::vars(antisigma)) +
    ggplot2::geom_vline(xintercept = 2, linetype = "dashed", colour = "grey") +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion()
    ) +
    ggplot2::scale_y_log10(
      labels       = scales::label_log(),
      minor_breaks = minor_breaks,
      limits       = c(1, 1000000)
    ) +
    ggplot2::scale_color_manual(values = tox_colors) +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("RLU/OD"[600]),
      color = "aTc (ng/ml) + arabinose (%)"
    ) +
    theme_tox() +
    ggplot2::annotation_logticks(sides = "l")
}

