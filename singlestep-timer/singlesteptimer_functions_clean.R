# =============================================================================
# singlesteptimer_functions.R
# Functions for single-step timer luminescence plate reader analysis
# Author: Rebecca Wolters
# =============================================================================

# --- Shared plot theme -------------------------------------------------------
# Defined once here and reused across all plot functions.
# Override individual elements by adding + theme(...)

theme_timer <- function(base_size = 8) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position   = "bottom",
      panel.grid.major  = ggplot2::element_line(color = "transparent"),
      panel.grid.minor  = ggplot2::element_line(color = "transparent"),
      axis.ticks.length = ggplot2::unit(-0.1, "cm"),
      panel.spacing.x   = ggplot2::unit(0.5, "cm")
    )
}

# Shared log y-axis scale for luminescence plots
scale_y_lum <- function(limits = c(10, 300000), minor_breaks = NULL) {
  ggplot2::scale_y_log10(
    labels       = function(x) parse(text = paste0("10^", log10(x))),
    minor_breaks = minor_breaks,
    limits       = limits
  )
}

# Shared log y-axis scale for fold-change / dynamic range plots
scale_y_fc <- function(limits = c(1, 10000)) {
  ggplot2::scale_y_log10(
    labels = function(x) parse(text = paste0("10^", log10(x))),
    limits = limits
  )
}

# Shared logticks annotation
logticks_left <- function() {
  ggplot2::annotation_logticks(
    sides = "l",
    short = ggplot2::unit(0.025, "cm"),
    mid   = ggplot2::unit(0.05,  "cm"),
    long  = ggplot2::unit(0.1,   "cm")
  )
}


# =============================================================================
# DATA I/O
# =============================================================================

#' Read raw plate reader data files from a folder
#'
#' @param folder Path to folder containing raw .csv files
#' @param pattern File name pattern to match (default "RW.csv")
#' @return Tibble with raw plate reader data, time converted from seconds to hours
read_in_plates <- function(folder, pattern = "RW.csv") {
  
  list.files(folder, pattern = pattern, full.names = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(file = value) |>
    dplyr::mutate(
      data = purrr::map(file, wellr::plate_read_biotek),
      exp  = dplyr::row_number()
    ) |>
    dplyr::select(-file) |>
    tidyr::unnest(data) |>
    dplyr::mutate(time = time / 60 / 60)  # seconds -> hours
}

#' Read metainfo files from a folder (plater format)
#'
#' @param folder Path to folder containing metainfo .csv files
#' @param pattern File name pattern to match (default "metainfo.csv")
#' @return Tibble with plate layout metadata
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

#' Plot blank readings before induction for QC
#'
#' @param .data Raw joined dataframe
#' @param blank_variable Column name of the variable to plot (string)
#' @param preinduction_time Time point of induction (hours)
#' @param variable_to_color Column name to use for colour (string)
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
    theme_timer()
}

#' Get average blank value before induction
#'
#' @param .data Raw dataframe
#' @param preinduction_time Time point of induction (hours)
#' @param blank_lum_or_od600 Bare column name to average (e.g. od600)
get_blank_average <- function(.data, preinduction_time, blank_lum_or_od600) {
  
  .data |>
    dplyr::filter(strain == "blank", time < preinduction_time) |>
    dplyr::summarise(blank_mean = mean({{ blank_lum_or_od600 }}, na.rm = TRUE))
}

#' Background correct OD600 and calculate lum/OD when NO blank well on plate
#'
#' @param data Raw joined dataframe
#' @param average_blank_od Average blank OD value (from get_blank_average)
#' @param pc_factor Pathlength correction factor (96-well: 4.898, 384-well: 3.09)
blank_plate <- function(data, average_blank_od, pc_factor) {
  
  data |>
    tidyr::drop_na(strain) |>
    dplyr::mutate(
      od600_blank = od600 * pc_factor - average_blank_od,
      lum_blank   = lum / od600_blank
    )
}

#' Background correct OD600 when blank well IS included on the plate
#'
#' @param .data Raw joined dataframe
#' @param pc_factor Pathlength correction factor
#' @param group_vars Character vector of columns to group by for blank subtraction
blank_on_plate <- function(.data, pc_factor, group_vars) {
  
  .data |>
    tidyr::drop_na(strain) |>
    dplyr::mutate(
      od600_blank = (od600 - mean(od600[strain == "blank"])) * pc_factor,
      lum_blank   = lum / od600_blank,
      .by         = all_of(group_vars)
    ) |>
    dplyr::filter(strain != "blank")
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

#' Average growth rate data across replicates
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


# =============================================================================
# BACKGROUND LUMINESCENCE
# =============================================================================

#' Calculate pre-induction background luminescence per condition
#'
#' Uses the pre-induction window (starttime to inductiontime) to estimate
#' baseline luminescence. Set wildtype = TRUE to use only the SV01 strain
#' as the global background reference; FALSE (default) uses each strain's own baseline.
#'
#' @param .data Averaged dataframe (df_av)
#' @param value Bare column name of luminescence variable (e.g. lum_blank_mean)
#' @param starttime Start of baseline window (hours)
#' @param inductiontime End of baseline window / induction time (hours)
#' @param group_vars Character vector of grouping columns
#' @param wildtype Logical. TRUE = use SV01 as reference; FALSE = per-strain baseline
get_background_lum <- function(.data,
                                value,
                                starttime,
                                inductiontime,
                                group_vars,
                                wildtype = FALSE) {
  
  df_pre <- .data |>
    dplyr::filter(dplyr::between(time, starttime, inductiontime)) |>
    dplyr::filter(!is.na({{ value }}))
  
  if (wildtype) {
    df_pre <- df_pre |> dplyr::filter(strain == "SV01")
  }
  
  df_pre |>
    dplyr::summarise(
      background    = mean({{ value }}, na.rm = TRUE),
      background_sd = sd({{ value }},   na.rm = TRUE),
      .by           = all_of(group_vars)
    )
}


# =============================================================================
# TIME DELAY
# =============================================================================

#' Calculate time delay using discrete measured time points (original method)
#'
#' Returns the first measured time point after induction where fold-change
#' exceeds the threshold. Note: resolution is limited to the measurement interval.
#' Use get_timedelay_splinefun() for interpolated estimates.
#'
#' @param .data Averaged dataframe (df_av)
#' @param bg_expr Background expression table from get_background_lum()
#' @param group_vars Character vector of grouping columns
#' @param foldchange Fold-change threshold (default 4)
#' @param inductiontime Time of induction (hours)
get_timedelay <- function(.data, bg_expr, group_vars, foldchange, inductiontime) {
  
  .data |>
    dplyr::left_join(bg_expr, by = group_vars) |>
    dplyr::mutate(
      fc_mean = lum_blank_mean / background,
      fc_sd   = sqrt((lum_blank_sd / lum_blank_mean)^2 +
                     (background_sd / background)^2) * fc_mean,
      .by     = group_vars
    ) |>
    dplyr::filter(time > inductiontime) |>
    dplyr::summarise(
      delay = dplyr::if_else(
        any(fc_mean > foldchange),
        time[which.max(fc_mean > foldchange)],
        NA_real_
      ),
      .by = group_vars
    )
}

#' Calculate time delay using cubic spline interpolation (recommended method)
#'
#' Fits a monotone cubic spline (Fritsch-Carlson) through post-induction
#' luminescence values and uses Brent's root-finding method (uniroot) to
#' determine the exact time at which fold-change crosses the threshold.
#' This is gives a more accurate evaluation than using the discrete measurement intervals.
#'
#' @param .data Averaged dataframe (df_av)
#' @param bg_expr Background expression table from get_background_lum()
#' @param group_vars Character vector of grouping columns
#' @param foldchange Fold-change threshold (default 4)
#' @param inductiontime Time of induction (hours)
get_timedelay_splinefun <- function(.data, bg_expr, group_vars, foldchange, inductiontime) {
  
  .data |>
    dplyr::left_join(bg_expr, by = group_vars) |>
    dplyr::mutate(
      fc_mean = lum_blank_mean / background,
      fc_sd   = sqrt((lum_blank_sd / lum_blank_mean)^2 +
                     (background_sd / background)^2) * fc_mean,
      .by     = group_vars
    ) |>
    dplyr::filter(time > inductiontime) |>
    dplyr::reframe(
      delay = {
        t  <- time
        fc <- fc_mean
        
        if (length(t) < 3 || !any(fc > foldchange, na.rm = TRUE)) {
          NA_real_
        } else {
          spline_fn <- splinefun(t, fc, method = "monoH.FC")
          cross_idx <- which(fc[-length(fc)] <= foldchange & fc[-1] > foldchange)[1]
          
          if (is.na(cross_idx)) {
            NA_real_
          } else {
            uniroot(
              f        = function(t_val) spline_fn(t_val) - foldchange,
              interval = c(t[cross_idx], t[cross_idx + 1]),
              tol      = 1e-6
            )$root
          }
        }
      },
      .by = group_vars
    )
}


# =============================================================================
# DYNAMIC RANGE / FOLD-CHANGE
# =============================================================================

#' Calculate fold-change over background at the final time point per condition
#'
#' Joins background values, calculates per-replicate fold-change, then
#' summarises across replicates. Negative fold-changes are floored at 1
#' (= no induction above baseline).
#'
#' @param .data Blanked dataframe (df_blank)
#' @param bg_expr Background expression table from get_background_lum()
#' @param group_vars Character vector of grouping columns
#' @param endpoint_hours Time point to use as endpoint (hours). If NULL, uses max(time).
average_dynamic_range <- function(.data, bg_expr, group_vars, endpoint_hours = NULL) {
  
  df <- .data |>
    dplyr::left_join(bg_expr, by = group_vars) |>
    dplyr::mutate(
      fold_change = pmax(lum_blank / background, 1),  # floor at 1 (no negative FC)
      .by         = all_of(c(group_vars, "rep", "exp"))
    )
  
  # Filter to endpoint
  if (!is.null(endpoint_hours)) {
    df <- df |>
      dplyr::filter(
        abs(time - endpoint_hours) == min(abs(time - endpoint_hours)),
        .by = all_of(c(group_vars, "rep", "exp"))
      )
  } else {
    df <- df |>
      dplyr::filter(
        time == max(time),
        .by = all_of(c(group_vars, "rep", "exp"))
      )
  }
  
  df |>
    dplyr::summarise(
      fold_change_mean = mean(fold_change, na.rm = TRUE),
      fold_change_sd   = sd(fold_change,   na.rm = TRUE),
      n                = dplyr::n(),
      replicate        = paste(well, collapse = " "),
      .by              = all_of(group_vars)
    )
}

#' Calculate fold-change at all time points (for heatmaps / delay vs FC plots)
#'
#' @param .data Averaged dataframe (df_av)
#' @param bg_expr Background expression table from get_background_lum()
#' @param group_vars Character vector of grouping columns
get_dynamic_range_all <- function(.data, bg_expr, group_vars) {
  
  .data |>
    dplyr::left_join(bg_expr, by = group_vars) |>
    dplyr::mutate(
      fc_mean = pmax(lum_blank_mean / background, 1),
      fc_sd   = sqrt((lum_blank_sd / lum_blank_mean)^2 +
                     (background_sd / background)^2) * fc_mean,
      .by     = group_vars
    )
}


# =============================================================================
# SIGNIFICANCE TESTING
# =============================================================================

#' Test whether threshold induction significantly changes fold-change per strain
#'
#' Runs a paired t-test (on log10 scale) per strain comparing maximum fold-change
#' at threshold ON (DHBA > 0 or aTc > 0) vs threshold OFF (both = 0),
#' collapsing across arabinose concentrations.
#'
#' @param .data Blanked dataframe (df_blank)
#' @param bg_expr Background expression table from get_background_lum()
#' @param group_vars Character vector of grouping columns
#' @param strain_group Character vector of strain IDs to test
#' @param threshold_inducer String: "aTc" or "DHBA"
test_threshold_significance <- function(.data,
                                         bg_expr,
                                         group_vars,
                                         strain_group,
                                         threshold_inducer) {
  
  .data |>
    dplyr::filter(strain %in% strain_group) |>
    dplyr::left_join(bg_expr, by = group_vars) |>
    dplyr::mutate(
      fold_change = pmax(lum_blank / background, 1)
    ) |>
    dplyr::mutate(
      threshold = dplyr::if_else(.data[[threshold_inducer]] > 0, "ON", "OFF")
    ) |>
    dplyr::filter(
      fold_change == max(fold_change, na.rm = TRUE),
      .by = c(antisigma, strain, all_of(threshold_inducer), rep, exp)
    ) |>
    dplyr::mutate(log_fc = log10(fold_change)) |>
    dplyr::reframe(
      {
        on_vals  <- log_fc[threshold == "ON"]
        off_vals <- log_fc[threshold == "OFF"]
        
        if (length(on_vals) < 2 || length(off_vals) < 2) {
          data.frame(
            t_statistic = NA_real_,
            df          = NA_real_,
            p_value     = NA_real_,
            mean_diff   = NA_real_,
            fold_diff   = NA_real_
          )
        } else {
          result <- t.test(on_vals, off_vals,
                           paired = length(on_vals) == length(off_vals))
          data.frame(
            t_statistic = result$statistic,
            df          = result$parameter,
            p_value     = result$p.value,
            mean_diff   = result$estimate[[1]],
            fold_diff   = 10^result$estimate[[1]]
          )
        }
      },
      .by = c(strain, antisigma)
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
# GROWTH RATE
# =============================================================================

#' Fit linear model to log(OD600) to extract growth rate and doubling time
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
# PLOTTING — QC & OVERVIEW
# =============================================================================

#' Plot data in plate layout format for QC
#'
#' @param .data Blanked dataframe (df_blank)
#' @param variable_to_plot Column name to plot on y-axis (string)
#' @param variable_to_color Column name to colour by (string)
#' @param inductiontime Time of induction (hours)
#' @param number_of_experiment Experiment number to plot (integer)
#' @param number_of_colors Number of colours for palette
plot_plate <- function(.data,
                       variable_to_plot,
                       variable_to_color,
                       inductiontime,
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
      y     = .data[[variable_to_plot]],
      color = as.factor(.data[[variable_to_color]])
    )) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(col), rows = ggplot2::vars(row)) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid", colour = "grey") +
    ggplot2::scale_color_manual(
      values = MetBrewer::met.brewer("Renoir", number_of_colors, type = "continuous")
    ) +
    ggplot2::scale_y_log10(minor_breaks = minor_breaks) +
    ggplot2::labs(color = variable_to_color) +
    theme_timer()
}

#' Plot raw replicate luminescence kinetics for QC
#'
#' @param .data Blanked dataframe (df_blank), pre-filtered to strain of interest
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label for colour legend (string)
#' @param inductiontime Time of induction (hours)
#' @param background_lum Background luminescence value for reference line
#' @param grid_col Bare column name for facet columns
#' @param grid_row Bare column name for facet rows
plot_replicates_lum <- function(.data,
                                 variable_to_colour,
                                 unit_of_colour,
                                 inductiontime,
                                 background_lum,
                                 grid_col,
                                 grid_row) {
  
  .data |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = as.numeric(lum_blank),
      color = as.factor(.data[[variable_to_colour]])
    )) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::facet_grid(
      cols = ggplot2::vars({{ grid_col }}),
      rows = ggplot2::vars({{ grid_row }})
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid", colour = "grey") +
    ggplot2::geom_hline(yintercept = background_lum, linetype = "dashed", colour = "tomato") +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0,10) 
    ) +
    scale_y_lum(limits = c(1, 300000)) +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("LUM/OD"[600]),
      color = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}


# =============================================================================
# PLOTTING — AVERAGED KINETICS
# =============================================================================

#' Plot averaged OD600 kinetics
#'
#' @param .data Averaged dataframe (df_av)
#' @param strain_to_plot Strain ID to filter to (string)
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label for colour legend (string)
#' @param inductiontime Time of induction (hours)
#' @param grid_col Bare column name for facet columns
#' @param grid_row Bare column name for facet rows
plot_average_od <- function(.data,
                             strain_to_plot,
                             variable_to_colour,
                             unit_of_colour,
                             inductiontime,
                             grid_col,
                             grid_row) {
  
  .data |>
    dplyr::filter(strain == strain_to_plot) |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = as.numeric(od600_blank_mean),
      color = as.factor(.data[[variable_to_colour]])
    )) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin   = od600_blank_mean - od600_blank_sd,
        ymax   = od600_blank_mean + od600_blank_sd,
        colour = as.factor(.data[[variable_to_colour]])
      ),
      alpha = 0.4
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars({{ grid_col }}),
      rows = ggplot2::vars({{ grid_row }})
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid", colour = "grey") +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0,10)
    ) +
    ggplot2::scale_y_log10(
      labels       = scales::label_number(),
      minor_breaks = minor_breaks,
      limits       = c(0.01, 2.0)
    ) +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("OD"[600]),
      color = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}

#' Plot averaged luminescence kinetics (all strains, no time delay)
#'
#' @param .data Averaged dataframe (df_av)
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label for colour legend (string)
#' @param inductiontime Time of induction (hours)
#' @param background_lum Background luminescence value for reference line
#' @param grid_col Bare column name for facet columns
#' @param grid_row Bare column name for facet rows
plot_average_lum <- function(.data,
                              variable_to_colour,
                              unit_of_colour,
                              inductiontime,
                              background_lum,
                              grid_col,
                              grid_row) {
  
  .data |>
    ggplot2::ggplot(ggplot2::aes(
      x     = time,
      y     = as.numeric(lum_blank_mean),
      color = as.factor(.data[[variable_to_colour]])
    )) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin   = pmax(lum_blank_mean - lum_blank_sd, 50),
        ymax   = lum_blank_mean + lum_blank_sd,
        colour = as.factor(.data[[variable_to_colour]])
      ),
      alpha = 0.4
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars({{ grid_col }}),
      rows = ggplot2::vars({{ grid_row }})
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid", colour = "grey") +
    ggplot2::geom_hline(yintercept = background_lum, linetype = "dashed", colour = "tomato") +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0,10)
    ) +
    scale_y_lum() +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("LUM/OD"[600]),
      color = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}

#' Plot luminescence kinetics with time delay markers and spline fit
#'
#' Shows measured data points, the spline fit used for interpolation,
#' dashed vertical lines at the interpolated delay time per condition,
#' and delay labels in minutes.
#'
#' @param .data Averaged dataframe (df_av)
#' @param strain_to_plot Strain ID to filter to (string)
#' @param df_timedelay Time delay table from get_timedelay_splinefun()
#' @param inductiontime Time of induction (hours)
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label for colour legend (string)
#' @param background_lum Background luminescence value for reference line
#' @param grid_col Facet column variable (string)
#' @param grid_row Facet row variable (string)
#' @param group_vars Character vector of grouping columns (needed for spline faceting)
plot_lumkinetic_with_timedelay <- function(.data,
                                            strain_to_plot,
                                            df_timedelay,
                                            inductiontime,
                                            variable_to_colour,
                                            unit_of_colour,
                                            background_lum,
                                            grid_col,
                                            grid_row,
                                            group_vars) {
  
  data_strain  <- .data        |> dplyr::filter(strain == strain_to_plot)
  delay_strain <- df_timedelay |> dplyr::filter(strain == strain_to_plot)
  
  # Generate spline fit on fine grid for each condition
  spline_df <- data_strain |>
    dplyr::filter(time > inductiontime) |>
    dplyr::reframe(
      {
        t  <- time
        fc <- lum_blank_mean
        if (length(t) >= 3 && !all(is.na(fc))) {
          spline_fn <- splinefun(t, fc, method = "monoH.FC")
          t_fine    <- seq(min(t), max(t), length.out = 500)
          data.frame(time = t_fine, lum_spline = spline_fn(t_fine))
        } else {
          data.frame(time = numeric(0), lum_spline = numeric(0))
        }
      },
      .by = all_of(c(group_vars, "strain"))
    )
  
  ggplot2::ggplot(
    data_strain,
    ggplot2::aes(
      x      = time,
      y      = lum_blank_mean,
      colour = as.factor(.data[[variable_to_colour]])
    )
  ) +
    ggplot2::geom_point(size = 1, alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(lum_blank_mean - lum_blank_sd, 10),
        ymax = lum_blank_mean + lum_blank_sd
      ),
      alpha = 0.2
    ) +
    # Spline fit line
    ggplot2::geom_line(
      data      = spline_df,
      ggplot2::aes(x = time, y = lum_spline,
                   colour = as.factor(.data[[variable_to_colour]])),
      linewidth = 0.5,
      alpha     = 0.6
    ) +
    # Induction time
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid",  colour = "grey") +
    # Background line
    ggplot2::geom_hline(yintercept = background_lum, linetype = "solid", colour = "steelblue") +
    # Delay lines
    ggplot2::geom_vline(
      data     = delay_strain,
      ggplot2::aes(xintercept = delay,
                   colour     = as.factor(.data[[variable_to_colour]])),
      linetype = "dashed"
    ) +
    # Delay labels in minutes
    ggplot2::geom_text(
      data        = delay_strain,
      ggplot2::aes(
        x      = delay,
        label  = paste0(round((delay - inductiontime) * 60, 0), " min"),
        colour = as.factor(.data[[variable_to_colour]])
      ),
      y           = log10(15),
      angle       = 90,
      vjust       = -0.4,
      hjust       = 0,
      size        = 2.5,
      show.legend = FALSE
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data[[grid_col]]),
      rows = ggplot2::vars(.data[[grid_row]])
    ) +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 8),
      expand = ggplot2::expansion(),
      limits = c(0,10)
    ) +
    scale_y_lum() +
    ggplot2::labs(
      x      = "time (h)",
      y      = expression("RLU/OD"[600]),
      colour = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}

#' Plot luminescence kinetics with 0 vs max aTc overlaid at max arabinose
#'
#' Designed to be combined with plot_lumkinetic_with_timedelay() via patchwork
#' to address reviewer comment on direct aTc comparison.
#'
#' @param .data Averaged dataframe (df_av)
#' @param strain_to_plot Strain ID (string)
#' @param df_timedelay Time delay table from get_timedelay_splinefun()
#' @param inductiontime Time of induction (hours)
#' @param background_lum Background luminescence for reference line
#' @param group_vars Character vector of grouping columns
#' @param atc_low Lower aTc concentration (default 0)
#' @param atc_high Higher aTc concentration (default 2.5)
#' @param grid_row Facet row variable (string, default "antisigma")
plot_lumkinetic_atc_overlay <- function(.data,
                                        strain_to_plot,
                                        df_timedelay,
                                        inductiontime,
                                        background_lum,
                                        group_vars,
                                        threshold_inducer = "aTc",  # "aTc" or "DHBA"
                                        inducer_low,                 # e.g. 0
                                        inducer_high,                # e.g. 2.5 or 250
                                        grid_row = "antisigma") {
  
  max_arab <- max(.data$arabinose, na.rm = TRUE)
  
  data_overlay <- .data |>
    dplyr::filter(
      strain    == strain_to_plot,
      arabinose == max_arab,
      .data[[threshold_inducer]] %in% c(inducer_low, inducer_high)
    )
  
  delay_overlay <- df_timedelay |>
    dplyr::filter(
      strain    == strain_to_plot,
      arabinose == max_arab,
      .data[[threshold_inducer]] %in% c(inducer_low, inducer_high)
    )
  
  spline_overlay <- data_overlay |>
    dplyr::filter(time > inductiontime) |>
    dplyr::reframe(
      {
        t  <- time
        fc <- lum_blank_mean
        if (length(t) >= 3 && !all(is.na(fc))) {
          spline_fn <- splinefun(t, fc, method = "monoH.FC")
          t_fine    <- seq(min(t), max(t), length.out = 500)
          data.frame(time = t_fine, lum_spline = spline_fn(t_fine))
        } else {
          data.frame(time = numeric(0), lum_spline = numeric(0))
        }
      },
      .by = all_of(c("strain", "antisigma", "arabinose", threshold_inducer))
    )
  
  # Build inducer label for legend
  inducer_unit  <- dplyr::case_when(
    threshold_inducer == "aTc"  ~ "ng/ml",
    threshold_inducer == "DHBA" ~ "µM"
  )
  inducer_label <- paste0(threshold_inducer, " (", inducer_unit, ")")
  
  ggplot2::ggplot(
    data_overlay,
    ggplot2::aes(
      x        = time,
      y        = lum_blank_mean,
      colour   = as.factor(.data[[threshold_inducer]]),
      linetype = as.factor(.data[[threshold_inducer]])
    )
  ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(lum_blank_mean - lum_blank_sd, 10),
        ymax = lum_blank_mean + lum_blank_sd
      ),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = spline_overlay,
      ggplot2::aes(
        x        = time,
        y        = lum_spline,
        colour   = as.factor(.data[[threshold_inducer]]),
        linetype = as.factor(.data[[threshold_inducer]])
      ),
      linewidth = 0.8, alpha = 0.8
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid", colour = "grey") +
    ggplot2::geom_hline(yintercept = background_lum, linetype = "solid", colour = "steelblue") +
    ggplot2::geom_vline(
      data     = delay_overlay,
      ggplot2::aes(
        xintercept = delay,
        colour     = as.factor(.data[[threshold_inducer]])
      ),
      linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data        = delay_overlay,
      ggplot2::aes(
        x      = delay,
        label  = paste0(round((delay - inductiontime) * 60, 0), " min"),
        colour = as.factor(.data[[threshold_inducer]])
      ),
      y = log10(15), angle = 90, vjust = -0.4, hjust = 0,
      size = 2.5, show.legend = FALSE
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data[[grid_row]])) +
    ggplot2::scale_color_manual(
      values = c("cornflowerblue", "darkorchid3"),
      labels = c(paste0(inducer_low,  " ", inducer_unit),
                 paste0(inducer_high, " ", inducer_unit))
    ) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed"),
      labels = c(paste0(inducer_low,  " ", inducer_unit),
                 paste0(inducer_high, " ", inducer_unit))
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(),
      limits = c(0, 10)
    ) +
    scale_y_lum() +
    ggplot2::labs(
      x        = "time (h)",
      y        = NULL,
      colour   = inducer_label,
      linetype = inducer_label,
      title    = paste0("arabinose: ", max_arab, "%")
    ) +
    theme_timer() +
    logticks_left()
}
#' Plot luminescence with/without antisigma overlay (max arabinose only)
#'
#' @param .data Averaged dataframe (df_av)
#' @param strain_to_plot Strain ID (string)
#' @param df_timedelay Time delay table
#' @param inductiontime Time of induction (hours)
#' @param background_lum Background luminescence value
plot_lum_with_withoutAS_overlay <- function(.data,
                                             strain_to_plot,
                                             df_timedelay,
                                             inductiontime,
                                             background_lum, 
                                             variable_to_colour, 
                                             unit_of_colour) {
  
  delay_strain <- df_timedelay |> dplyr::filter(strain == strain_to_plot)
  
  .data |>
    dplyr::filter(strain == strain_to_plot) |>
    ggplot2::ggplot(ggplot2::aes(
      x      = time,
      y      = lum_blank_mean,
      colour = as.factor(.data[[variable_to_colour]])
    )) +
    ggplot2::geom_point(size = 1, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid",  colour = "grey") +
    ggplot2::geom_vline(
      data     = delay_strain,
      ggplot2::aes(xintercept = delay),
      linetype = "dashed",
      colour   = "grey"
    ) +
    ggplot2::geom_hline(yintercept = background_lum, linetype = "solid", colour = "steelblue") +
    ggplot2::geom_text(
      data        = delay_strain,
      ggplot2::aes(
        x     = delay,
        label = paste0(round((delay - inductiontime) * 60, 0), " min")
      ),
      y           = log10(15),
      angle       = 90,
      vjust       = -0.4,
      hjust       = 0,
      size        = 2.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin   = pmax(lum_blank_mean - lum_blank_sd, 10),
        ymax   = lum_blank_mean + lum_blank_sd,
        colour = as.factor(aTc)
      ),
      alpha = 0.2
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(arabinose), rows = ggplot2::vars(antisigma)) +
    ggplot2::scale_color_manual(values = c("cornflowerblue", "darkorchid3") ) + # For 384-well c("lightskyblue1", "cornflowerblue","plum1", "darkorchid3") 
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 8),
      expand = ggplot2::expansion(),
      limits = c(0,10)
    ) +
    scale_y_lum() +
    ggplot2::labs(
      x      = "time (h)",
      y      = expression("RLU/OD"[600]),
      colour = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}


# =============================================================================
# PLOTTING — HEATMAPS
# =============================================================================

#' Plot time delay as a heatmap over inducer concentration grid
#'
#' @param .data Time delay dataframe from get_timedelay_splinefun()
#' @param fill_variable Bare column name for fill values (e.g. delay)
#' @param x_data Column name for x-axis (string)
#' @param y_data Column name for y-axis (string)
#' @param x_unit Unit label for x-axis (string)
#' @param y_unit Unit label for y-axis (string)
#' @param inductiontime Time of induction (hours) — used to convert delay to minutes
plot_timedelay_heatmap <- function(.data,
                                    fill_variable,
                                    x_data,
                                    y_data,
                                    x_unit,
                                    y_unit,
                                    inductiontime) {
  
  .data |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(.data[[x_data]]),
      y    = as.factor(.data[[y_data]]),
      fill = ({{ fill_variable }} - inductiontime) * 60
    )) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = round(({{ fill_variable }} - inductiontime) * 60, 1)),
      vjust = 1.25,
      color = "black",
      size  = 3
    ) +
    ggplot2::scale_fill_gradient(low = "cornflowerblue", high = "brown2", na.value = "white") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      x    = paste0(x_data, " (", x_unit, ")"),
      y    = paste0(y_data, " (", y_unit, ")"),
      fill = "delay (min)"
    ) +
    theme_timer()
}

#' Plot dynamic range (fold-change) as a heatmap over inducer concentration grid
#'
#' @param .data Dynamic range dataframe
#' @param fill_variable Bare column name for fill values
#' @param x_data Column name for x-axis (string)
#' @param y_data Column name for y-axis (string)
#' @param x_unit Unit label for x-axis (string)
#' @param y_unit Unit label for y-axis (string)
plot_dynamic_range_heatmap <- function(.data,
                                        fill_variable,
                                        x_data,
                                        y_data,
                                        x_unit,
                                        y_unit) {
  
  .data |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(.data[[x_data]]),
      y    = as.factor(.data[[y_data]]),
      fill = {{ fill_variable }}
    )) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = round({{ fill_variable }}, 1)),
      vjust = 1.25,
      color = "black",
      size  = 3
    ) +
    ggplot2::scale_fill_gradient(low = "cornflowerblue", high = "brown2", na.value = "white") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(
      x    = paste0(x_data, " (", x_unit, ")"),
      y    = paste0(y_data, " (", y_unit, ")"),
      fill = "fold change"
    ) +
    theme_timer()
}


# =============================================================================
# PLOTTING — DYNAMIC RANGE BAR PLOTS
# =============================================================================

#' Plot maximum fold-change per strain as paired bar plot (threshold ON vs OFF)
#'
#' @param .data Fold-change summary dataframe from average_dynamic_range()
#' @param strain_group Character vector of strain IDs to include
#' @param threshold_inducer String: "aTc" or "DHBA"
plot_dynamic_range_averages <- function(.data, strain_group, threshold_inducer) {
  
  inducer_label <- dplyr::case_when(
    threshold_inducer == "aTc"  ~ "aTc (ng/ml)",
    threshold_inducer == "DHBA" ~ "DHBA (µM)"
  )
  
  .data |>
    dplyr::filter(strain %in% strain_group) |>
    dplyr::filter(
      fold_change_mean == max(fold_change_mean),
      .by = c(antisigma, strain, all_of(threshold_inducer))
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x    = antisigma , #reorder(antisigma, dplyr::desc(fold_change_mean)),
      y    = fold_change_mean,
      fill = as.factor(.data[[threshold_inducer]])
    )) +
    ggplot2::geom_bar(
      stat     = "identity",
      width    = 0.6,
      position = ggplot2::position_dodge(width = 0.7)
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(fold_change_mean - fold_change_sd, 1), #
        ymax = fold_change_mean + fold_change_sd
      ),
      width    = 0.2,
      position = ggplot2::position_dodge(width = 0.7),
      colour   = "black"
    ) +
    ggplot2::scale_fill_manual(values = c("cornflowerblue", "darkorchid3")) + # For 384-well c("lightskyblue1", "cornflowerblue","plum1", "darkorchid3")
    scale_y_fc() +
    ggplot2::labs(
      x    = expression("anti-" * sigma),
      y    = "Max. fold induction",
      fill = inducer_label
    ) +
    theme_timer() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    ) +
    logticks_left()
}



plot_dynamic_range_with_annotation <- function(.data,
                                               strain_group,
                                               threshold_inducer,
                                               sig_results) {
  # --- Prepare plot data (same as plot_dynamic_range_averages) ---
  plot_df <- .data |>
    dplyr::filter(strain %in% strain_group) |>
    dplyr::filter(
      fold_change_mean == max(fold_change_mean),
      .by = c(antisigma, strain, all_of(threshold_inducer))
    )
  
  # Determine x-axis order from the bar plot (descending fold change)
  antisigma_order <- plot_df |>
    dplyr::filter(.data[[threshold_inducer]] == 0) |>
    dplyr::arrange(dplyr::desc(fold_change_mean)) |>
    dplyr::pull(antisigma) |>
    unique()
  
  plot_df <- plot_df |>
    dplyr::mutate(antisigma = factor(antisigma, levels = antisigma_order))
  
  # --- Main bar plot ---
  p_bars <- plot_df |>
    ggplot2::ggplot(ggplot2::aes(
      x    = antisigma,
      y    = fold_change_mean,
      fill = as.factor(.data[[threshold_inducer]])
    )) +
    ggplot2::geom_bar(
      stat     = "identity",
      width    = 0.6,
      position = ggplot2::position_dodge(width = 0.7)
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(fold_change_mean - fold_change_sd, 0.05),
        ymax = fold_change_mean + fold_change_sd
      ),
      width    = 0.2,
      position = ggplot2::position_dodge(width = 0.7),
      colour   = "black"
    ) +
    ggplot2::scale_fill_manual(values = c("cornflowerblue", "darkorchid3")) +
    scale_y_fc() +
    ggplot2::labs(
      x    = NULL,          # remove x label — annotation panel provides context
      y    = "Max. fold induction",
      fill = dplyr::case_when(
        threshold_inducer == "aTc"  ~ "aTc (ng/ml)",
        threshold_inducer == "DHBA" ~ "DHBA (µM)"
      )
    ) +
    theme_timer() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),  # hide — annotation panel shows labels
      axis.ticks.x = ggplot2::element_blank()
    ) +
    logticks_left()
  
  # --- Annotation data ---
  # Join significance results and define inducibility
  # "Inducible" = fold change > 4 without AS (baseline expressible)
  # "Significant" = p < 0.05 from test_threshold_significance()
  annot_df <- plot_df |>
    dplyr::filter(.data[[threshold_inducer]] == 0) |>
    dplyr::select(antisigma, strain, fold_change_mean) |>
    dplyr::left_join(
      sig_results |> dplyr::select(strain, p_value, significance),
      by = "strain"
    ) |>
    dplyr::mutate(
      antisigma   = factor(antisigma, levels = antisigma_order),
      inducible   = fold_change_mean > 4,
      significant = !is.na(p_value) & p_value < 0.05
    ) |>
    tidyr::pivot_longer(
      cols      = c(inducible, significant),
      names_to  = "property",
      values_to = "value"
    ) |>
    dplyr::mutate(
      property = factor(property,
                        levels   = c("inducible", "significant"),
                        labels   = c("ECF inducible", paste0("AS effect ", significance_label(threshold_inducer))))
    )
  
  # --- Annotation tile panel ---
  p_annot <- annot_df |>
    ggplot2::ggplot(ggplot2::aes(
      x    = antisigma,
      y    = property,
      fill = value
    )) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#4daf4a", "FALSE" = "#e0e0e0"),
      labels = c("TRUE" = "yes", "FALSE" = "no"),
      name   = NULL
    ) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = expression("anti-" * sigma), y = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.ticks       = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      legend.position  = "none",
      panel.border     = ggplot2::element_rect(colour = "grey80")
    )
  
  # --- Combine with patchwork ---
  # Heights: bar plot gets ~75% of space, annotation ~25%
  p_bars / p_annot +
    patchwork::plot_layout(heights = c(3, 1))
}

# Helper for legend label
significance_label <- function(threshold_inducer) {
  dplyr::case_when(
    threshold_inducer == "aTc"  ~ "(aTc)",
    threshold_inducer == "DHBA" ~ "(DHBA)"
  )
}






#' Plot delay vs fold-change scatter plot
#'
#' @param .data Joined delays + fold-change dataframe (df_delays_fc)
#' @param inductiontime Time of induction (hours)
#' @param as_inducer Bare column name of antisigma inducer to colour by
#' @param legend Logical. TRUE shows full legend; FALSE plots all points black
plot_delay_vs_fc <- function(.data, inductiontime, as_inducer, legend = TRUE) {
  
  base <- .data |>
    dplyr::filter(fc_mean > 4) |>
    ggplot2::ggplot(ggplot2::aes(
      x = (delay - inductiontime) * 60,
      y = fc_mean
    )) +
    ggplot2::scale_y_log10(
      minor_breaks = minor_breaks,
      limits       = c(1, 1000),
      expand       = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 500), expand = c(0, 0)) +
    ggplot2::labs(x = "time delay (min)", y = "fold-induction") +
    theme_timer() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 8))
  
  if (legend) {
    base +
      ggplot2::geom_point(ggplot2::aes(
        color = as.factor({{ as_inducer }}),
        shape = as.factor(antisigma),
        size  = arabinose
      )) +
      ggplot2::scale_color_viridis_d(direction = -1) +
      ggplot2::labs(
        color = "AS inducer",
        shape = expression("anti-" * sigma),
        size  = "arabinose (%)"
      )
  } else {
    base +
      ggplot2::geom_point(size = 1, colour = "black")
  }
}

#' Plot dose-response curves (luminescence at final time point vs arabinose)
#'
#' @param .data Averaged dataframe (df_av), pre-filtered to strain of interest
plot_doseresponse <- function(.data) {
  
  .data |>
    dplyr::filter(time == max(time)) |>
    ggplot2::ggplot(ggplot2::aes(x = arabinose, y = lum_blank_mean)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(lum_blank_mean - lum_blank_sd, 1),
        ymax = lum_blank_mean + lum_blank_sd
      ),
      width = 0.05
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(aTc)) +
    ggplot2::scale_y_log10(
      labels       = function(x) parse(text = paste0("10^", log10(x))),
      minor_breaks = minor_breaks,
      limits       = c(1, 300000)
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_log(base = 10),
      limits = c(0.000001, 0.0002),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "arabinose (%)",
      y = expression("RLU/OD"[600])
    ) +
    theme_timer() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )
}
