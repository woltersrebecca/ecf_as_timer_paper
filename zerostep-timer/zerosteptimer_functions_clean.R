# =============================================================================
# zerosteptimer_functions_clean.R
# Functions for zero-step timer luminescence plate reader analysis
# Author: Rebecca Wolters
#
# The zero-step timer characterises inducible promoters (Ptet, Pbad, PpcaU)
# driving luciferase directly — no antisigma cascade. Key metadata columns:
#   promoter    — inducible promoter name (e.g. "Ptet", "Pbad", "PpcaU")
#   copy        — plasmid copy number ("sc" = single copy, "mc" = multi copy)
#   concentration — inducer concentration (aTc in ng/ml, arabinose in %, DHBA in µM)
#   strain      — strain ID (e.g. "GFC0200")

# =============================================================================

# --- Shared plot theme -------------------------------------------------------

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

scale_y_lum <- function(limits = c(10, 200000), minor_breaks = NULL) {
  ggplot2::scale_y_log10(
    labels       = function(x) parse(text = paste0("10^", log10(x))),
    minor_breaks = minor_breaks,
    limits       = limits
  )
}

scale_y_fc <- function(limits = c(0.5, 10000)) {
  ggplot2::scale_y_log10(
    labels = function(x) parse(text = paste0("10^", log10(x))),
    limits = limits
  )
}

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
    dplyr::mutate(time = time / 60 / 60)
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

#' Plot blank readings before induction for QC
#'
#' @param .data Raw joined dataframe
#' @param blank_variable Column name to plot (string)
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
    theme_timer()
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

#' Background correct OD600 when NO blank well on plate
#'
#' @param data Raw joined dataframe
#' @param average_blank_od Average blank OD (from get_blank_average)
#' @param pc_factor Pathlength correction factor (96-well: 4.898, 384-well: 3.09)
blank_plate <- function(data, average_blank_od, pc_factor) {
  data |>
    tidyr::drop_na(strain) |>
    dplyr::mutate(
      od600_blank = od600 * pc_factor - average_blank_od,
      lum_blank   = lum / od600_blank
    )
}

#' Background correct OD600 when blank well IS on the plate
#'
#' @param .data Raw joined dataframe
#' @param pc_factor Pathlength correction factor
#' @param group_vars Character vector of columns to group blank subtraction by
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
#' @param .data Averaged dataframe (df_av)
#' @param value Bare column name of luminescence variable (e.g. lum_blank_mean)
#' @param starttime Start of baseline window (hours)
#' @param inductiontime End of baseline window / induction time (hours)
#' @param group_vars Character vector of grouping columns
#' @param wildtype Logical. TRUE = use blank strain; FALSE = per-condition baseline
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
    df_pre <- df_pre |> dplyr::filter(strain == "blank")
  }

  df_pre |>
    dplyr::summarise(
      background    = mean({{ value }}, na.rm = TRUE),
      background_sd = sd({{ value }},   na.rm = TRUE),
      .by           = all_of(group_vars)
    )
}


# =============================================================================
# TIME DELAY — INTERPOLATED (recommended)
# =============================================================================

#' Calculate time delay using cubic spline interpolation
#'
#' Fits a monotone cubic spline through post-induction luminescence values
#' and uses Brent's root-finding (uniroot) to find the exact time at which
#' luminescence crosses foldchange × pre-induction baseline.
#' This removes quantisation error from discrete measurement intervals.
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
#' @param .data Blanked dataframe (df_blank)
#' @param bg_expr Background expression table from get_background_lum()
#' @param group_vars Character vector of grouping columns
#' @param endpoint_hours Time point to use as endpoint (hours). NULL = max(time).
average_dynamic_range <- function(.data, bg_expr, group_vars, endpoint_hours = NULL) {
  df <- .data |>
    dplyr::left_join(bg_expr, by = group_vars) |>
    dplyr::mutate(
      fold_change = pmax(lum_blank / background, 1),
      .by         = all_of(c(group_vars, "exp"))
    )

  if (!is.null(endpoint_hours)) {
    df <- df |>
      dplyr::filter(
        abs(time - endpoint_hours) == min(abs(time - endpoint_hours)),
        .by = all_of(c(group_vars,"exp"))
      )
  } else {
    df <- df |>
      dplyr::filter(
        time == max(time),
        .by = all_of(c(group_vars,"exp"))
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

#' Calculate fold-change at all time points (for delay vs FC scatter plots)
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
    dplyr::mutate(gr = estimate, dt = log(2) / gr) |>
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
    ggplot2::scale_color_viridis_d() +
    ggplot2::scale_y_log10(minor_breaks = minor_breaks) +
    ggplot2::labs(color = variable_to_color) +
    theme_timer()
}


# =============================================================================
# PLOTTING — AVERAGED KINETICS
# =============================================================================

#' Plot averaged OD600 kinetics
#'
#' Zero-step specific: pre-filtered to a strain subset before calling —
#' the strain_to_plot argument is used only for the plot title.
#'
#' @param .data Averaged dataframe (df_av), pre-filtered to strains of interest
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label for colour legend (string)
#' @param strain_to_plot Label for plot title (string)
#' @param inductiontime Time of induction (hours)
#' @param grid_col Column name for facet columns (string)
#' @param grid_row Column name for facet rows (string)
plot_average_od <- function(.data,
                             variable_to_colour,
                             unit_of_colour,
                             strain_to_plot = NULL,
                             inductiontime,
                             grid_col,
                             grid_row) {
  .data |>
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
      cols = ggplot2::vars(.data[[grid_col]]),
      rows = ggplot2::vars(.data[[grid_row]])
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid", colour = "grey") +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0, 10)
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(),
      limits = c(0.05, 1)
    ) +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("OD"[600]),
      title = strain_to_plot,
      color = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}

#' Plot averaged luminescence kinetics
#'
#' @param .data Averaged dataframe (df_av), pre-filtered to strains of interest
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label for colour legend (string)
#' @param inductiontime Time of induction (hours)
#' @param background_lum Background luminescence for reference line
#' @param grid_col Column name for facet columns (string)
#' @param grid_row Column name for facet rows (string)
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
        ymin   = pmax(lum_blank_mean - lum_blank_sd, 1),
        ymax   = lum_blank_mean + lum_blank_sd,
        colour = as.factor(.data[[variable_to_colour]])
      ),
      alpha = 0.4
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data[[grid_col]]),
      rows = ggplot2::vars(.data[[grid_row]])
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid",  colour = "grey") +
    ggplot2::geom_hline(yintercept = background_lum, linetype = "dashed", colour = "tomato") +
    ggplot2::scale_color_viridis_d(direction = -1) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(),
      expand = ggplot2::expansion(),
      limits = c(0, 10)
    ) +
    scale_y_lum() +
    ggplot2::labs(
      x     = "time (h)",
      y     = expression("LUM/OD"[600]),
      color = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    theme_timer()
}

#' Plot luminescence kinetics with interpolated time delay and spline fit
#'
#' Zero-step specific: filters by strain %in% strain_to_plot (vector accepted),
#' colours by inducer concentration, facets by copy number × promoter.
#'
#' @param .data Averaged dataframe (df_av)
#' @param strain_to_plot Strain ID or vector of IDs (string or character vector)
#' @param df_timedelay Time delay table from get_timedelay_splinefun()
#' @param inductiontime Time of induction (hours)
#' @param variable_to_colour Column to colour by (string)
#' @param unit_of_colour Unit label (string)
#' @param background_lum Background luminescence for reference line
#' @param grid_col Facet column variable (string)
#' @param grid_row Facet row variable (string)
#' @param group_vars Character vector of grouping columns
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

  data_strain  <- .data        |> dplyr::filter(strain %in% strain_to_plot)
  delay_strain <- df_timedelay |> dplyr::filter(strain %in% strain_to_plot)

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
    # Spline fit
    ggplot2::geom_line(
      data      = spline_df,
      ggplot2::aes(x = time, y = lum_spline,
                   colour = as.factor(.data[[variable_to_colour]])),
      linewidth = 0.5, alpha = 0.6
    ) +
    ggplot2::geom_vline(xintercept = inductiontime, linetype = "solid",  colour = "grey") +
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
      limits = c(0, 10)
    ) +
    scale_y_lum() +
    ggplot2::labs(
      x      = "time (h)",
      y      = expression("RLU/OD"[600]),
      colour = paste0(variable_to_colour, " (", unit_of_colour, ")")
    ) +
    ggplot2::theme(
      legend.position   = "right",
      panel.grid.major  = ggplot2::element_line(color = "transparent"),
      panel.grid.minor  = ggplot2::element_line(color = "transparent"),
      axis.ticks.length = ggplot2::unit(-0.1, "cm"),
      panel.spacing.x   = ggplot2::unit(0.5, "cm")
    )
}


# =============================================================================
# PLOTTING — DYNAMIC RANGE
# =============================================================================

#' Plot maximum fold-change as bar plot
#'
#' Zero-step specific: x-axis is inducer concentration, bars filled by
#' promoter, faceted by copy number. Shows fold-change at the highest
#' concentration for each promoter × copy combination.
#'
#' @param .data Fold-change summary dataframe from average_dynamic_range()
#' @param promoter_filter Character vector of promoter names to include.
#'   NULL includes all promoters.
plot_dynamic_range_averages <- function(.data, promoter_filter = NULL) {

  if (!is.null(promoter_filter)) {
    .data <- .data |> dplyr::filter(promoter %in% promoter_filter)
  }

  .data |>
    dplyr::filter(
      fold_change_mean == max(fold_change_mean),
      .by = c(promoter, concentration, copy, strain)
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x    = as.factor(concentration),
      y    = fold_change_mean,
      fill = as.factor(concentration)
    )) +
    ggplot2::geom_col(
      stat     = "identity",
      width    = 0.5,
      position = ggplot2::position_dodge(width = 0.5),
      colour   = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(fold_change_mean - fold_change_sd, 0.5),
        ymax = fold_change_mean + fold_change_sd
      ),
      width    = 0.2,
      position = ggplot2::position_dodge(width = 0.5),
      colour   = "black"
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(copy)) +
    ggplot2::scale_fill_viridis_d(direction = -1) +
    scale_y_fc() +
    ggplot2::labs(
      x    = "inducer concentration",
      y    = "Max. fold induction",
      fill = "concentration"
    ) +
    ggplot2::theme(
      legend.position   = "right",
      panel.grid.major  = ggplot2::element_line(color = "transparent"),
      panel.grid.minor  = ggplot2::element_line(color = "transparent"),
      axis.ticks.length = ggplot2::unit(-0.1, "cm"),
      panel.spacing.x   = ggplot2::unit(0.5, "cm")
    )+
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0.5)
    ) +
    logticks_left()
}

#' Plot delay vs fold-change scatter plot
#'
#' Zero-step specific: coloured by copy number, shaped by promoter.
#' inductiontime is passed explicitly to convert delay to minutes.
#'
#' @param .data Joined delays + fold-change dataframe (df_delays_fc)
#' @param inductiontime Time of induction (hours)
#' @param legend Logical. TRUE shows full legend; FALSE = minimal styling
plot_delay_vs_fc <- function(.data, inductiontime, legend = TRUE) {

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
    ggplot2::theme(
      axis.ticks.length = ggplot2::unit(-0.1, "cm"),
      axis.text         = ggplot2::element_text(size = 8),
      axis.title        = ggplot2::element_text(size = 8)
    )

  if (legend) {
    base +
      ggplot2::geom_point(
        ggplot2::aes(
          color = as.factor(copy),
          shape = as.factor(promoter),
          size  = concentration
        ),
        alpha = 0.7
      ) +
      ggplot2::scale_color_manual(values = c("cornflowerblue", "tomato")) +
      ggplot2::labs(
        color = "copy number",
        shape = "promoter",
        size  = "concentration"
      )
  } else {
    base +
      ggplot2::geom_point(
        ggplot2::aes(color = as.factor(copy)),
        size  = 1,
        alpha = 0.7
      ) +
      ggplot2::scale_color_manual(values = c("black", "grey40")) +
      ggplot2::theme(legend.position = "none")
  }
}

#' Plot dose-response curves (luminescence at final time point vs concentration)
#'
#' Zero-step specific: filters by promoter, facets by copy number.
#' Uses log10 x-axis for concentration and summarises replicates with
#' geom_point + geom_line via stat_summary.
#'
#' @param .data Averaged dataframe (df_av)
#' @param promoter_name Promoter to filter to (string, e.g. "Ptet")
#' @param inducer Label for x-axis (string, e.g. "aTc")
#' @param unit_of_inducer Unit label for x-axis (string, e.g. "ng/ml")
plot_doseresponse <- function(.data, promoter_name, inducer, unit_of_inducer) {
  .data |>
    dplyr::filter(
      time    == max(time),
      promoter == promoter_name
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = concentration, y = lum_blank_mean)) +
    ggplot2::stat_summary(fun = mean, geom = "point") +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(lum_blank_mean - lum_blank_sd, 1),
        ymax = lum_blank_mean + lum_blank_sd
      ),
      width = 0.05
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(copy)) +
    ggplot2::scale_y_log10(
      labels       = function(x) parse(text = paste0("10^", log10(x))),
      minor_breaks = minor_breaks,
      limits       = c(1, 300000)
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_comma(),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = paste0(inducer, " (", unit_of_inducer, ")"),
      y = expression("RLU/OD"[600])
    ) +
    theme_timer() +
    ggplot2::theme(
      panel.spacing.x = ggplot2::unit(1, "cm"),
      axis.text.x     = ggplot2::element_text(angle = 45, vjust = 0.5, hjust = 0.5)
    )
}
