## ============================================================
## R/simulation_suite.R
## Simulation wrapper, validation metrics, and ggplot2 diagnostics
## ============================================================

# ================================================================
# 1.  Internal metric helpers
# ================================================================

#' Root Mean Square Error over all elements of two conformable matrices
#' @noRd
.rmse <- function(X_hat, X_true)
  sqrt(mean((as.numeric(X_hat) - as.numeric(X_true))^2))

#' Lin's Concordance Correlation Coefficient  (Lin 1989, Biometrics 45:255)
#'
#' Uses population (MLE) variances, vectorising over all elements.
#' rho_c = 2*cov(x,y) / (var(x) + var(y) + (mu_x - mu_y)^2)
#' @noRd
.ccc <- function(x_true, x_hat) {
  x    <- as.numeric(x_true)
  y    <- as.numeric(x_hat)
  mu_x <- mean(x);   mu_y <- mean(y)
  s2_x <- mean((x - mu_x)^2)
  s2_y <- mean((y - mu_y)^2)
  s_xy <- mean((x - mu_x) * (y - mu_y))
  denom <- s2_x + s2_y + (mu_x - mu_y)^2
  if (denom < .Machine$double.eps) return(NA_real_)
  2 * s_xy / denom
}

#' Per-cell metrics: RMSE and CCC computed row-wise
#' @noRd
.cell_metrics <- function(X_hat, X_true) {
  n <- nrow(X_true)
  data.frame(
    rmse_cell = vapply(seq_len(n),
                       function(i) .rmse(X_hat[i, ], X_true[i, ]),
                       numeric(1)),
    ccc_cell  = vapply(seq_len(n),
                       function(i) .ccc(X_true[i, ], X_hat[i, ]),
                       numeric(1))
  )
}

#' Per-year metrics: RMSE and CCC computed column-wise
#' @noRd
.year_metrics <- function(X_hat, X_true) {
  T_yr <- ncol(X_true)
  data.frame(
    year      = seq_len(T_yr),
    rmse_year = vapply(seq_len(T_yr),
                       function(t) .rmse(X_hat[, t], X_true[, t]),
                       numeric(1)),
    ccc_year  = vapply(seq_len(T_yr),
                       function(t) .ccc(X_true[, t], X_hat[, t]),
                       numeric(1))
  )
}


# ================================================================
# 2.  run_simulation_suite()
# ================================================================

#' Run a Hyperparameter Grid Search for Spatio-Temporal Deconvolution
#'
#' @description
#' Sweeps a grid of \eqn{(\lambda_t, \lambda_s)} values, fits the ADMM
#' solver for each combination, and returns RMSE and Lin's Concordance
#' Correlation Coefficient (CCC) against the known ground-truth latent
#' signal.
#'
#' @details
#' ## Efficiency
#'
#' The ADMM system matrix
#' \deqn{A = \tilde{M}^\top W \tilde{M}
#'           + \rho (D_t^\top D_t + D_s^\top D_s) + \delta I}
#' does **not** depend on \eqn{\lambda_t} or \eqn{\lambda_s}.  Therefore
#' `run_simulation_suite()` calls [ADMMSolver]`$setup()` **once** (one
#' sparse Cholesky factorisation) and then re-runs the ADMM loop for each
#' grid point, reusing the factorisation.  This makes grid searches
#' substantially faster than naive repeated calls to [solve_deconv()].
#'
#' ## Metrics
#'
#' * **RMSE** -- root mean squared error over all \eqn{n \times T} elements.
#' * **CCC** -- Lin's (1989) concordance correlation coefficient, which
#'   simultaneously penalises bias (location shift), scale differences,
#'   and correlation:
#'   \deqn{\rho_c = \frac{2\,\text{Cov}(X, \hat{X})}
#'         {\text{Var}(X) + \text{Var}(\hat{X}) +
#'          (\bar{X} - \bar{\hat{X}})^2}.}
#'   Ranges in \eqn{[-1, 1]}; perfect recovery gives \eqn{\rho_c = 1}.
#'
#' Per-cell and per-year breakdowns of each metric are included in the
#' returned `$detail` list.
#'
#' @param gen A [SyntheticDataGenerator] object that has already been
#'   simulated (i.e., `$simulate()` called).
#' @param lambda_t_grid Numeric vector of \eqn{\lambda_t} values to try.
#' @param lambda_s_grid Numeric vector of \eqn{\lambda_s} values to try.
#' @param D   (T-k)  x  T temporal difference matrix from
#'   [build_temporal_penalty()].
#' @param Ls  n  x  n spatial Laplacian from [build_spatial_laplacian()].
#' @param W   (nS)  x  (nS) diagonal precision matrix from
#'   [build_precision_weights()].  Pass `NULL` for unweighted LS.
#' @param rho Positive ADMM penalty parameter (default `1`).
#' @param ridge Small ridge for numerical stability (default `1e-8`).
#' @param max_iter Maximum ADMM iterations per grid point (default `1000`).
#' @param tol_abs,tol_rel Convergence tolerances (defaults `1e-4`, `1e-3`).
#' @param keep_best Logical.  If `TRUE` (default), store `$X_hat` for the
#'   grid point achieving the lowest RMSE.
#' @param verbose Logical.  Print per-combination progress (default `FALSE`).
#'
#' @return An S3 object of class `"SimulationSuite"` -- a named list with:
#' \describe{
#'   \item{`$results`}{`data.frame` with columns `lambda_t`, `lambda_s`,
#'     `rmse`, `ccc`, `converged`, `n_iter`.}
#'   \item{`$best`}{Named list: `lambda_t`, `lambda_s`, `rmse`, `ccc`,
#'     `X_hat` for the best (lowest-RMSE) grid point.  `NULL` if
#'     `keep_best = FALSE`.}
#'   \item{`$detail`}{Named list with `$cell` (per-cell metrics for the
#'     best run) and `$year` (per-year metrics for the best run).}
#'   \item{`$timing_sec`}{Wall-clock seconds for setup and total sweep.}
#'   \item{`$params`}{Copy of key input parameters for reproducibility.}
#' }
#'
#' @examples
#' \dontrun{
#' gen <- SyntheticDataGenerator$new(n_row=6, n_col=6, seed=1L)
#' gen$simulate()
#' D  <- build_temporal_penalty(gen$params$T_years, k = 2)
#' Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
#' W  <- build_precision_weights(matrix(0.4, 36, 16))
#'
#' suite <- run_simulation_suite(
#'   gen           = gen,
#'   lambda_t_grid = c(0.01, 0.05, 0.1, 0.5),
#'   lambda_s_grid = c(0.01, 0.05, 0.1),
#'   D = D, Ls = Ls, W = W,
#'   verbose = TRUE
#' )
#' print(suite)
#' plot_lambda_heatmap(suite)
#' plot_cell_signals(gen, suite$best$X_hat)
#' }
#'
#' @export
run_simulation_suite <- function(gen,
                                 lambda_t_grid,
                                 lambda_s_grid,
                                 D,
                                 Ls,
                                 W        = NULL,
                                 rho      = 1.0,
                                 ridge    = 1e-8,
                                 max_iter = 1000L,
                                 tol_abs  = 1e-4,
                                 tol_rel  = 1e-3,
                                 keep_best = TRUE,
                                 verbose  = FALSE) {

  if (!inherits(gen, "SyntheticDataGenerator"))
    stop("`gen` must be a SyntheticDataGenerator object.")
  if (is.null(gen$X) || is.null(gen$Y))
    stop("`gen` must be fully simulated first (call `gen$simulate()`).")
  if (any(lambda_t_grid < 0) || any(lambda_s_grid < 0))
    stop("All lambda values must be non-negative.")

  lambda_t_grid <- sort(unique(lambda_t_grid))
  lambda_s_grid <- sort(unique(lambda_s_grid))
  param_grid    <- expand.grid(lambda_t = lambda_t_grid,
                               lambda_s = lambda_s_grid,
                               KEEP.OUT.ATTRS = FALSE,
                               stringsAsFactors = FALSE)
  n_combos <- nrow(param_grid)

  # -- One-time setup (single Cholesky factorisation) ------------
  if (verbose)
    cat(sprintf("Setting up ADMM system (%d cells \u00D7 %d years)...\n",
                gen$params$n_cells, gen$params$T_years))
  t_setup_start <- proc.time()[["elapsed"]]

  solver <- ADMMSolver$new(
    Y = gen$Y, M = gen$M, W = W,
    D = D, Ls = Ls,
    lambda_t = param_grid$lambda_t[1L],
    lambda_s = param_grid$lambda_s[1L],
    rho = rho, ridge = ridge
  )
  solver$setup()

  t_setup <- proc.time()[["elapsed"]] - t_setup_start
  if (verbose)
    cat(sprintf("  Setup done in %.2f s. Sweeping %d \u03bb combinations...\n",
                t_setup, n_combos))

  # -- Pre-allocate result storage -------------------------------
  res <- data.frame(
    lambda_t  = param_grid$lambda_t,
    lambda_s  = param_grid$lambda_s,
    rmse      = NA_real_,
    ccc       = NA_real_,
    converged = NA,
    n_iter    = NA_integer_
  )

  best_rmse  <- Inf
  best_X_hat <- NULL
  best_idx   <- NA_integer_

  t_sweep_start <- proc.time()[["elapsed"]]

  for (i in seq_len(n_combos)) {

    # Update lambda without recomputing Cholesky
    solver$lambda_t <- param_grid$lambda_t[i]
    solver$lambda_s <- param_grid$lambda_s[i]

    suppressWarnings(
      solver$run(max_iter   = as.integer(max_iter),
                 tol_abs    = tol_abs,
                 tol_rel    = tol_rel,
                 verbose    = FALSE,
                 warm_start = FALSE)
    )

    rmse_i <- .rmse(solver$X_hat, gen$X)
    ccc_i  <- .ccc(gen$X, solver$X_hat)

    res$rmse[i]      <- rmse_i
    res$ccc[i]       <- ccc_i
    res$converged[i] <- solver$converged
    res$n_iter[i]    <- nrow(solver$history)

    if (keep_best && rmse_i < best_rmse) {
      best_rmse  <- rmse_i
      best_X_hat <- solver$X_hat
      best_idx   <- i
    }

    if (verbose)
      cat(sprintf("  [%3d/%d] \u03bb_t=%-8.4g \u03bb_s=%-8.4g"
                  , i, n_combos,
                  param_grid$lambda_t[i], param_grid$lambda_s[i]),
          sprintf("RMSE=%6.4f  CCC=%6.4f  %s  (%d iter)\n",
                  rmse_i, ccc_i,
                  if (solver$converged) "conv" else "----",
                  nrow(solver$history)))
  }

  t_total <- proc.time()[["elapsed"]] - t_setup_start

  # -- Per-cell / per-year detail for best solution --------------
  detail <- if (keep_best && !is.na(best_idx)) {
    list(
      cell = cbind(
        sf::st_drop_geometry(gen$grid)[, "cell_id", drop = FALSE],
        .cell_metrics(best_X_hat, gen$X)
      ),
      year = .year_metrics(best_X_hat, gen$X)
    )
  } else {
    list(cell = NULL, year = NULL)
  }

  best <- if (keep_best && !is.na(best_idx)) {
    list(
      lambda_t = res$lambda_t[best_idx],
      lambda_s = res$lambda_s[best_idx],
      rmse     = best_rmse,
      ccc      = res$ccc[best_idx],
      X_hat    = best_X_hat,
      idx      = best_idx
    )
  } else NULL

  structure(
    list(
      results       = res,
      best          = best,
      detail        = detail,
      timing_sec    = c(setup = t_setup, total = t_total),
      params        = list(rho = rho, ridge = ridge, max_iter = max_iter,
                           tol_abs = tol_abs, tol_rel = tol_rel,
                           n_cells  = gen$params$n_cells,
                           T_years  = gen$params$T_years,
                           n_combos = n_combos)
    ),
    class = "SimulationSuite"
  )
}


# ================================================================
# 3.  S3 methods for SimulationSuite
# ================================================================

#' Print method for SimulationSuite
#' @param x A `"SimulationSuite"` object.
#' @param ... Ignored.
#' @export
print.SimulationSuite <- function(x, ...) {
  r <- x$results
  cat("\u2500\u2500 SimulationSuite (deconvCore) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  Grid     : %d \u03bb_t \u00d7 %d \u03bb_s = %d combinations\n",
              length(unique(r$lambda_t)), length(unique(r$lambda_s)), nrow(r)))
  cat(sprintf("  Problem  : %d cells \u00d7 %d years\n",
              x$params$n_cells, x$params$T_years))
  cat(sprintf("  Timing   : setup %.2f s  |  total %.2f s\n",
              x$timing_sec["setup"], x$timing_sec["total"]))
  cat(sprintf("  RMSE     : %.4f \u2013 %.4f  (best %.4f at \u03bb_t=%.4g, \u03bb_s=%.4g)\n",
              min(r$rmse, na.rm = TRUE), max(r$rmse, na.rm = TRUE),
              x$best$rmse, x$best$lambda_t, x$best$lambda_s))
  cat(sprintf("  CCC      : %.4f \u2013 %.4f  (best %.4f)\n",
              min(r$ccc,  na.rm = TRUE), max(r$ccc,  na.rm = TRUE),
              x$best$ccc))
  cat(sprintf("  Converged: %d / %d runs\n",
              sum(r$converged, na.rm = TRUE), nrow(r)))
  cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  invisible(x)
}


# ================================================================
# 4.  plot_cell_signals()
# ================================================================

#' Plot True, Observed, and Recovered Signals for One Grid Cell
#'
#' @description
#' Produces a [ggplot2][ggplot2::ggplot2-package] figure overlaying three
#' time series for a single spatial cell:
#'
#' * **True annual signal** \eqn{X_{i,:}} -- the latent ground truth.
#' * **Recovered annual signal** \eqn{\hat{X}_{i,:}} -- the ADMM estimate.
#' * **Observed 5-year moving average** \eqn{Y_{i,:}} -- shown as horizontal
#'   segments spanning each averaging window, emphasising the temporal
#'   smearing that deconvolution reverses.
#'
#' @param gen     A simulated [SyntheticDataGenerator] object.
#' @param X_hat   n  x  T matrix of recovered signals (e.g., `suite$best$X_hat`
#'   or the output of [solve_deconv()]).
#' @param cell_id Integer scalar. Which cell to plot (1-based row index of
#'   `gen$X`). Defaults to a random cell.
#' @param title   Character string for the plot title. Auto-generated if
#'   `NULL`.
#'
#' @return A [ggplot2::ggplot] object (invisible).
#'
#' @examples
#' \dontrun{
#' plot_cell_signals(gen, suite$best$X_hat, cell_id = 42)
#' }
#' @export
plot_cell_signals <- function(gen, X_hat, cell_id = NULL, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_cell_signals().")

  n  <- gen$params$n_cells
  T  <- gen$params$T_years
  w  <- gen$params$window
  S  <- gen$params$n_periods

  if (is.null(cell_id))
    cell_id <- sample.int(n, 1L)
  cell_id <- as.integer(cell_id)
  if (cell_id < 1L || cell_id > n)
    stop("`cell_id` must be between 1 and ", n, ".")

  # -- Extract the three time series for the selected cell -------
  x_true <- as.numeric(gen$X[cell_id, ])
  x_hat  <- as.numeric(X_hat[cell_id, ])
  y_obs  <- as.numeric(gen$Y[cell_id, ])   # length S

  # True and recovered: one point per year
  df_annual <- data.frame(
    year   = seq_len(T),
    True   = x_true,
    Recovered = x_hat
  )

  # Reshape to long format for ggplot2
  df_long <- data.frame(
    year   = rep(seq_len(T), 2L),
    value  = c(x_true, x_hat),
    Series = rep(c("True annual", "Recovered annual"), each = T)
  )

  # Y obs: each period s spans years [s, s + w - 1]; midpoint = s + (w-1)/2
  df_obs <- data.frame(
    xstart = seq_len(S),
    xend   = seq_len(S) + w - 1L,
    y      = y_obs,
    Series = "Observed 5-yr MA"
  )

  # Colour palette: accessible, distinct
  cols <- c("True annual"      = "#2c7bb6",
            "Recovered annual" = "#d7191c",
            "Observed 5-yr MA" = "#fdae61")

  lty  <- c("True annual"      = "solid",
            "Recovered annual" = "dashed",
            "Observed 5-yr MA" = "solid")

  cell_label <- if (!is.null(rownames(gen$X)))
    rownames(gen$X)[cell_id]
  else
    paste0("cell_", cell_id)

  if (is.null(title))
    title <- sprintf("Deconvolution \u2014 %s  (%d\u00d7%d grid, T=%d, w=%d)",
                     cell_label,
                     gen$params$n_row, gen$params$n_col,
                     T, w)

  rmse_cell <- .rmse(x_hat, x_true)
  ccc_cell  <- .ccc(x_true, x_hat)
  subtitle  <- sprintf("RMSE = %.4f    CCC = %.4f", rmse_cell, ccc_cell)

  p <- ggplot2::ggplot() +
    # Observed Y: horizontal segments (one per period)
    ggplot2::geom_segment(
      data = df_obs,
      ggplot2::aes(x    = .data$xstart, xend = .data$xend,
                   y    = .data$y,      yend = .data$y,
                   colour = .data$Series),
      linewidth = 2.0, alpha = 0.75
    ) +
    # True and recovered annual lines
    ggplot2::geom_line(
      data = df_long,
      ggplot2::aes(x = .data$year, y = .data$value,
                   colour = .data$Series, linetype = .data$Series),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = df_long,
      ggplot2::aes(x = .data$year, y = .data$value,
                   colour = .data$Series),
      size = 1.8, shape = 19
    ) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_linetype_manual(values = lty,  name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, T, by = max(1L, T %/% 10L)),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "Year",
      y        = "Signal value"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      plot.title       = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )

  invisible(p)
}


# ================================================================
# 5.  plot_lambda_heatmap()
# ================================================================

#' Plot RMSE and CCC Heatmaps Over the Lambda Grid
#'
#' @description
#' Creates a two-panel [ggplot2][ggplot2::ggplot2-package] heatmap from a
#' `"SimulationSuite"` object, showing how RMSE and CCC vary across the
#' \eqn{(\lambda_t, \lambda_s)} grid.  The optimal grid point is
#' highlighted with a cross.
#'
#' @param suite A `"SimulationSuite"` object returned by
#'   [run_simulation_suite()].
#' @param metric Character.  One of `"both"` (default), `"rmse"`, or
#'   `"ccc"`.
#'
#' @return A [ggplot2::ggplot] object (invisible).
#'
#' @examples
#' \dontrun{
#' plot_lambda_heatmap(suite)
#' plot_lambda_heatmap(suite, metric = "rmse")
#' }
#' @export
plot_lambda_heatmap <- function(suite, metric = c("both", "rmse", "ccc")) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot_lambda_heatmap().")
  metric <- match.arg(metric)

  r <- suite$results
  r$log_lambda_t <- log10(r$lambda_t)
  r$log_lambda_s <- log10(r$lambda_s)

  best_row <- suite$best$idx
  best_pt  <- data.frame(
    log_lambda_t = log10(suite$best$lambda_t),
    log_lambda_s = log10(suite$best$lambda_s)
  )

  # Helper: one heatmap tile panel
  .tile <- function(data, fill_var, fill_lab, low, high, best_pt,
                    is_reverse = FALSE) {
    ggplot2::ggplot(data,
                    ggplot2::aes(x = .data$log_lambda_t,
                                 y = .data$log_lambda_s,
                                 fill = .data[[fill_var]])) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
      ggplot2::geom_point(data = best_pt,
                          ggplot2::aes(x = .data$log_lambda_t,
                                       y = .data$log_lambda_s),
                          inherit.aes = FALSE,
                          shape = 3, size = 5, stroke = 1.8,
                          colour = "black") +
      ggplot2::scale_fill_gradient(
        low = low, high = high, name = fill_lab,
        guide = ggplot2::guide_colourbar(barwidth = 8, barheight = 0.6,
                                         title.position = "top",
                                         title.hjust = 0.5)
      ) +
      ggplot2::scale_x_continuous(
        name   = expression(log[10](lambda[t])),
        breaks = unique(r$log_lambda_t),
        labels = function(x) formatC(10^x, format = "g", digits = 2)
      ) +
      ggplot2::scale_y_continuous(
        name   = expression(log[10](lambda[s])),
        breaks = unique(r$log_lambda_s),
        labels = function(x) formatC(10^x, format = "g", digits = 2)
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position  = "bottom",
        panel.grid       = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }

  p_rmse <- .tile(r, "rmse", "RMSE", "#fee08b", "#d73027", best_pt) +
    ggplot2::ggtitle("RMSE (lower \u2190 better)")

  p_ccc  <- .tile(r, "ccc",  "CCC",  "#d73027", "#1a9850", best_pt) +
    ggplot2::ggtitle("CCC (higher \u2192 better)")

  if (metric == "rmse") return(invisible(p_rmse))
  if (metric == "ccc")  return(invisible(p_ccc))

  # Both panels side-by-side via patchwork (if available) or cowplot
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- patchwork::wrap_plots(p_rmse, p_ccc, ncol = 2L) +
      patchwork::plot_annotation(
        title    = "Hyperparameter grid search",
        subtitle = sprintf(
          "Best: \u03bb_t=%.4g, \u03bb_s=%.4g  (RMSE=%.4f, CCC=%.4f)",
          suite$best$lambda_t, suite$best$lambda_s,
          suite$best$rmse,     suite$best$ccc),
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(face = "bold", size = 13),
          plot.subtitle = ggplot2::element_text(colour = "grey40")
        )
      )
  } else {
    message("Install 'patchwork' for side-by-side panels; returning RMSE heatmap only.")
    p <- p_rmse
  }
  invisible(p)
}
