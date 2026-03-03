## ============================================================
## R/synthetic_data.R
## Synthetic spatio-temporal data generator for TWDeConv
## ============================================================

# ---- helpers (not exported) ------------------------------------

#' Build an exponential spatial covariance matrix
#'
#' @param coords Numeric matrix of centroid coordinates (n x 2).
#' @param range  Positive scalar. Spatial range parameter controlling
#'   how quickly correlation decays with distance.
#' @param nugget Small positive scalar added to the diagonal for numerical
#'   stability (default `1e-6`).
#' @return A symmetric positive-definite matrix of dimension n x n.
#' @noRd
.spatial_cov <- function(coords, range, nugget = 1e-6) {
  D     <- as.matrix(dist(coords))
  Sigma <- exp(-D / range)
  diag(Sigma) <- diag(Sigma) + nugget
  Sigma
}

#' Draw spatially correlated standard-normal samples
#'
#' @param n     Number of spatial locations.
#' @param k     Number of independent realisations to draw.
#' @param L_chol Lower-triangular Cholesky factor of the spatial covariance
#'   (`t(chol(Sigma))`).
#' @return Matrix n x k; each column is one spatial realisation.
#' @noRd
.draw_spatial <- function(n, k, L_chol) {
  Z <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  L_chol %*% Z
}

# ---- R6 class --------------------------------------------------

#' Synthetic Spatio-Temporal Data Generator
#'
#' @description
#' An [R6::R6Class] that simulates a latent annual raster signal, convolves it
#' with a 5-year moving-average kernel, and returns the "observed" blurred
#' data alongside the ground-truth objects needed for benchmarking
#' deconvolution methods.
#'
#' @details
#' ## Data-generating model
#'
#' Let \eqn{n = n_{\text{row}} \times n_{\text{col}}} be the number of grid
#' cells and \eqn{T} be the number of latent years.  The model is
#'
#' \deqn{Y = X M^\top + \varepsilon, \qquad
#'       \varepsilon_{i,s} \stackrel{\text{iid}}{\sim}
#'       \mathcal{N}(0, \sigma^2_\varepsilon)}
#'
#' where
#' * \eqn{X} (\eqn{n \times T}) - latent annual signal with spatial and
#'   temporal autocorrelation.
#' * \eqn{M} (\eqn{(T - w + 1) \times T}) - moving-average convolution matrix
#'   with window \eqn{w} (default 5).
#' * \eqn{Y} (\eqn{n \times (T - w + 1)}) - observed (smoothed) signal.
#'
#' ## Latent signal construction
#'
#' For each of the \eqn{T} years an *innovation* vector
#' \eqn{z_t \sim \mathcal{N}(0, \Sigma_s)} is drawn from a spatially
#' correlated Gaussian (exponential covariance kernel with range parameter
#' `spatial_range`).  The annual signals are then linked through an AR(1)
#' filter with coefficient `temporal_rho`:
#'
#' \deqn{X_{:,t} = \rho \, X_{:,t-1} +
#'       \sqrt{1 - \rho^2} \, z_t + \delta_t}
#'
#' where \eqn{\delta_t = \beta \cdot t} is a shared linear trend scaled by
#' `trend_sd`.
#'
#' @examples
#' \dontrun{
#' gen <- SyntheticDataGenerator$new(n_row = 10, n_col = 10,
#'                                   T_years = 20, seed = 2024)
#' gen$simulate()
#' gen$print()
#'
#' # Latent signal matrix  (100 cells x 20 years)
#' dim(gen$X)
#'
#' # Convolution matrix    (16 periods x 20 years)
#' dim(gen$M)
#'
#' # Observed data matrix  (100 cells x 16 periods)
#' dim(gen$Y)
#'
#' # sf object with X columns attached
#' sf_out <- gen$as_sf()
#' }
#'
#' @export
SyntheticDataGenerator <- R6::R6Class(
  classname = "SyntheticDataGenerator",

  # ---- public fields ------------------------------------------
  public = list(

    #' @field grid [`sf`][sf::sf] object containing the rectangular grid
    #'   polygons plus centroid coordinates.
    grid = NULL,

    #' @field X Numeric matrix \eqn{n \times T} of latent annual signals.
    #'   Rows are cells (`cell_1`, ...); columns are years (`year_1`, ...).
    X = NULL,

    #' @field M Numeric matrix \eqn{(T - w + 1) \times T}: the moving-average
    #'   convolution kernel.  Each row sums to 1.
    M = NULL,

    #' @field Y Numeric matrix \eqn{n \times (T - w + 1)} of observed
    #'   (smoothed + noisy) data.  Rows are cells; columns are ACS periods.
    Y = NULL,

    #' @field params Named list of all simulation hyperparameters.
    params = NULL,

    # ---- initialize -------------------------------------------
    #' @description
    #' Create a new `SyntheticDataGenerator`.
    #'
    #' @param n_row       Integer. Number of grid rows (default `10`).
    #' @param n_col       Integer. Number of grid columns (default `10`).
    #' @param T_years     Integer. Number of latent annual time points
    #'   (default `20`).
    #' @param window      Integer. Moving-average window length (default `5`
    #'   for ACS 5-year).
    #' @param spatial_range Positive numeric. Range of the exponential spatial
    #'   covariance kernel (in grid-cell units, default `3.0`).
    #' @param temporal_rho  Numeric in \eqn{(-1, 1)}.  AR(1) coefficient for
    #'   temporal autocorrelation (default `0.7`).
    #' @param trend_sd    Non-negative numeric.  Standard deviation of the
    #'   linear trend component (default `1.0`).
    #' @param noise_sd    Non-negative numeric.  Standard deviation of the
    #'   observation noise \eqn{\varepsilon} (default `0.5`).
    #' @param seed        Integer random seed for reproducibility (default `42`).
    initialize = function(n_row        = 10L,
                          n_col        = 10L,
                          T_years      = 20L,
                          window       = 5L,
                          spatial_range = 3.0,
                          temporal_rho = 0.7,
                          trend_sd     = 1.0,
                          noise_sd     = 0.5,
                          seed         = 42L) {
      stopifnot(
        is.numeric(spatial_range), spatial_range > 0,
        is.numeric(temporal_rho),  abs(temporal_rho) < 1,
        is.numeric(noise_sd),      noise_sd >= 0,
        is.numeric(trend_sd),      trend_sd >= 0,
        as.integer(T_years) > as.integer(window)
      )
      self$params <- list(
        n_row         = as.integer(n_row),
        n_col         = as.integer(n_col),
        n_cells       = as.integer(n_row) * as.integer(n_col),
        T_years       = as.integer(T_years),
        window        = as.integer(window),
        n_periods     = as.integer(T_years) - as.integer(window) + 1L,
        spatial_range = spatial_range,
        temporal_rho  = temporal_rho,
        trend_sd      = trend_sd,
        noise_sd      = noise_sd,
        seed          = as.integer(seed)
      )
      set.seed(self$params$seed)
      invisible(self)
    },

    # ---- generate_grid ----------------------------------------
    #' @description
    #' Build a regular rectangular grid using [sf::st_make_grid()].
    #'
    #' Cells are unit squares; the bounding box spans
    #' \eqn{[0, n_{\text{col}}] \times [0, n_{\text{row}}]}.
    #' Centroid coordinates are stored as extra columns `centroid_x` and
    #' `centroid_y` for downstream distance calculations.
    #'
    #' @return `self` (invisibly) -- modifies `self$grid` in place.
    generate_grid = function() {
      p    <- self$params
      bbox <- sf::st_bbox(c(xmin = 0, ymin = 0,
                            xmax = p$n_col, ymax = p$n_row),
                          crs = sf::NA_crs_)
      polys <- sf::st_make_grid(sf::st_as_sfc(bbox),
                                n      = c(p$n_col, p$n_row),
                                square = TRUE)
      grid  <- sf::st_sf(cell_id = seq_along(polys),
                         geometry = polys)
      ctr   <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(grid)))
      grid$centroid_x <- ctr[, 1]
      grid$centroid_y <- ctr[, 2]
      self$grid <- grid
      invisible(self)
    },

    # ---- generate_latent_signal --------------------------------
    #' @description
    #' Simulate the latent annual signal matrix **X** (\eqn{n \times T}).
    #'
    #' @details
    #' 1. Compute an \eqn{n \times n} exponential covariance matrix
    #'    \eqn{\Sigma_s} from centroid distances.
    #' 2. Obtain the lower Cholesky factor \eqn{L = \text{chol}(\Sigma_s)^\top}.
    #' 3. At each year \eqn{t}, draw a spatially correlated innovation
    #'    \eqn{z_t \sim \mathcal{N}(0, \Sigma_s)} via \eqn{L z}, \eqn{z \sim
    #'    \mathcal{N}(0, I)}.
    #' 4. Apply an AR(1) filter across years with coefficient `temporal_rho`.
    #' 5. Superimpose a cell-shared linear trend scaled by `trend_sd`.
    #'
    #' Calls `generate_grid()` automatically if `self$grid` is `NULL`.
    #'
    #' @return `self` (invisibly) -- modifies `self$X` in place.
    generate_latent_signal = function() {
      p <- self$params
      if (is.null(self$grid)) self$generate_grid()

      # --- Spatial covariance & Cholesky factor -------------------
      coords <- cbind(self$grid$centroid_x, self$grid$centroid_y)
      Sigma_s <- .spatial_cov(coords, p$spatial_range)
      L_chol  <- t(chol(Sigma_s))          # lower-triangular factor

      # --- Draw spatially correlated innovations (n x T) ----------
      Z_innov <- .draw_spatial(p$n_cells, p$T_years, L_chol)

      # --- Shared linear trend ------------------------------------
      t_seq   <- seq_len(p$T_years)
      trend   <- p$trend_sd * (t_seq / p$T_years)  # normalised 0\u2192trend_sd

      # --- AR(1) temporal filter ----------------------------------
      rho <- p$temporal_rho
      X   <- matrix(0, nrow = p$n_cells, ncol = p$T_years)
      X[, 1L] <- Z_innov[, 1L]
      for (t in seq.int(2L, p$T_years)) {
        X[, t] <- rho * X[, t - 1L] +
                  sqrt(1 - rho^2) * Z_innov[, t]
      }

      # --- Add linear trend (broadcast across rows) ---------------
      X <- X + matrix(trend, nrow = p$n_cells, ncol = p$T_years,
                      byrow = TRUE)

      dimnames(X) <- list(
        paste0("cell_",  self$grid$cell_id),
        paste0("year_",  t_seq)
      )
      self$X <- X
      invisible(self)
    },

    # ---- build_convolution_matrix ------------------------------
    #' @description
    #' Construct the moving-average convolution matrix **M**.
    #'
    #' @details
    #' **M** has dimensions \eqn{(T - w + 1) \times T}.  Row \eqn{s}
    #' contains \eqn{1/w} in columns \eqn{s, s+1, \ldots, s+w-1} and
    #' zero elsewhere, encoding the \eqn{s}-th \eqn{w}-year average:
    #' \deqn{M_{s, \cdot} = \frac{1}{w}
    #'       \bigl[\underbrace{0,\ldots,0}_{s-1},
    #'             \underbrace{1,\ldots,1}_{w},
    #'             \underbrace{0,\ldots,0}_{T-s-w+1}\bigr].}
    #'
    #' @return `self` (invisibly) -- modifies `self$M` in place.
    build_convolution_matrix = function() {
      p   <- self$params
      T   <- p$T_years
      w   <- p$window
      n_p <- p$n_periods                    # T - w + 1

      M <- matrix(0, nrow = n_p, ncol = T)
      for (s in seq_len(n_p)) {
        M[s, s:(s + w - 1L)] <- 1 / w
      }
      dimnames(M) <- list(
        paste0("period_", seq_len(n_p)),
        paste0("year_",   seq_len(T))
      )
      self$M <- M
      invisible(self)
    },

    # ---- generate_observed_data --------------------------------
    #' @description
    #' Generate the observed data matrix **Y** (\eqn{n \times (T-w+1)}).
    #'
    #' @details
    #' \deqn{Y = X M^\top + \varepsilon, \quad
    #'       \varepsilon_{i,s} \stackrel{\text{iid}}{\sim}
    #'       \mathcal{N}(0, \sigma^2_\varepsilon).}
    #'
    #' Observation noise `noise_sd` is intentionally kept mild relative to
    #' the signal so that deconvolution methods have a meaningful but
    #' tractable target.
    #'
    #' Calls `generate_latent_signal()` and `build_convolution_matrix()`
    #' automatically if they have not yet been called.
    #'
    #' @return `self` (invisibly) -- modifies `self$Y` in place.
    generate_observed_data = function() {
      p <- self$params
      if (is.null(self$X)) self$generate_latent_signal()
      if (is.null(self$M)) self$build_convolution_matrix()

      # Y = X %*% t(M)  ->  (n x T) %*% (T x n_periods)  =  (n x n_periods)
      Y_true  <- self$X %*% t(self$M)

      epsilon <- matrix(stats::rnorm(p$n_cells * p$n_periods,
                                     sd = p$noise_sd),
                        nrow = p$n_cells, ncol = p$n_periods)

      self$Y <- Y_true + epsilon
      dimnames(self$Y) <- list(rownames(self$X), rownames(self$M))
      invisible(self)
    },

    # ---- simulate (pipeline) -----------------------------------
    #' @description
    #' Run the full four-step simulation pipeline:
    #' `generate_grid()` -> `generate_latent_signal()` ->
    #' `build_convolution_matrix()` -> `generate_observed_data()`.
    #'
    #' @return `self` (invisibly).
    simulate = function() {
      self$generate_grid()
      self$generate_latent_signal()
      self$build_convolution_matrix()
      self$generate_observed_data()
      invisible(self)
    },

    # ---- as_sf --------------------------------------------------
    #' @description
    #' Return an [`sf`][sf::sf] object with the latent annual signal columns
    #' (`year_1`, ..., `year_T`) appended to the grid geometries.
    #'
    #' Useful for spatial visualisation with **ggplot2** + **geom_sf**.
    #'
    #' @return An `sf` data frame with `n_cells` rows and `T_years + 2`
    #'   columns (`cell_id`, `year_*`, `geometry`).
    as_sf = function() {
      if (is.null(self$X))
        stop("Call simulate() or generate_latent_signal() first.")
      meta <- sf::st_drop_geometry(self$grid)[, "cell_id", drop = FALSE]
      sf::st_sf(cbind(meta, as.data.frame(self$X)),
                geometry = sf::st_geometry(self$grid))
    },

    #' @description
    #' Return a long-format `data.frame` of the observed data **Y** for
    #' easy plotting with **ggplot2**.
    #'
    #' @return A `data.frame` with columns `cell_id`, `period`, and `value`.
    observed_long = function() {
      if (is.null(self$Y))
        stop("Call simulate() or generate_observed_data() first.")
      df        <- as.data.frame(self$Y)
      df$cell_id <- rownames(self$Y)
      long      <- stats::reshape(df,
                                  varying   = setdiff(names(df), "cell_id"),
                                  v.names   = "value",
                                  timevar   = "period",
                                  times     = colnames(self$Y),
                                  direction = "long")
      rownames(long) <- NULL
      long[, c("cell_id", "period", "value")]
    },

    # ---- print --------------------------------------------------
    #' @description Print a concise summary of the generator state.
    #' @param ... Unused; accepted for compatibility with the generic
    #'   `print()` method.
    print = function(...) {
      p    <- self$params
      fmt  <- function(mat) if (!is.null(mat))
                paste(dim(mat), collapse = " \u00d7 ") else "\u2014"
      cat("\u2500\u2500 SyntheticDataGenerator (TWDeConv) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
      cat(sprintf("  Grid     : %d \u00d7 %d = %d cells\n",
                  p$n_row, p$n_col, p$n_cells))
      cat(sprintf("  Window   : %d-year moving average\n", p$window))
      cat(sprintf("  Latent T : %d years  |  Observed periods: %d\n",
                  p$T_years, p$n_periods))
      cat(sprintf("  X  %-5s : %s  (latent annual signal)\n",
                  "", fmt(self$X)))
      cat(sprintf("  M  %-5s : %s  (convolution matrix)\n",
                  "", fmt(self$M)))
      cat(sprintf("  Y  %-5s : %s  (observed smoothed data)\n",
                  "", fmt(self$Y)))
      cat(sprintf("  Params   : spatial_range=%.1f  temporal_rho=%.2f\n",
                  p$spatial_range, p$temporal_rho))
      cat(sprintf("             trend_sd=%.1f  noise_sd=%.1f  seed=%d\n",
                  p$trend_sd, p$noise_sd, p$seed))
      cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
      invisible(self)
    }
  )  # end public
)  # end R6Class
