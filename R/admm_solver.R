## ============================================================
## R/admm_solver.R
## ADMM solver for the spatio-temporal generalised lasso
##
## Vectorisation convention (throughout this file):
##   x = as.vector(t(X))   where X is n  x  T  ("time-within-cell")
##   y = as.vector(t(Y))   where Y is n  x  S
##
## Kronecker structure:
##   M_kron  = I_n \u2297 M          (nS  x  nT)   observation operator
##   Dt_kron = I_n \u2297 D          (n(T-k)  x  nT) temporal penalty
##   Ds_kron = L_s \u2297 I_T        (nT  x  nT)   spatial penalty
## ============================================================

# ---- internal helpers ------------------------------------------

#' Element-wise soft-thresholding operator S_kappa(v)
#' @noRd
.soft_thresh <- function(v, kappa) {
  sign(v) * pmax(abs(v) - kappa, 0.0)
}

#' Build Kronecker products (pure sparse Matrix arithmetic)
#' Returns a list with M_kron, Dt_kron, Ds_kron.
#' @noRd
.build_kron <- function(M_sp, D_sp, Ls_sp, n, T_yr) {
  I_n <- Matrix::Diagonal(n)
  I_T <- Matrix::Diagonal(T_yr)
  list(
    M_kron  = Matrix::kronecker(I_n,   M_sp),   # (nS) \u00D7 (nT)
    Dt_kron = Matrix::kronecker(I_n,   D_sp),   # n(T-k) \u00D7 (nT)
    Ds_kron = Matrix::kronecker(Ls_sp, I_T)     # (nT)  \u00D7 (nT)
  )
}

# ================================================================
# 1.  R6 class:  ADMMSolver
# ================================================================

#' ADMM Solver for the Spatio-Temporal Generalised Lasso
#'
#' @description
#' Solves the Kronecker-structured generalised lasso problem
#'
#' \deqn{\min_X \;
#'   \tfrac{1}{2} \| W^{1/2} \bigl(\operatorname{vec}(Y)
#'       - (I_n \otimes M)\operatorname{vec}(X^\top)\bigr) \|_2^2
#'   + \lambda_t \|(I_n \otimes D)\operatorname{vec}(X^\top)\|_1
#'   + \lambda_s \|(L_s \otimes I_T)\operatorname{vec}(X^\top)\|_1}
#'
#' via the Alternating Direction Method of Multipliers (ADMM).
#'
#' @details
#' ## Algorithm (scaled dual form, Boyd et al. 2010, Sec. 3.1.1)
#'
#' Introducing auxiliary variables \eqn{z_t = D_t x} and \eqn{z_s = D_s x},
#' the augmented Lagrangian is split into three closed-form steps per
#' iteration:
#'
#' \describe{
#'   \item{x-update}{Solve one positive-definite linear system
#'     \eqn{A x^{k+1} = b^k} where
#'     \deqn{A = \tilde{M}^\top W \tilde{M}
#'               + \rho (D_t^\top D_t + D_s^\top D_s) + \delta I}
#'     is **precomputed and Cholesky-factorised once** (before the loop).
#'     Only the right-hand side \eqn{b^k} changes each iteration.}
#'   \item{z-updates}{Element-wise soft-thresholding (proximal operator
#'     of the \eqn{\ell_1} norm):
#'     \eqn{z_t^{k+1} = S_{\lambda_t/\rho}(D_t x^{k+1} + u_t^k)},
#'     similarly for \eqn{z_s}.}
#'   \item{u-updates}{Dual variable accumulation:
#'     \eqn{u_t^{k+1} = u_t^k + D_t x^{k+1} - z_t^{k+1}}.}
#' }
#'
#' ## Efficiency note
#'
#' Because \eqn{A} does **not** depend on \eqn{\lambda_t} or
#' \eqn{\lambda_s}, a single call to `$setup()` (which factors \eqn{A})
#' supports a full grid search over \eqn{(\lambda_t, \lambda_s)} by
#' re-running `$run()` only -- no re-factorisation needed.
#'
#' ## Convergence criterion
#'
#' Stopping uses the primal and dual residual norms (Boyd et al. 2010,
#' Sec. 3.3.1) with absolute tolerance `tol_abs` and relative tolerance
#' `tol_rel`.
#'
#' @references
#' Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2010).
#' *Distributed Optimization and Statistical Learning via the Alternating
#' Direction Method of Multipliers.* Foundations and Trends in Machine
#' Learning, 3(1), 1-122.
#'
#' @examples
#' \dontrun{
#' gen <- SyntheticDataGenerator$new(n_row=6, n_col=6, seed=1L)
#' gen$simulate()
#' D  <- build_temporal_penalty(gen$params$T_years, k = 2)
#' Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
#' W  <- build_precision_weights(matrix(0.5, 36, 16))
#'
#' solver <- ADMMSolver$new(gen$Y, gen$M, W, D, Ls,
#'                          lambda_t = 0.1, lambda_s = 0.05)
#' solver$setup()
#' solver$run(verbose = TRUE)
#' solver$print()
#' dim(solver$X_hat)   # 36 x 20
#' solver$compute_objective()
#' }
#'
#' @export
ADMMSolver <- R6::R6Class(
  classname = "ADMMSolver",

  public = list(

    # ---- inputs ------------------------------------------------
    #' @field Y n  x  S matrix of observed smoothed data.
    Y        = NULL,
    #' @field M S  x  T convolution matrix (sparse).
    M        = NULL,
    #' @field W (nS)  x  (nS) sparse diagonal precision matrix.
    W        = NULL,
    #' @field D (T-k)  x  T temporal difference operator (sparse).
    D        = NULL,
    #' @field Ls n  x  n spatial graph Laplacian (sparse, symmetric).
    Ls       = NULL,
    #' @field lambda_t Temporal regularisation parameter >= 0.
    lambda_t = NULL,
    #' @field lambda_s Spatial regularisation parameter >= 0.
    lambda_s = NULL,
    #' @field rho ADMM augmented-Lagrangian penalty parameter > 0.
    rho      = NULL,
    #' @field ridge Ridge added to system matrix diagonal for PD guarantee.
    ridge    = NULL,

    # ---- dimensions -------------------------------------------
    #' @field n_cells Number of spatial cells (n).
    n_cells  = NULL,
    #' @field T_years Number of latent annual time points (T).
    T_years  = NULL,
    #' @field n_periods Number of observed smoothed periods (S = ncol(Y)).
    n_periods = NULL,

    # ---- Kronecker operators (populated by setup()) ------------
    #' @field M_kron  (nS)  x  (nT) sparse matrix \eqn{I_n \otimes M}.
    M_kron   = NULL,
    #' @field Dt_kron n(T-k)  x  (nT) sparse matrix \eqn{I_n \otimes D}.
    Dt_kron  = NULL,
    #' @field Ds_kron (nT)  x  (nT) sparse matrix \eqn{L_s \otimes I_T}.
    Ds_kron  = NULL,
    #' @field A_sys (nT)  x  (nT) sparse SPD system matrix.
    A_sys    = NULL,

    # ---- results ----------------------------------------------
    #' @field X_hat n  x  T matrix of recovered latent annual signals.
    X_hat    = NULL,
    #' @field history `data.frame` with columns `iter`, `r_primal`, `r_dual`.
    history  = NULL,
    #' @field converged Logical; TRUE if stopping criteria were satisfied.
    converged = FALSE,

    # ---- private cache ----------------------------------------
    #' @field .chol_A CHMfactor (sparse Cholesky of `A_sys`). Populated by
    #'   `$setup()`.
    .chol_A  = NULL,
    #' @field .q_data nT-vector \eqn{M^\top W y} (fixed data term in RHS).
    #'   Populated by `$setup()`.
    .q_data  = NULL,

    # ===========================================================
    # initialize
    # ===========================================================
    #' @description
    #' Create a new `ADMMSolver`.
    #'
    #' @param Y         n  x  S numeric matrix of observed data.
    #' @param M         S  x  T convolution matrix.
    #' @param W         (nS)  x  (nS) sparse diagonal precision matrix;
    #'   `NULL` for unweighted LS; a positive scalar for homoscedastic LS.
    #' @param D         (T-k)  x  T temporal difference matrix from
    #'   [build_temporal_penalty()].
    #' @param Ls        n  x  n spatial Laplacian from
    #'   [build_spatial_laplacian()].
    #' @param lambda_t  Non-negative scalar temporal regularisation strength.
    #' @param lambda_s  Non-negative scalar spatial regularisation strength.
    #' @param rho       ADMM penalty parameter (default 1).  A good rule of
    #'   thumb is `rho ~= sqrt(lambda_t * lambda_s)` or `rho = 1` to start.
    #' @param ridge     Small positive scalar added to diag(A_sys) to
    #'   guarantee strict positive definiteness (default 1e-8).
    initialize = function(Y, M, W = NULL, D, Ls,
                          lambda_t,
                          lambda_s,
                          rho   = 1.0,
                          ridge = 1e-8) {

      .sp <- function(x)
        if (inherits(x, "Matrix")) x else Matrix::Matrix(x, sparse = TRUE)

      self$Y         <- if (is.matrix(Y)) Y else as.matrix(Y)
      self$M         <- .sp(M)
      self$D         <- .sp(D)
      self$Ls        <- .sp(Ls)
      self$lambda_t  <- lambda_t
      self$lambda_s  <- lambda_s
      self$rho       <- rho
      self$ridge     <- ridge
      self$n_cells   <- nrow(Y)
      self$n_periods <- ncol(Y)
      self$T_years   <- ncol(M)

      # W: NULL -> identity; scalar -> scaled identity; matrix -> as-is
      nS <- self$n_cells * self$n_periods
      self$W <- if (is.null(W)) {
        Matrix::Diagonal(nS, 1.0)
      } else if (length(W) == 1L) {
        Matrix::Diagonal(nS, as.numeric(W))
      } else {
        .sp(W)
      }

      # Dimension guards
      if (nrow(self$M) != self$n_periods)
        stop("nrow(M) must equal ncol(Y) (= ", self$n_periods, ").")
      if (nrow(self$W) != nS)
        stop("W must be a ", nS, " x ", nS, " matrix.")
      if (ncol(self$D) != self$T_years)
        stop("ncol(D) must equal ncol(M) (= ", self$T_years, ").")
      if (!all(dim(self$Ls) == self$n_cells))
        stop("Ls must be ", self$n_cells, " x ", self$n_cells, ".")
      if (lambda_t < 0 || lambda_s < 0)
        stop("lambda_t and lambda_s must be non-negative.")
      if (rho <= 0)
        stop("rho must be strictly positive.")

      invisible(self)
    },

    # ===========================================================
    # setup  -- build Kronecker products & factorise A_sys  ONCE
    # ===========================================================
    #' @description
    #' Precompute Kronecker operators, assemble the system matrix \eqn{A},
    #' and compute its sparse Cholesky factorisation.
    #'
    #' @details
    #' \eqn{A} is assembled as
    #' \deqn{A = \tilde{M}^\top W \tilde{M}
    #'           + \rho (D_t^\top D_t + D_s^\top D_s) + \delta I,}
    #' where \eqn{\delta} = `ridge`.  Because \eqn{A} is independent of
    #' \eqn{\lambda_t} and \eqn{\lambda_s}, this factorisation can be
    #' reused for any regularisation grid.
    #'
    #' @return `self` invisibly.
    setup = function() {
      n  <- self$n_cells
      T  <- self$T_years
      nT <- n * T

      # -- Kronecker operators -----------------------------------
      kron        <- .build_kron(self$M, self$D, self$Ls, n, T)
      self$M_kron  <- kron$M_kron
      self$Dt_kron <- kron$Dt_kron
      self$Ds_kron <- kron$Ds_kron

      # -- Fixed data term in RHS: q = M_kron^T W y -------------
      y_vec         <- as.numeric(t(self$Y))        # vec(Y^T), length nS
      MtW           <- Matrix::crossprod(self$M_kron, self$W)  # (nT \u00D7 nS) * W
      self$.q_data  <- as.numeric(MtW %*% y_vec)   # nT-vector

      # -- System matrix A_sys (assembled once, factorised once) -
      A_data    <- MtW %*% self$M_kron                       # M^T W M
      A_temp    <- Matrix::crossprod(self$Dt_kron)            # D_t^T D_t
      A_spat    <- Matrix::crossprod(self$Ds_kron)            # D_s^T D_s
      self$A_sys <- Matrix::forceSymmetric(
        A_data + self$rho * (A_temp + A_spat) +
          self$ridge * Matrix::Diagonal(nT)
      )

      # -- Sparse Cholesky with AMD permutation (the expensive step) -
      self$.chol_A <- Matrix::Cholesky(self$A_sys, perm = TRUE, LDL = FALSE)

      invisible(self)
    },

    # ===========================================================
    # run  -- ADMM loop
    # ===========================================================
    #' @description
    #' Execute the ADMM iterations.  Calls `setup()` automatically if not
    #' already invoked.
    #'
    #' @param max_iter   Positive integer. Maximum iterations (default 1000).
    #' @param tol_abs    Absolute convergence tolerance (default 1e-4).
    #' @param tol_rel    Relative convergence tolerance (default 1e-3).
    #' @param verbose    Logical. Print progress every 100 iterations
    #'   (default `FALSE`).
    #' @param warm_start Logical. If `TRUE` and `$X_hat` already exists,
    #'   initialise `x` from the previous solution (default `FALSE`).
    #'
    #' @return `self` invisibly.  Results stored in `$X_hat` and `$history`.
    run = function(max_iter    = 1000L,
                   tol_abs     = 1e-4,
                   tol_rel     = 1e-3,
                   verbose     = FALSE,
                   warm_start  = FALSE) {

      if (is.null(self$.chol_A)) self$setup()

      rho    <- self$rho
      lt_rho <- self$lambda_t / rho     # soft-threshold for temporal
      ls_rho <- self$lambda_s / rho     # soft-threshold for spatial

      Dt  <- self$Dt_kron
      Ds  <- self$Ds_kron
      DtT <- Matrix::t(Dt)
      DsT <- Matrix::t(Ds)

      nT   <- self$n_cells * self$T_years
      n_zt <- nrow(Dt)          # n * (T - k)
      n_zs <- nrow(Ds)          # n * T

      # -- initialise primal / dual variables --------------------
      x <- if (warm_start && !is.null(self$X_hat))
        as.numeric(t(self$X_hat))
      else
        numeric(nT)
      z_t <- numeric(n_zt)
      z_s <- numeric(n_zs)
      u_t <- numeric(n_zt)
      u_s <- numeric(n_zs)

      # -- pre-allocate convergence history ----------------------
      r_pri_hist  <- numeric(max_iter)
      r_dual_hist <- numeric(max_iter)

      converged <- FALSE
      n_iter    <- as.integer(max_iter)

      for (k in seq_len(max_iter)) {

        # -- x-update: one sparse triangular solve --------------
        b     <- self$.q_data +
                 rho * (as.numeric(DtT %*% (z_t - u_t)) +
                        as.numeric(DsT %*% (z_s - u_s)))
        x_new <- as.numeric(Matrix::solve(self$.chol_A, b))

        # -- z-updates: element-wise soft-thresholding ----------
        Dt_x    <- as.numeric(Dt %*% x_new)
        Ds_x    <- as.numeric(Ds %*% x_new)
        z_t_new <- .soft_thresh(Dt_x + u_t, lt_rho)
        z_s_new <- .soft_thresh(Ds_x + u_s, ls_rho)

        # -- u-updates: dual accumulation -----------------------
        r_t <- Dt_x - z_t_new        # primal residual (temporal)
        r_s <- Ds_x - z_s_new        # primal residual (spatial)
        u_t <- u_t + r_t
        u_s <- u_s + r_s

        # -- convergence diagnostics (Boyd et al. 2010 \u00A73.3.1) --
        r_primal <- sqrt(sum(r_t^2) + sum(r_s^2))

        # dual residuals: rho * D^T (z_new - z_old)
        s_t    <- rho * as.numeric(DtT %*% (z_t_new - z_t))
        s_s    <- rho * as.numeric(DsT %*% (z_s_new - z_s))
        r_dual <- sqrt(sum(s_t^2) + sum(s_s^2))

        # feasibility norms for relative stopping
        p_norm   <- sqrt(sum(Dt_x^2)   + sum(Ds_x^2))
        z_norm   <- sqrt(sum(z_t_new^2) + sum(z_s_new^2))
        u_norm   <- sqrt(
          sum(as.numeric(DtT %*% u_t)^2) +
          sum(as.numeric(DsT %*% u_s)^2)
        )

        eps_pri  <- sqrt(n_zt + n_zs) * tol_abs +
                    tol_rel * max(p_norm, z_norm)
        eps_dual <- sqrt(nT) * tol_abs +
                    tol_rel * rho * u_norm

        # commit z after computing dual residuals
        z_t <- z_t_new
        z_s <- z_s_new
        x   <- x_new

        r_pri_hist[k]  <- r_primal
        r_dual_hist[k] <- r_dual

        if (verbose && (k %% 100L == 0L || k == 1L))
          cat(sprintf("  [ADMM %4d]  r_pri=%.3e  r_dual=%.3e"
                      , k, r_primal, r_dual),
              sprintf("  eps_pri=%.3e  eps_dual=%.3e\n",
                      eps_pri, eps_dual))

        if (r_primal < eps_pri && r_dual < eps_dual) {
          converged <- TRUE
          n_iter    <- k
          if (verbose)
            cat(sprintf("  Converged at iteration %d.\n", k))
          break
        }
      }

      if (!converged)
        warning("ADMM did not converge in ", max_iter,
                " iterations.  Consider increasing max_iter or tuning rho.",
                call. = FALSE)

      # -- reconstruct X from x = vec(X^T) ---------------------
      # x is time-within-cell: x[((i-1)*T + t)] = X[i, t]
      # So: X = t(matrix(x, nrow = T, ncol = n))
      self$X_hat <- t(matrix(x, nrow = self$T_years, ncol = self$n_cells))
      if (!is.null(rownames(self$Y)))
        rownames(self$X_hat) <- rownames(self$Y)
      colnames(self$X_hat)   <- paste0("year_", seq_len(self$T_years))

      self$converged <- converged
      self$history   <- data.frame(
        iter     = seq_len(n_iter),
        r_primal = r_pri_hist[seq_len(n_iter)],
        r_dual   = r_dual_hist[seq_len(n_iter)]
      )

      invisible(self)
    },

    # ===========================================================
    # compute_objective
    # ===========================================================
    #' @description
    #' Evaluate the three-term objective at the current solution `$X_hat`.
    #'
    #' @return Named numeric vector with elements `data`, `temporal`,
    #'   `spatial`, and `total`.
    compute_objective = function() {
      if (is.null(self$X_hat))
        stop("Call run() first.")
      x         <- as.numeric(t(self$X_hat))
      y_vec     <- as.numeric(t(self$Y))
      r         <- y_vec - as.numeric(self$M_kron %*% x)
      data_term <- 0.5 * as.numeric(crossprod(r, as.numeric(self$W %*% r)))
      temp_term <- self$lambda_t * sum(abs(as.numeric(self$Dt_kron %*% x)))
      spat_term <- self$lambda_s * sum(abs(as.numeric(self$Ds_kron %*% x)))
      c(data     = data_term,
        temporal = temp_term,
        spatial  = spat_term,
        total    = data_term + temp_term + spat_term)
    },

    # ===========================================================
    # print
    # ===========================================================
    #' @description Print a concise solver summary.
    #' @param ... Unused; accepted for compatibility with the generic
    #'   `print()` method.
    print = function(...) {
      cat("\u2500\u2500 ADMMSolver (TWDeConv) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
      cat(sprintf("  Problem  : %d cells \u00d7 %d latent years"
                  , self$n_cells, self$T_years))
      cat(sprintf(" | %d observed periods\n", self$n_periods))
      cat(sprintf("  Params   : \u03bb_t=%.4g  \u03bb_s=%.4g  \u03c1=%.4g"
                  , self$lambda_t, self$lambda_s, self$rho))
      cat(sprintf("  ridge=%.2g\n", self$ridge))
      cat(sprintf("  Setup    : %s\n",
                  if (!is.null(self$.chol_A)) "done" else "pending"))
      if (!is.null(self$history)) {
        nr <- nrow(self$history)
        cat(sprintf("  Iters    : %d  |  Converged: %s\n",
                    nr, self$converged))
        cat(sprintf("  Final    : r_primal=%.3e  r_dual=%.3e\n",
                    self$history$r_primal[nr],
                    self$history$r_dual[nr]))
        obj <- self$compute_objective()
        cat(sprintf("  Objective: total=%.4g  (data=%.4g"
                    , obj["total"], obj["data"]))
        cat(sprintf("  temp=%.4g  spat=%.4g)\n",
                    obj["temporal"], obj["spatial"]))
      } else {
        cat("  Status   : not yet run\n")
      }
      cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
      invisible(self)
    }
  )
)


# ================================================================
# 2.  User-facing wrapper: solve_deconv()
# ================================================================

#' Solve the Spatio-Temporal Deconvolution Problem
#'
#' @description
#' Recovers the latent annual signal matrix \eqn{\hat{X}} from observed
#' multi-year moving-average data \eqn{Y} by solving the Kronecker-structured
#' generalised lasso
#'
#' \deqn{\min_X \;
#'   \tfrac{1}{2} \| W^{1/2} \bigl(\operatorname{vec}(Y)
#'     - (I_n \otimes M)\operatorname{vec}(X^\top)\bigr) \|_2^2
#'   + \lambda_t \|(I_n \otimes D)\operatorname{vec}(X^\top)\|_1
#'   + \lambda_s \|(L_s \otimes I_T)\operatorname{vec}(X^\top)\|_1}
#'
#' via ADMM with a precomputed sparse Cholesky factorisation.  See
#' [ADMMSolver] for full algorithmic details.
#'
#' @param Y            n  x  S numeric matrix of observed smoothed data.
#' @param M            S  x  T convolution matrix (from
#'   `SyntheticDataGenerator$M` or [build_temporal_penalty()]).
#' @param W            (nS)  x  (nS) sparse diagonal precision matrix from
#'   [build_precision_weights()].  Pass `NULL` for ordinary (unweighted)
#'   least squares; pass a positive scalar for homoscedastic WLS.
#' @param D            (T-k)  x  T temporal difference matrix from
#'   [build_temporal_penalty()].
#' @param Ls           n  x  n spatial graph Laplacian from
#'   [build_spatial_laplacian()].
#' @param lambda_t     Non-negative scalar.  Temporal regularisation
#'   strength (\eqn{\lambda_t}).
#' @param lambda_s     Non-negative scalar.  Spatial regularisation
#'   strength (\eqn{\lambda_s}).
#' @param rho          Positive scalar.  ADMM augmented-Lagrangian penalty
#'   parameter (default `1`).  A useful heuristic is
#'   `rho = sqrt(lambda_t * lambda_s)` or to set it near the expected
#'   signal-to-noise ratio.
#' @param ridge        Small positive ridge added to the system matrix
#'   diagonal to guarantee strict positive definiteness (default `1e-8`).
#' @param max_iter     Maximum ADMM iterations (default `1000`).
#' @param tol_abs      Absolute convergence tolerance (default `1e-4`).
#' @param tol_rel      Relative convergence tolerance (default `1e-3`).
#' @param verbose      Logical.  Print convergence progress (default `FALSE`).
#' @param return_solver Logical.  If `TRUE`, return the full [ADMMSolver]
#'   object (which carries `$X_hat`, `$history`, `$A_sys`, etc.) instead
#'   of just `$X_hat` (default `FALSE`).
#'
#' @return
#' * If `return_solver = FALSE` (default): an n  x  T numeric matrix
#'   \eqn{\hat{X}} with row names from `Y` and column names
#'   `year_1`, ..., `year_T`.
#' * If `return_solver = TRUE`: the [ADMMSolver] R6 object.
#'
#' @examples
#' \dontrun{
#' # -- Simulate ground-truth data --------------------------------
#' gen <- SyntheticDataGenerator$new(n_row = 8, n_col = 8,
#'                                   T_years = 20, seed = 42L)
#' gen$simulate()
#'
#' # -- Build penalty operators ----------------------------------
#' D  <- build_temporal_penalty(gen$params$T_years, k = 2)
#' Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
#' W  <- build_precision_weights(matrix(0.5, nrow = 64, ncol = 16))
#'
#' # -- Deconvolve -----------------------------------------------
#' X_hat <- solve_deconv(
#'   Y = gen$Y, M = gen$M, W = W, D = D, Ls = Ls,
#'   lambda_t = 0.1, lambda_s = 0.05,
#'   verbose  = TRUE
#' )
#' dim(X_hat)  # 64 x 20
#'
#' # -- Recovery quality -----------------------------------------
#' cat("RMSE:", sqrt(mean((X_hat - gen$X)^2)), "\n")
#' cat("Baseline (Y vs X):",
#'     sqrt(mean((gen$Y - gen$X[, 1:ncol(gen$Y)])^2)), "\n")
#'
#' # -- Return full solver object for diagnostics ----------------
#' solver <- solve_deconv(
#'   Y = gen$Y, M = gen$M, W = W, D = D, Ls = Ls,
#'   lambda_t = 0.1, lambda_s = 0.05,
#'   return_solver = TRUE
#' )
#' solver$compute_objective()
#' plot(solver$history$r_primal, type = "l", log = "y",
#'      ylab = "Primal residual", xlab = "Iteration")
#' }
#'
#' @seealso [ADMMSolver], [build_temporal_penalty()],
#'   [build_spatial_laplacian()], [build_precision_weights()],
#'   [SyntheticDataGenerator]
#' @export
solve_deconv <- function(Y, M, W = NULL, D, Ls,
                         lambda_t,
                         lambda_s,
                         rho           = 1.0,
                         ridge         = 1e-8,
                         max_iter      = 1000L,
                         tol_abs       = 1e-4,
                         tol_rel       = 1e-3,
                         verbose       = FALSE,
                         return_solver = FALSE) {
  solver <- ADMMSolver$new(
    Y = Y, M = M, W = W, D = D, Ls = Ls,
    lambda_t = lambda_t, lambda_s = lambda_s,
    rho = rho, ridge = ridge
  )
  solver$setup()
  solver$run(max_iter   = as.integer(max_iter),
             tol_abs    = tol_abs,
             tol_rel    = tol_rel,
             verbose    = verbose)

  if (return_solver) solver else solver$X_hat
}
