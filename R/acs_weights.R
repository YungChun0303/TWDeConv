## ============================================================
## R/acs_weights.R
## Precision weight matrix from ACS replicate estimates
##
## Methodology
##   The ACS Successive Difference Replication (SDR) variance estimator
##   computes the variance of a tract-level estimate theta^_i as:
##     V(theta^_i) = (4/80) . Sigma_{r=1}^{80} (theta^_{r,i} - theta^_i)^2
##
##   We express this using a formal survey replicate design
##   (srvyr::as_survey_rep / survey::svrepdesign) so that downstream
##   code can treat the variance estimates as coming from an auditable,
##   reproducible survey-statistical framework.
##
##   For each tract i the design is built by converting the 80
##   pre-computed replicate ESTIMATES (theta^_{r,i}) into replicate FACTORS:
##     repwt[i, r] = theta^_{r,i} / theta^_i
##   With base weight 1 and combined.weights = FALSE, the domain total
##   in replicate r is repwt[i,r] . theta^_i = theta^_{r,i}, so the SDR
##   variance over the single-observation domain {i} reduces to the
##   textbook formula above.
##
##   When VRE replicates are unavailable the function falls back to
##   MOE-derived SE already stored in the ACSData object.
## ============================================================


# ================================================================
# build_acs_precision_weights()   (exported)
# ================================================================

#' Build a Precision Weight Matrix from an ACSData Object
#'
#' @description
#' Converts ACS 5-year standard errors (either from the 80 Variance
#' Replicate Estimates via SDR or from the published margin of error)
#' into the diagonal inverse-variance weight matrix \eqn{W} used in the
#' Weighted Least Squares deconvolution objective.
#'
#' When the **srvyr** and **survey** packages are installed and the
#' `ACSData` object contains VRE replicate columns
#' (`vre_source == "sdr"`), the function constructs a formal
#' *successive-difference replicate* survey design via
#' `srvyr::as_survey_rep()` and computes per-tract variances using
#' domain estimation (`survey::svyby()`).
#' Otherwise it uses the pre-computed SE already stored in
#' `acs_data$se` (which is the SDR SE when replicates were available
#' or MOE-derived SE otherwise).
#'
#' @details
#' ## Survey replicate design
#'
#' For each selected variable the function builds an \eqn{n}-row data
#' frame (one row per tract) and creates a replicate design with
#' 80 replicate weight columns equal to
#' \eqn{\hat{\theta}_{r,i} / \hat{\theta}_i}.
#' With base weight 1, the domain total in replicate \eqn{r} for
#' tract \eqn{i} is exactly \eqn{\hat{\theta}_{r,i}}, so
#' `survey::svyby(~estimate, ~GEOID, design, svytotal)` returns the
#' SDR standard error
#' \eqn{\sqrt{4/80 \cdot \sum_r (\hat{\theta}_{r,i} -
#' \hat{\theta}_i)^2}}.
#'
#' Tracts with a zero or missing base estimate cannot form a replicate
#' ratio and fall back to MOE-derived SE.
#'
#' ## Output ordering
#'
#' The diagonal of \eqn{W} follows **survey-within-cell** ordering,
#' matching the ADMM \eqn{y}-vector convention
#' (\eqn{y = \operatorname{vec}(Y^\top)}):
#' \deqn{W_{(i-1)p+j,\,(i-1)p+j} = 1 / \hat{V}(\hat{\theta}_{i,j})}
#' where \eqn{i} indexes tracts and \eqn{j} indexes variables.
#'
#' @param acs_data An object of class `"ACSData"` returned by
#'   [fetch_acs_vre_data()].
#' @param variables Character vector of variable names (must match
#'   `acs_data$metadata$var_names`). `NULL` (default) uses all
#'   variables.
#' @param use_srvyr Logical. If `TRUE` (default), use
#'   `srvyr::as_survey_rep()` when both **srvyr** / **survey** are
#'   installed and VRE replicates are available.  Set to `FALSE` to
#'   always use the pre-computed SE from `acs_data$se`.
#' @param zero_var_action Character. How to handle zero or non-finite
#'   variance: `"max_weight"` (default) caps precision at the maximum
#'   observed value; `"drop"` sets precision to 0 (observation
#'   excluded from WLS); `"error"` stops with an error message.
#' @param quiet Logical. Suppress progress messages (default `FALSE`).
#'
#' @return A sparse diagonal matrix of class
#'   [`ddiMatrix`][Matrix::ddiMatrix-class] with dimension
#'   \eqn{(n \cdot p) \times (n \cdot p)}, where \eqn{n} is the
#'   number of tracts and \eqn{p} the number of selected variables.
#'   Row/column names are `"GEOID:variable"` in survey-within-cell
#'   order.  Pass directly to [solve_deconv()] as the `W` argument.
#'
#' @examples
#' \dontrun{
#' acs <- fetch_acs_vre_data(
#'   state     = "IL",
#'   county    = "Cook",
#'   year_end  = 2022,
#'   variables = c(med_inc = "B19013_001"),
#'   vre       = TRUE
#' )
#'
#' W <- build_acs_precision_weights(acs)
#' dim(W)   # n_tracts  x  n_tracts (one variable)
#'
#' # Plug into deconvolution pipeline
#' Ls <- build_spatial_laplacian(acs$sf)
#' D  <- build_temporal_penalty(T = 10, k = 2)
#' # (stack Y across ACS windows first -- see fetch_acs_vre_data docs)
#' }
#'
#' @seealso [fetch_acs_vre_data()], [build_precision_weights()],
#'   [compute_sdr_se()], [solve_deconv()]
#' @export
build_acs_precision_weights <- function(acs_data,
                                         variables       = NULL,
                                         use_srvyr       = TRUE,
                                         zero_var_action = c("max_weight",
                                                             "drop",
                                                             "error"),
                                         quiet           = FALSE) {

  # -- 0. Validate input --------------------------------------------
  if (!inherits(acs_data, "ACSData"))
    stop("`acs_data` must be an ACSData object from fetch_acs_vre_data().")

  zero_var_action <- match.arg(zero_var_action)

  m         <- acs_data$metadata
  est_mat   <- acs_data$estimates
  se_mat    <- acs_data$se
  reps_list <- acs_data$replicates
  n_tracts  <- nrow(est_mat)

  # -- 1. Variable selection ----------------------------------------
  var_sel <- if (is.null(variables)) {
    m$var_names
  } else {
    v <- intersect(as.character(variables), m$var_names)
    if (length(v) == 0L)
      stop("None of the requested variables found in ACSData. ",
           "Available: ", paste(m$var_names, collapse = ", "))
    v
  }
  p <- length(var_sel)

  # -- 2. Check for srvyr / survey ----------------------------------
  has_survey_pkgs <- use_srvyr &&
    requireNamespace("srvyr",   quietly = TRUE) &&
    requireNamespace("survey",  quietly = TRUE)

  # -- 3. Compute variance per (tract  x  variable) -------------------
  # var_mat[i, j] = estimated variance for tract i, variable j
  var_mat <- matrix(NA_real_, nrow = n_tracts, ncol = p,
                    dimnames = list(rownames(est_mat), var_sel))

  for (j in seq_len(p)) {
    vj    <- var_sel[j]
    est_j <- est_mat[, vj]
    se_j  <- se_mat[, vj]
    rep_j <- reps_list[[vj]]

    use_design <- !is.null(rep_j) &&
                  ncol(rep_j) == 80L &&
                  has_survey_pkgs

    if (use_design) {
      # -- Path A: formal survey replicate design (srvyr + survey) --
      # Replicate factors: repwt[i,r] = R_ir / y_i
      # => svytotal for domain {i} in replicate r = repwt[i,r]*y_i = R_ir
      # => SDR var_i = (4/80) * sum_r (R_ir - y_i)^2
      if (!quiet)
        message(sprintf("  [%s] Using survey replicate design (SDR).", vj))

      var_j <- se_j^2    # initialise with MOE-derived fallback for all rows
      valid <- is.finite(est_j) & est_j != 0

      if (sum(valid) > 0L) {
        est_v    <- est_j[valid]
        rep_v    <- rep_j[valid, , drop = FALSE]
        geoids_v <- rownames(est_mat)[valid]

        # Build replicate ratio matrix (n_valid  x  80)
        rep_ratios <- sweep(rep_v, 1L, est_v, FUN = "/")
        rep_ratios[!is.finite(rep_ratios)] <- 1   # guard for edge cases

        # Construct design data frame with embedded replicate columns
        rep_col_nms <- paste0(".R", seq_len(80L))
        df_design   <- data.frame(
          GEOID       = geoids_v,
          estimate    = est_v,
          base_weight = 1,
          stringsAsFactors = FALSE
        )
        df_design[, rep_col_nms] <- as.data.frame(rep_ratios)

        # Create survey replicate design via srvyr
        design <- srvyr::as_survey_rep(
          .data            = df_design,
          repweights       = rep_col_nms,
          weights          = "base_weight",
          type             = "successive-diff",
          scale            = 4 / 80,
          combined.weights = FALSE
        )

        # Per-tract domain SE via survey::svyby
        # Each GEOID is a single-observation domain, so:
        #   svytotal(~estimate)_domain_i = y_i  (base)
        #   svytotal(~estimate)_domain_i_rep_r  = R_ir
        #   SE_i = sqrt((4/80) * sum_r (R_ir - y_i)^2)
        by_tract <- survey::svyby(
          formula = ~estimate,
          by      = ~GEOID,
          design  = design,
          FUN     = survey::svytotal,
          covmat  = FALSE
        )

        # SE column: survey names it "se.estimate" (or "se" if single var)
        se_col_nm <- grep("^se", names(by_tract), value = TRUE)[1L]
        idx       <- match(by_tract[["GEOID"]], rownames(est_mat))
        var_j[idx] <- by_tract[[se_col_nm]]^2
      }

      var_mat[, j] <- var_j

    } else {
      # -- Path B: use pre-computed SE from ACSData ------------------
      if (!quiet) {
        src <- if (m$vre_source == "sdr" && !is.null(rep_j))
          "pre-computed SDR" else "MOE-derived"
        message(sprintf("  [%s] Using %s SE.", vj, src))
      }
      var_mat[, j] <- se_j^2
    }
  }

  # -- 4. Flatten in survey-within-cell order -----------------------
  # Matches the ADMM y-vector: y = as.vector(t(Y)), so
  # element (i-1)*p+j corresponds to tract i, variable j.
  # as.vector(t(var_mat)) where var_mat is n x p gives column-major
  # of p x n = (tract1,v1),(tract1,v2),...,(tract2,v1),... OK
  var_vec <- as.vector(t(var_mat))    # length n_tracts * p

  # -- 5. Handle zero / non-finite variances -----------------------
  finite_pos <- is.finite(var_vec) & var_vec > 0

  if (!all(finite_pos)) {
    n_bad <- sum(!finite_pos)
    msg   <- sprintf(
      "%d tract-variable combination(s) have zero/NA/non-positive variance.",
      n_bad
    )
    if (zero_var_action == "error") {
      stop(msg)
    } else if (zero_var_action == "drop") {
      warning(msg, " Precision set to 0 (excluded from WLS).", call. = FALSE)
      var_vec[!finite_pos] <- Inf    # 1/Inf = 0
    } else {   # "max_weight"
      warning(msg, " Precision capped at max observed value.", call. = FALSE)
      max_prec             <- max(1 / var_vec[finite_pos], na.rm = TRUE)
      var_vec[!finite_pos] <- 1 / max_prec
    }
  }

  prec_vec <- ifelse(is.infinite(var_vec), 0, 1 / var_vec)

  # -- 6. Build and name the diagonal matrix -----------------------
  W <- Matrix::Diagonal(x = prec_vec)

  geoids <- rownames(est_mat)
  # Names in survey-within-cell order: (tract1:v1),(tract1:v2),...
  # as.vector(t(outer(geoids, var_sel))) =
  #   column-major of p x n = [g1:v1,g1:v2,...,g1:vp, g2:v1,...,gn:vp]
  nms <- as.vector(t(outer(geoids, var_sel, paste, sep = ":")))
  dimnames(W) <- list(nms, nms)

  W
}
