## ============================================================
## R/penalty_matrices.R
## Penalty operators for the spatio-temporal generalised lasso
## ============================================================

# ---- internal helpers (not exported) ---------------------------

#' First-order difference matrix of size (n-1) x n
#' @noRd
.D1 <- function(n) {
  n <- as.integer(n)
  Matrix::sparseMatrix(
    i    = rep(seq_len(n - 1L), each = 2L),
    j    = c(rbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)),
    x    = rep(c(-1, 1), times = n - 1L),
    dims = c(n - 1L, n)
  )
}

#' Convert an spdep neighbour list to a sparse binary adjacency matrix
#' @noRd
.nb_to_sparse_adj <- function(nb, n) {
  pairs <- do.call(rbind, lapply(seq_len(n), function(i) {
    js <- nb[[i]]
    if (identical(js, 0L)) return(NULL)   # isolated cell
    cbind(i = i, j = js)
  }))
  if (is.null(pairs)) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                dims = c(n, n)))
  }
  Matrix::sparseMatrix(
    i    = pairs[, "i"],
    j    = pairs[, "j"],
    x    = 1,
    dims = c(n, n)
  )
}

# ================================================================
# 1.  Temporal difference penalty matrix D
# ================================================================

#' Build the k-th order temporal difference penalty matrix
#'
#' @description
#' Constructs the sparse matrix \eqn{D_k} that is the standard penalty
#' operator for 1-D trend-filtering / fused-lasso problems on \eqn{T}
#' equidistant time points.
#'
#' @details
#' The first-order difference matrix \eqn{D_1} of size \eqn{(T-1) \times T}
#' has the stencil \eqn{[-1,\; 1]} shifted along the diagonal:
#' \deqn{(D_1 x)_t = x_{t+1} - x_t, \quad t = 1,\ldots,T-1.}
#' Higher-order matrices are formed recursively:
#' \deqn{D_k = D_1^{(T-k+1)} \, D_{k-1},}
#' so that \eqn{D_k} has dimensions \eqn{(T-k) \times T}.
#'
#' Common uses:
#' \describe{
#'   \item{k = 1}{Fused lasso -- penalises \eqn{\sum_t |x_{t+1}-x_t|}.}
#'   \item{k = 2}{Trend filtering -- penalises \eqn{\sum_t |x_{t+1}-2x_t+x_{t-1}|};
#'     the discrete analogue of a second-derivative penalty.}
#' }
#'
#' @param T Positive integer. Number of latent time periods.
#' @param k Positive integer. Differencing order; must satisfy `k < T`.
#'   Default `1`.
#'
#' @return A sparse matrix of class [`dgCMatrix`][Matrix::dgCMatrix-class]
#'   with dimensions \eqn{(T-k) \times T}.  Row names are `"d<k>_1"`, ...,
#'   `"d<k>_<T-k>"`; column names are `"year_1"`, ..., `"year_T"`.
#'
#' @examples
#' D1 <- build_temporal_penalty(T = 10, k = 1)
#' dim(D1)          # 9 x 10
#' as.matrix(D1[1:4, 1:5])
#'
#' D2 <- build_temporal_penalty(T = 10, k = 2)
#' dim(D2)          # 8 x 10
#' as.matrix(D2[1:4, 1:6])
#'
#' @export
build_temporal_penalty <- function(T, k = 1L) {
  T <- as.integer(T)
  k <- as.integer(k)
  if (k < 1L)  stop("`k` must be >= 1.")
  if (T <= k)  stop("`T` must be strictly greater than `k`.")

  D <- .D1(T)                          # (T-1) x T

  if (k >= 2L) {
    for (ord in seq_len(k - 1L)) {
      D <- .D1(T - ord) %*% D          # shrinks one row each step
    }
  }

  dimnames(D) <- list(
    paste0("d", k, "_", seq_len(nrow(D))),
    paste0("year_", seq_len(T))
  )
  D
}


# ================================================================
# 2.  Spatial graph Laplacian L_s
# ================================================================

#' Build the spatial graph Laplacian from an sf grid
#'
#' @description
#' Given an [`sf`][sf::sf] object of polygon geometries, computes the
#' contiguity graph among cells using [spdep::poly2nb()] and returns the
#' combinatorial graph Laplacian
#' \deqn{L_s = D_s - A,}
#' where \eqn{A} is the binary adjacency matrix and
#' \eqn{D_s = \operatorname{diag}(A \mathbf{1})} is the degree matrix.
#'
#' @details
#' **Contiguity rule.**  `queen = TRUE` (default) uses queen contiguity
#' (shared edge *or* vertex), giving up to 8 neighbours on a regular grid.
#' `queen = FALSE` uses rook contiguity (shared edge only), giving up to 4
#' neighbours -- the standard discrete Laplacian on a lattice.
#'
#' **Properties of \eqn{L_s}.**
#' \eqn{L_s} is symmetric positive semi-definite.  Its null space is spanned
#' by the constant vector, so the smallest eigenvalue is 0 (one per connected
#' component). The penalty \eqn{x^\top L_s x = \sum_{(i,j) \in \mathcal{E}}
#' (x_i - x_j)^2} discourages spatially non-smooth solutions.
#'
#' **Isolated cells.** Cells with no neighbours (e.g., non-contiguous islands
#' in a cropped grid) receive a zero row/column and are flagged with a warning.
#'
#' @param sf_obj An [`sf`][sf::sf] object whose geometries are polygons.
#'   Typically the `$grid` field from a [SyntheticDataGenerator] instance.
#' @param queen Logical. Use queen contiguity? Default `TRUE`.
#'
#' @return A sparse symmetric matrix of class
#'   [`dsCMatrix`][Matrix::dsCMatrix-class] with dimensions
#'   \eqn{n \times n} (one row/column per cell).  Row and column names
#'   match `sf_obj$cell_id` when present, otherwise `"cell_1"`, ....
#'
#' @examples
#' \dontrun{
#' gen <- SyntheticDataGenerator$new(n_row = 5, n_col = 4)
#' gen$generate_grid()
#'
#' Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
#' dim(Ls)       # 20 x 20
#' # Laplacian property: row sums are zero
#' all(abs(Matrix::rowSums(Ls)) < 1e-12)
#' }
#'
#' @export
build_spatial_laplacian <- function(sf_obj, queen = TRUE) {
  if (!inherits(sf_obj, "sf"))
    stop("`sf_obj` must be an sf object.")

  n  <- nrow(sf_obj)
  nb <- spdep::poly2nb(sf_obj, queen = queen)

  # Identify isolated cells
  n_isolated <- sum(vapply(nb, function(x) identical(x, 0L), logical(1)))
  if (n_isolated > 0L)
    warning(n_isolated, " isolated cell(s) with no neighbours detected.",
            call. = FALSE)

  # Binary adjacency matrix (sparse)
  A <- .nb_to_sparse_adj(nb, n)

  # Degree matrix (diagonal)
  deg <- Matrix::rowSums(A)
  Ds  <- Matrix::Diagonal(n = n, x = deg)

  # Laplacian: symmetric, so store as dsCMatrix
  Ls <- Matrix::forceSymmetric(Ds - A)

  # Dimnames: prefer cell_id column if present
  cell_nms <- if ("cell_id" %in% names(sf_obj))
    paste0("cell_", sf_obj$cell_id)
  else
    paste0("cell_", seq_len(n))
  dimnames(Ls) <- list(cell_nms, cell_nms)

  Ls
}


# ================================================================
# 3.  Empirical Spatial Laplacian from ACS tract geometries
# ================================================================

#' Build the Empirical Spatial Graph Laplacian from ACS Census Tract Geometries
#'
#' @description
#' Given an [`sf`][sf::sf] object of Census tract polygon geometries (typically
#' the `$sf` component of an `ACSData` object returned by
#' [fetch_acs_vre_data()]), computes the Queen contiguity graph among tracts
#' using [spdep::poly2nb()] and returns the combinatorial graph Laplacian
#' \deqn{L_s = D_s - A,}
#' where \eqn{A} is the binary adjacency matrix and
#' \eqn{D_s = \operatorname{diag}(A \mathbf{1})} is the degree matrix.
#'
#' Row and column names of the returned matrix are set to the `GEOID` column
#' of `acs_sf_data`, so \eqn{L_s} aligns directly with the weight matrix
#' \eqn{W} produced by [build_acs_precision_weights()].
#'
#' @details
#' **Contiguity rule.** `queen = TRUE` (default) applies Queen contiguity
#' (shared edge *or* vertex), which is recommended for irregular Census tract
#' polygons to avoid spurious disconnections at near-touching boundaries.
#'
#' **Properties of \eqn{L_s}.** \eqn{L_s} is symmetric positive semi-definite.
#' The penalty \eqn{x^\top L_s x = \sum_{(i,j)\in\mathcal{E}}(x_i-x_j)^2}
#' discourages spatially non-smooth deconvolution solutions.
#'
#' **Island tracts.** Tracts with no contiguous neighbours receive a zero
#' row/column in \eqn{L_s} and trigger a warning. They are not penalised
#' spatially but remain in the deconvolution system.
#'
#' @param acs_sf_data An [`sf`][sf::sf] object of Census tract polygons.
#'   Typically the `$sf` component of an `ACSData` object returned by
#'   [fetch_acs_vre_data()]. Should contain a `GEOID` column; if absent,
#'   row names fall back to `"tract_1"`, `"tract_2"`, ....
#' @param queen Logical. Use Queen contiguity (shared edge or vertex)?
#'   Default `TRUE`. Set to `FALSE` for Rook contiguity (shared edge only).
#'
#' @return A sparse symmetric matrix of class
#'   [`dsCMatrix`][Matrix::dsCMatrix-class] with dimensions
#'   \eqn{n \times n} where \eqn{n} is the number of tracts. Row and column
#'   names match the `GEOID` values of `acs_sf_data`.
#'
#' @examples
#' \dontrun{
#' acs <- fetch_acs_vre_data(
#'   state     = "IL",
#'   county    = "Cook",
#'   year_end  = 2022,
#'   variables = c(med_inc = "B19013_001")
#' )
#'
#' Ls <- build_empirical_laplacian(acs$sf)
#' dim(Ls)   # n_tracts x n_tracts
#'
#' # Laplacian property: row sums are zero
#' all(abs(Matrix::rowSums(Ls)) < 1e-12)
#'
#' # Plug into ADMM solver alongside W
#' W <- build_acs_precision_weights(acs)
#' }
#'
#' @seealso [fetch_acs_vre_data()], [build_acs_precision_weights()],
#'   [build_spatial_laplacian()], [spdep::poly2nb()]
#' @export
build_empirical_laplacian <- function(acs_sf_data, queen = TRUE) {
  if (!inherits(acs_sf_data, "sf"))
    stop("`acs_sf_data` must be an sf object ",
         "(e.g., the `$sf` field of an ACSData object).")

  n  <- nrow(acs_sf_data)
  nb <- spdep::poly2nb(acs_sf_data, queen = queen)

  # Identify and warn about island tracts
  n_islands <- sum(vapply(nb, function(x) identical(x, 0L), logical(1L)))
  if (n_islands > 0L)
    warning(n_islands, " island tract(s) with no contiguous neighbours ",
            "detected. These receive zero rows/columns in L_s and are ",
            "not spatially penalised.", call. = FALSE)

  # Binary adjacency matrix (sparse)
  A  <- .nb_to_sparse_adj(nb, n)

  # Degree matrix and Laplacian: L_s = D_s - A
  deg <- Matrix::rowSums(A)
  Ds  <- Matrix::Diagonal(n = n, x = deg)
  Ls  <- Matrix::forceSymmetric(Ds - A)

  # Dimnames: use GEOID column to align with build_acs_precision_weights()
  tract_nms <- if ("GEOID" %in% names(acs_sf_data))
    as.character(acs_sf_data[["GEOID"]])
  else
    paste0("tract_", seq_len(n))
  dimnames(Ls) <- list(tract_nms, tract_nms)

  Ls
}


# ================================================================
# 4.  Diagonal precision weight matrix W
# ================================================================

#' Build the diagonal precision (inverse-variance) weight matrix
#'
#' @description
#' Converts a matrix of standard errors \eqn{\hat{\sigma}_{is}} for the
#' observed data \eqn{Y} into the diagonal weight matrix
#' \deqn{W = \operatorname{diag}\!\bigl(1/\hat{\sigma}_{is}^2\bigr),}
#' used in the Weighted Least Squares (WLS) fitting criterion
#' \deqn{\min_X \; \|\, W^{1/2}(\operatorname{vec}(Y) - \tilde{M}
#'       \operatorname{vec}(X))\,\|_2^2 + \text{penalty}.}
#'
#' @details
#' The input matrix `Y_se` has the same shape as \eqn{Y}: \eqn{n} rows
#' (cells) and \eqn{S} columns (ACS periods).  It is column-major vectorised
#' to form the diagonal of \eqn{W}:
#' \deqn{W_{(i-1)S+s,\,(i-1)S+s} = 1 / \hat{\sigma}_{is}^2.}
#'
#' **Edge cases:**
#' \describe{
#'   \item{`NA` standard errors}{Weight is set to 0 (observation excluded
#'     from the fit).}
#'   \item{Zero standard errors}{Treatment controlled by `zero_se_action`:
#'     replace with the maximum finite weight (`"max_weight"`, default),
#'     raise a warning (`"warn"`), or throw an error (`"error"`).}
#' }
#'
#' @param Y_se Numeric matrix of standard errors with the same dimensions as
#'   \eqn{Y} (\eqn{n_{\text{cells}} \times n_{\text{periods}}}).  Must be
#'   non-negative.
#' @param zero_se_action Character scalar. Action when a zero SE is
#'   encountered: `"max_weight"` (default), `"warn"`, or `"error"`.
#'
#' @return A sparse diagonal matrix of class
#'   [`ddiMatrix`][Matrix::ddiMatrix-class] with dimension
#'   \eqn{(n \cdot S) \times (n \cdot S)}.  Names are `"cell_i:period_s"`
#'   when dimnames are present on `Y_se`.
#'
#' @examples
#' # Homoscedastic case: all SEs equal -> W is a scaled identity
#' Y_se_hom <- matrix(0.5, nrow = 100, ncol = 16)
#' W_hom <- build_precision_weights(Y_se_hom)
#' dim(W_hom)          # 1600 x 1600
#' range(Matrix::diag(W_hom))   # all 4 (= 1/0.25)
#'
#' # Heteroscedastic case
#' set.seed(1)
#' Y_se_het <- matrix(runif(100 * 16, 0.3, 1.2), nrow = 100)
#' W_het <- build_precision_weights(Y_se_het)
#'
#' @export
build_precision_weights <- function(Y_se,
                                    zero_se_action = c("max_weight",
                                                       "warn",
                                                       "error")) {
  zero_se_action <- match.arg(zero_se_action)

  if (!is.matrix(Y_se)) Y_se <- as.matrix(Y_se)
  if (any(Y_se < 0, na.rm = TRUE))
    stop("`Y_se` must contain only non-negative values.")

  # Transpose to (n_periods  x  n_cells) so that as.vector() produces the
  # survey-within-cell order that matches the ADMM y-vector:
  #   y = as.vector(t(Y))  =>  y[(i-1)*S+s] = Y[i,s]
  # With Y_se_t (S x n), as.vector(1/var_t)[(i-1)*S+s] = 1/sigma^2[i,s]  OK
  Y_se_t  <- t(Y_se)          # S \u00D7 n
  var_t   <- Y_se_t^2

  zero_idx <- which(var_t == 0 & !is.na(var_t))
  na_idx   <- which(is.na(var_t))

  # Handle zero SEs
  if (length(zero_idx) > 0L) {
    if (zero_se_action == "error") {
      stop(length(zero_idx), " zero SE value(s) found; ",
           "cannot compute inverse variance. ",
           "Set zero_se_action = 'max_weight' to cap these.")
    }
    if (zero_se_action == "warn") {
      warning(length(zero_idx),
              " zero SE value(s) replaced with the maximum finite weight.",
              call. = FALSE)
    }
    # Temporarily mark as NA so 1/var is computed only for valid entries
    var_t[zero_idx] <- NA_real_
  }

  # Compute weights in survey-within-cell order; NA SEs -> weight 0
  w_vec <- as.vector(1 / var_t)
  w_vec[is.na(w_vec)] <- 0

  # Back-fill zero-SE entries with max finite weight
  if (length(zero_idx) > 0L && zero_se_action %in% c("max_weight", "warn")) {
    max_w            <- max(w_vec, na.rm = TRUE)
    w_vec[zero_idx]  <- max_w
  }

  n_total <- length(w_vec)

  # Build sparse diagonal matrix
  W <- Matrix::Diagonal(n = n_total, x = w_vec)

  # Attach names: "cell_i:period_s" in survey-within-cell order
  rn <- rownames(Y_se)   # cell names (length n)
  cn <- colnames(Y_se)   # period names (length S)
  if (!is.null(rn) && !is.null(cn)) {
    # as.vector(t(outer(rn, cn))) -> survey-within-cell:
    #   (cell1:s1, cell1:s2, ..., cell1:sS, cell2:s1, ...)
    nms <- as.vector(t(outer(rn, cn, paste, sep = ":")))
    dimnames(W) <- list(nms, nms)
  }

  W
}


# ================================================================
# 5.  ACS convolution matrix M
# ================================================================

#' Build the ACS 5-year Convolution Matrix
#'
#' @description
#' Constructs the sparse S x T matrix M that encodes the ACS 5-year averaging
#' process: each row corresponds to one observed ACS window and has entries
#' of 1/5 in the five columns corresponding to the latent annual values it
#' averages.
#'
#' @details
#' For a sequence of ACS 5-year windows with end years
#' \code{start_year}, \code{start_year + 1}, ..., \code{end_year}:
#' \itemize{
#'   \item \eqn{S = \text{end\_year} - \text{start\_year} + 1} observed windows
#'   \item \eqn{T = S + 4} latent annual time points
#'   \item The first latent year is \code{start_year - 4}
#' }
#' Window \eqn{s} (1-indexed) covers latent columns \eqn{s} through \eqn{s+4},
#' each with weight 1/5:
#' \deqn{M_{s,t} = \frac{1}{5}, \quad t \in \{s, s+1, s+2, s+3, s+4\}}
#'
#' @param start_year Integer. End year of the first ACS 5-year window
#'   (e.g., 2013 for ACS 2009-2013). Must be >= 2013.
#' @param end_year Integer. End year of the last ACS 5-year window.
#'   Must be >= \code{start_year}.
#'
#' @return A sparse matrix of class [`dgCMatrix`][Matrix::dgCMatrix-class]
#'   with dimensions S  x  T where S = end_year - start_year + 1 and T = S + 4.
#'   Row names are `"window_<year>"` for each end year;
#'   column names are `"year_<year>"` for each latent year.
#'
#' @note
#'   The caller is responsible for ensuring `T > k` where `k` is the
#'   temporal difference order passed to [build_temporal_penalty()].
#'
#' @examples
#' # Single window
#' M1 <- build_acs_convolution_matrix(2013, 2013)
#' dim(M1)          # 1 x 5
#' as.matrix(M1)    # all 0.2
#'
#' # Three consecutive windows
#' M3 <- build_acs_convolution_matrix(2013, 2015)
#' dim(M3)                   # 3 x 7
#' Matrix::rowSums(M3)       # all 1
#'
#' @seealso [build_temporal_penalty()], [prep_acs_for_deconv()]
#' @export
build_acs_convolution_matrix <- function(start_year, end_year) {
  start_year <- as.integer(start_year)
  end_year   <- as.integer(end_year)

  if (end_year < start_year)
    stop("`end_year` must be >= `start_year`.")
  if (start_year < 2013L)
    stop("`start_year` must be >= 2013 (earliest ACS 5-year available).")

  S <- end_year - start_year + 1L
  T <- S + 4L

  i_idx <- rep(seq_len(S), each = 5L)
  j_idx <- unlist(lapply(seq_len(S), function(s) s + 0L:4L))

  M <- Matrix::sparseMatrix(
    i    = i_idx,
    j    = j_idx,
    x    = 1 / 5,
    dims = c(S, T)
  )

  rownames(M) <- paste0("window_", start_year:end_year)
  colnames(M) <- paste0("year_",   (start_year - 4L):end_year)

  M
}


# ================================================================
# 6.  Convenience inspector
# ================================================================

#' Summarise penalty matrix properties
#'
#' @description
#' Prints key properties of any of the three penalty matrices returned by
#' this module: dimensions, sparsity, and (for the Laplacian) whether row
#' sums are zero.
#'
#' @param mat A sparse matrix produced by [build_temporal_penalty()],
#'   [build_spatial_laplacian()], or [build_precision_weights()].
#' @param label Optional character label printed in the header.
#'
#' @return `mat` invisibly.
#' @export
summarise_penalty <- function(mat, label = deparse(substitute(mat))) {
  d       <- dim(mat)
  nnz     <- Matrix::nnzero(mat)
  density <- round(100 * nnz / prod(d), 3)
  cat(sprintf("\u2500\u2500 %s \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n", label))
  cat(sprintf("  Class      : %s\n", class(mat)))
  cat(sprintf("  Dimensions : %d \u00d7 %d\n", d[1], d[2]))
  cat(sprintf("  Non-zeros  : %d  (%.3f%% dense)\n", nnz, density))
  if (d[1] == d[2]) {
    rs <- Matrix::rowSums(mat)
    cat(sprintf("  Row-sum range : [%.6g, %.6g]\n", min(rs), max(rs)))
    cat(sprintf("  Symmetric  : %s\n",
                isSymmetric(mat, tol = 1e-12)))
  }
  cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  invisible(mat)
}
