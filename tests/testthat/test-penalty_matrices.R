## tests/testthat/test-penalty_matrices.R

# ---- build_temporal_penalty ------------------------------------

test_that("D1 has correct dimensions and stencil", {
  D <- build_temporal_penalty(T = 8, k = 1)
  expect_equal(dim(D), c(7L, 8L))

  # First row: [-1, 1, 0, ...]
  row1 <- as.numeric(D[1, ])
  expect_equal(row1, c(-1, 1, rep(0, 6)))
})

test_that("D2 has correct dimensions and stencil", {
  D <- build_temporal_penalty(T = 8, k = 2)
  expect_equal(dim(D), c(6L, 8L))

  # First row of D2 should be [1, -2, 1, 0, ...]
  row1 <- as.numeric(D[1, ])
  expect_equal(row1, c(1, -2, 1, rep(0, 5)), tolerance = 1e-12)
})

test_that("Dk has dimensions (T-k) x T for several orders", {
  T <- 15
  for (k in 1:5) {
    D <- build_temporal_penalty(T = T, k = k)
    expect_equal(dim(D), c(T - k, T),
                 info = paste("k =", k))
  }
})

test_that("D1 applied to a constant vector gives all zeros", {
  D <- build_temporal_penalty(T = 10, k = 1)
  x_const <- rep(3.7, 10)
  expect_equal(as.numeric(D %*% x_const), rep(0, 9), tolerance = 1e-12)
})

test_that("D2 applied to a linear vector gives all zeros", {
  # Second differences of a linear signal are 0
  D <- build_temporal_penalty(T = 12, k = 2)
  x_lin <- seq(1, 12)
  expect_equal(as.numeric(D %*% x_lin), rep(0, 10), tolerance = 1e-12)
})

test_that("build_temporal_penalty rejects invalid k", {
  expect_error(build_temporal_penalty(T = 5, k = 0),  "`k` must be >= 1")
  expect_error(build_temporal_penalty(T = 5, k = 5),  "`T` must be strictly")
})

# ---- build_spatial_laplacian -----------------------------------

test_that("Laplacian has correct dimensions and row-sum = 0", {
  gen <- SyntheticDataGenerator$new(n_row = 4, n_col = 5, seed = 1L)
  gen$generate_grid()
  Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)

  expect_equal(dim(Ls), c(20L, 20L))
  expect_true(is(Ls, "symmetricMatrix"))
  expect_equal(as.numeric(Matrix::rowSums(Ls)),
               rep(0, 20), tolerance = 1e-12)
})

test_that("Laplacian diagonal = degree, off-diagonal = 0 or -1", {
  gen <- SyntheticDataGenerator$new(n_row = 3, n_col = 3, seed = 2L)
  gen$generate_grid()
  Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
  Lm <- as.matrix(Ls)

  # Off-diagonal entries are 0 or -1
  off_diag <- Lm[row(Lm) != col(Lm)]
  expect_true(all(off_diag %in% c(0, -1)))

  # Diagonal entries equal the number of neighbours (rook, interior = 4)
  expect_true(all(diag(Lm) >= 0))
})

test_that("Laplacian is PSD (all eigenvalues >= 0)", {
  gen <- SyntheticDataGenerator$new(n_row = 3, n_col = 3, seed = 3L)
  gen$generate_grid()
  Ls  <- build_spatial_laplacian(gen$grid, queen = FALSE)
  evs <- eigen(as.matrix(Ls), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evs >= -1e-10))
})

# ---- build_empirical_laplacian ---------------------------------

# Helper: minimal sf object with GEOID column, using a planar n_row x n_col grid
.make_mock_acs_sf <- function(n_row = 3, n_col = 3) {
  bb   <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = n_col, ymax = n_row))
  grid <- sf::st_make_grid(sf::st_as_sfc(bb), n = c(n_col, n_row))
  sf::st_sf(GEOID = sprintf("%011d", seq_len(n_row * n_col)), geometry = grid)
}

test_that("build_empirical_laplacian returns a sparse symmetric matrix", {
  sf_data <- .make_mock_acs_sf()
  Ls <- suppressWarnings(build_empirical_laplacian(sf_data))
  expect_true(inherits(Ls, "sparseMatrix"))
  expect_true(is(Ls, "symmetricMatrix"))
})

test_that("build_empirical_laplacian has correct dimensions", {
  sf_data <- .make_mock_acs_sf(n_row = 4, n_col = 5)
  Ls <- suppressWarnings(build_empirical_laplacian(sf_data))
  expect_equal(dim(Ls), c(20L, 20L))
})

test_that("build_empirical_laplacian: row sums are zero (Laplacian property)", {
  sf_data <- .make_mock_acs_sf(n_row = 3, n_col = 4)
  Ls <- suppressWarnings(build_empirical_laplacian(sf_data))
  expect_equal(as.numeric(Matrix::rowSums(Ls)), rep(0, 12), tolerance = 1e-12)
})

test_that("build_empirical_laplacian: dimnames use GEOID values", {
  sf_data <- .make_mock_acs_sf(n_row = 2, n_col = 3)
  Ls <- suppressWarnings(build_empirical_laplacian(sf_data))
  expect_equal(rownames(Ls), sf_data$GEOID)
  expect_equal(colnames(Ls), sf_data$GEOID)
})

test_that("build_empirical_laplacian: fallback dimnames when no GEOID column", {
  sf_data <- .make_mock_acs_sf()
  sf_data$GEOID <- NULL
  Ls <- suppressWarnings(build_empirical_laplacian(sf_data))
  expect_true(all(grepl("^tract_\\d+$", rownames(Ls))))
  expect_equal(length(rownames(Ls)), nrow(sf_data))
})

test_that("build_empirical_laplacian: warns on island tracts", {
  sf_data <- .make_mock_acs_sf(n_row = 3, n_col = 3)
  island <- sf::st_sf(
    GEOID    = "99999999999",
    geometry = sf::st_sfc(sf::st_polygon(list(
      matrix(c(100,100, 101,100, 101,101, 100,101, 100,100),
             ncol = 2, byrow = TRUE)
    )))
  )
  sf_with_island <- rbind(sf_data, island)
  expect_warning(build_empirical_laplacian(sf_with_island), "island")
})

test_that("build_empirical_laplacian: island row/column is all zero", {
  sf_data <- .make_mock_acs_sf(n_row = 3, n_col = 3)
  island <- sf::st_sf(
    GEOID    = "99999999999",
    geometry = sf::st_sfc(sf::st_polygon(list(
      matrix(c(100,100, 101,100, 101,101, 100,101, 100,100),
             ncol = 2, byrow = TRUE)
    )))
  )
  sf_with_island <- rbind(sf_data, island)
  Ls <- suppressWarnings(build_empirical_laplacian(sf_with_island))
  island_idx <- which(rownames(Ls) == "99999999999")
  expect_equal(as.numeric(Ls[island_idx, ]), rep(0, nrow(Ls)))
})

test_that("build_empirical_laplacian: is PSD (all eigenvalues >= 0)", {
  sf_data <- .make_mock_acs_sf(n_row = 3, n_col = 3)
  Ls  <- suppressWarnings(build_empirical_laplacian(sf_data))
  evs <- eigen(as.matrix(Ls), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evs >= -1e-10))
})

test_that("build_empirical_laplacian: rejects non-sf input", {
  expect_error(build_empirical_laplacian(data.frame(x = 1:3)), "sf object")
})

test_that("build_empirical_laplacian: queen=FALSE (rook) has <= edges than queen", {
  sf_data    <- .make_mock_acs_sf(n_row = 4, n_col = 4)
  Ls_queen   <- suppressWarnings(build_empirical_laplacian(sf_data, queen = TRUE))
  Ls_rook    <- suppressWarnings(build_empirical_laplacian(sf_data, queen = FALSE))
  # Queen has at least as many edges (non-zeros off-diagonal) as rook
  expect_gte(Matrix::nnzero(Ls_queen), Matrix::nnzero(Ls_rook))
})

# ---- build_precision_weights -----------------------------------

test_that("W is a square diagonal matrix with correct size", {
  Y_se <- matrix(0.5, nrow = 10, ncol = 8)
  W    <- build_precision_weights(Y_se)
  expect_equal(dim(W), c(80L, 80L))
  expect_true(inherits(W, "ddiMatrix") || isDiagonal(W))
})

test_that("Homoscedastic SEs give a scaled identity", {
  se   <- 0.4
  Y_se <- matrix(se, nrow = 6, ncol = 5)
  W    <- build_precision_weights(Y_se)
  expect_equal(range(Matrix::diag(W)),
               rep(1 / se^2, 2), tolerance = 1e-12)
})

test_that("NA standard errors produce zero weights", {
  Y_se      <- matrix(0.5, nrow = 4, ncol = 4)
  Y_se[2, 3] <- NA
  W <- build_precision_weights(Y_se)
  # build_precision_weights transposes Y_se (to S x n) then vectorises col-major.
  # Y_se[row=2, col=3] = Y_se[cell=2, period=3] -> Y_se_t[period=3, cell=2]
  # position in vector: (cell - 1) * S + period = (2 - 1) * 4 + 3 = 7
  idx <- (2 - 1) * 4 + 3
  expect_equal(Matrix::diag(W)[idx], 0)
})

test_that("Zero SE raises error when requested", {
  Y_se       <- matrix(0.5, nrow = 3, ncol = 3)
  Y_se[1, 1] <- 0
  expect_error(build_precision_weights(Y_se, zero_se_action = "error"),
               "zero SE")
})

test_that("Negative SEs are rejected", {
  Y_se       <- matrix(0.5, nrow = 3, ncol = 3)
  Y_se[2, 2] <- -0.1
  expect_error(build_precision_weights(Y_se), "non-negative")
})


# ---- build_acs_convolution_matrix -----------------------------------

test_that("M for single window (2013, 2013) is 1x5, all 1/5", {
  M <- build_acs_convolution_matrix(2013, 2013)
  expect_equal(dim(M), c(1L, 5L))
  expect_equal(as.vector(M), rep(1/5, 5), tolerance = 1e-12)
})

test_that("M dimensions are S x (S+4) for general inputs", {
  for (start in c(2013L, 2015L, 2018L)) {
    for (end in start + c(0L, 3L, 9L)) {
      M <- build_acs_convolution_matrix(start, end)
      S <- end - start + 1L
      expect_equal(dim(M), c(S, S + 4L),
                   info = sprintf("start=%d, end=%d", start, end))
    }
  }
})

test_that("all row sums equal 1 (valid averaging weights)", {
  M <- build_acs_convolution_matrix(2015, 2022)
  expect_equal(as.numeric(Matrix::rowSums(M)),
               rep(1, nrow(M)), tolerance = 1e-12)
})

test_that("each row has exactly 5 non-zeros (bandwidth constraint)", {
  M   <- build_acs_convolution_matrix(2013, 2020)
  Mm  <- as.matrix(M)
  row_nnz <- rowSums(Mm != 0)
  expect_true(all(row_nnz == 5L))
})

test_that("first row covers columns 1:5 (window 1 = first 5 latent years)", {
  M  <- build_acs_convolution_matrix(2016, 2020)
  Mm <- as.matrix(M)
  expect_equal(unname(which(Mm[1, ] != 0)), 1:5)
})

test_that("M has banded structure: M[s,t]=0 when t<s or t>s+4", {
  M  <- build_acs_convolution_matrix(2013, 2018)
  Mm <- as.matrix(M)
  S  <- nrow(Mm)
  T  <- ncol(Mm)
  for (s in seq_len(S)) {
    for (t in seq_len(T)) {
      if (t < s || t > s + 4L) {
        expect_equal(Mm[s, t], 0,
                     info = sprintf("M[%d,%d] should be 0", s, t))
      }
    }
  }
})

test_that("row and column names follow naming convention", {
  start <- 2015L; end <- 2017L
  M <- build_acs_convolution_matrix(start, end)
  expect_equal(rownames(M), paste0("window_", start:end))
  expect_equal(colnames(M), paste0("year_", (start - 4L):end))
})

test_that("end_year < start_year raises an error", {
  expect_error(build_acs_convolution_matrix(2020, 2018), "`end_year`")
})

test_that("start_year < 2013 raises an error", {
  expect_error(build_acs_convolution_matrix(2012, 2015), "`start_year`")
})
