## tests/testthat/test-simulation_suite.R

# ---- shared fixture -----------------------------------------------
make_suite_fixture <- function(seed = 5L) {
  gen <- SyntheticDataGenerator$new(
    n_row = 4, n_col = 4, T_years = 12,
    noise_sd = 0.3, seed = seed
  )
  gen$simulate()
  D  <- build_temporal_penalty(12, k = 2)
  Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
  W  <- build_precision_weights(matrix(0.3, nrow = 16, ncol = 8))
  list(gen = gen, D = D, Ls = Ls, W = W)
}

# ---- metric helpers -----------------------------------------------

# Extract internal helpers from the closure environment (no package install needed)
.env_ss  <- environment(run_simulation_suite)
.rmse_fn <- get(".rmse", envir = .env_ss)
.ccc_fn  <- get(".ccc",  envir = .env_ss)

test_that(".rmse is zero for perfect recovery", {
  X <- matrix(rnorm(20), 4, 5)
  expect_equal(.rmse_fn(X, X), 0)
})

test_that(".ccc is 1 for perfect agreement", {
  X <- matrix(rnorm(20), 4, 5)
  expect_equal(.ccc_fn(X, X), 1, tolerance = 1e-12)
})

test_that(".ccc is < 1 under location shift (Pearson stays 1)", {
  x <- 1:20
  y <- x + 5
  expect_lt(.ccc_fn(x, y), 1)
  expect_gt(.ccc_fn(x, y), 0)
})

test_that(".ccc is symmetric", {
  x <- 1:10
  expect_equal(.ccc_fn(x, -x + mean(x) * 2),
               .ccc_fn(-x + mean(x) * 2, x),
               tolerance = 1e-12)
})

# ---- run_simulation_suite: structure ------------------------------

test_that("run_simulation_suite returns correct structure", {
  f <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = c(0.05, 0.2),
    lambda_s_grid = c(0.05, 0.2)
  )
  expect_s3_class(suite, "SimulationSuite")
  expect_true(all(c("results", "best", "detail",
                    "timing_sec", "params") %in% names(suite)))
})

test_that("results data.frame has correct rows and columns", {
  f <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = c(0.05, 0.1, 0.5),
    lambda_s_grid = c(0.05, 0.2)
  )
  expect_equal(nrow(suite$results), 3 * 2)   # 6 combos
  expect_true(all(c("lambda_t", "lambda_s", "rmse",
                    "ccc", "converged", "n_iter") %in% names(suite$results)))
})

test_that("all RMSE values are finite and non-negative", {
  f <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = c(0.1, 0.5),
    lambda_s_grid = c(0.1)
  )
  expect_true(all(is.finite(suite$results$rmse)))
  expect_true(all(suite$results$rmse >= 0))
})

test_that("CCC values are in [-1, 1]", {
  f <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = c(0.1, 0.5),
    lambda_s_grid = c(0.1)
  )
  expect_true(all(suite$results$ccc >= -1 - 1e-10))
  expect_true(all(suite$results$ccc <=  1 + 1e-10))
})

# ---- best recovery -----------------------------------------------

test_that("best field is consistent with results", {
  f <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = c(0.05, 0.2, 1.0),
    lambda_s_grid = c(0.05, 0.2)
  )
  expect_equal(suite$best$rmse,
               min(suite$results$rmse, na.rm = TRUE),
               tolerance = 1e-10)
  expect_equal(dim(suite$best$X_hat),
               c(f$gen$params$n_cells, f$gen$params$T_years))
})

test_that("higher lambda_t produces temporally smoother X_hat (roughness test)", {
  # Temporal roughness = mean second-difference sum per cell.
  # Strong temporal regularisation must reduce it, regardless of RMSE vs X_true
  # (RMSE depends on signal properties; roughness is a direct penalty outcome).
  f <- make_suite_fixture()
  roughness <- function(X)
    mean(apply(X, 1, function(xi) sum(diff(diff(xi))^2)))

  X_low  <- suppressWarnings(
    solve_deconv(f$gen$Y, f$gen$M, f$W, f$D, f$Ls,
                 lambda_t = 0.01, lambda_s = 0.05)
  )
  X_high <- solve_deconv(f$gen$Y, f$gen$M, f$W, f$D, f$Ls,
                         lambda_t = 20.0, lambda_s = 0.05)

  expect_lt(roughness(X_high), roughness(X_low))
})

# ---- detail metrics -----------------------------------------------

test_that("detail$cell has n_cells rows with rmse/ccc columns", {
  f     <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = 0.1, lambda_s_grid = 0.1
  )
  expect_equal(nrow(suite$detail$cell), f$gen$params$n_cells)
  expect_true(all(c("rmse_cell", "ccc_cell") %in% names(suite$detail$cell)))
  expect_true(all(suite$detail$cell$rmse_cell >= 0))
})

test_that("detail$year has T_years rows", {
  f     <- make_suite_fixture()
  suite <- run_simulation_suite(
    gen = f$gen, D = f$D, Ls = f$Ls, W = f$W,
    lambda_t_grid = 0.1, lambda_s_grid = 0.1
  )
  expect_equal(nrow(suite$detail$year), f$gen$params$T_years)
})

# ---- Cholesky reuse efficiency -----------------------------------

test_that("suite with 6 combos is faster than 6 separate solve_deconv calls", {
  f <- make_suite_fixture()
  lt <- c(0.05, 0.1, 0.5)
  ls <- c(0.05, 0.2)

  t_suite <- system.time(
    run_simulation_suite(f$gen, f$D, f$Ls, W = f$W,
                         lambda_t_grid = lt, lambda_s_grid = ls)
  )[["elapsed"]]

  t_naive <- system.time(
    suppressWarnings(
      for (lti in lt) for (lsi in ls)
        solve_deconv(f$gen$Y, f$gen$M, f$W, f$D, f$Ls,
                     lambda_t = lti, lambda_s = lsi)
    )
  )[["elapsed"]]

  # Suite should be faster; allow 20 % slack for timing noise
  expect_lt(t_suite, t_naive * 1.2)
})

# ---- plot_cell_signals: smoke test --------------------------------

test_that("plot_cell_signals returns a ggplot without error", {
  skip_if_not_installed("ggplot2")
  f     <- make_suite_fixture()
  suite <- run_simulation_suite(
    f$gen, f$D, f$Ls, W = f$W,
    lambda_t_grid = 0.1, lambda_s_grid = 0.1
  )
  p <- plot_cell_signals(f$gen, suite$best$X_hat, cell_id = 3L)
  expect_s3_class(p, "ggplot")
})

test_that("plot_cell_signals uses a random cell when cell_id is NULL", {
  skip_if_not_installed("ggplot2")
  f     <- make_suite_fixture(seed = 77L)
  suite <- run_simulation_suite(
    f$gen, f$D, f$Ls, W = f$W,
    lambda_t_grid = 0.1, lambda_s_grid = 0.1
  )
  set.seed(42)
  p <- plot_cell_signals(f$gen, suite$best$X_hat)
  expect_s3_class(p, "ggplot")
})

# ---- plot_lambda_heatmap: smoke test ------------------------------

test_that("plot_lambda_heatmap returns a ggplot", {
  skip_if_not_installed("ggplot2")
  f     <- make_suite_fixture()
  suite <- run_simulation_suite(
    f$gen, f$D, f$Ls, W = f$W,
    lambda_t_grid = c(0.05, 0.2),
    lambda_s_grid = c(0.05, 0.2)
  )
  p <- plot_lambda_heatmap(suite, metric = "rmse")
  expect_s3_class(p, "ggplot")
})
