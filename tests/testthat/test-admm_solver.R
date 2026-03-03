## tests/testthat/test-admm_solver.R

# ---- shared fixture ------------------------------------------------
make_problem <- function(n_row = 4, n_col = 4, T = 12, k = 2,
                         noise_sd = 0.3, seed = 7L) {
  gen <- SyntheticDataGenerator$new(
    n_row = n_row, n_col = n_col, T_years = T,
    noise_sd = noise_sd, seed = seed
  )
  gen$simulate()
  D  <- build_temporal_penalty(T, k = k)
  Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)
  W  <- build_precision_weights(
    matrix(noise_sd, nrow = n_row * n_col, ncol = T - 5 + 1)
  )
  list(gen = gen, D = D, Ls = Ls, W = W)
}

# ---- ADMMSolver: construction ------------------------------------

test_that("ADMMSolver initialises without error", {
  p <- make_problem()
  expect_no_error(
    ADMMSolver$new(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                   lambda_t = 0.1, lambda_s = 0.05)
  )
})

test_that("ADMMSolver rejects negative lambda", {
  p <- make_problem()
  expect_error(
    ADMMSolver$new(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                   lambda_t = -0.1, lambda_s = 0.05),
    "non-negative"
  )
})

test_that("ADMMSolver rejects non-positive rho", {
  p <- make_problem()
  expect_error(
    ADMMSolver$new(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                   lambda_t = 0.1, lambda_s = 0.05, rho = 0),
    "strictly positive"
  )
})

# ---- setup: Kronecker dimensions ---------------------------------

test_that("setup() produces Kronecker matrices with correct dimensions", {
  p  <- make_problem(n_row = 3, n_col = 3, T = 10, k = 1)
  n  <- 9; T <- 10; S <- 6; k <- 1
  s  <- ADMMSolver$new(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                       lambda_t = 0.1, lambda_s = 0.05)
  s$setup()

  expect_equal(dim(s$M_kron),  c(n * S, n * T))
  expect_equal(dim(s$Dt_kron), c(n * (T - k), n * T))
  expect_equal(dim(s$Ds_kron), c(n * T, n * T))
  expect_equal(dim(s$A_sys),   c(n * T, n * T))
})

test_that("A_sys is symmetric positive definite", {
  p <- make_problem()
  s <- ADMMSolver$new(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05)
  s$setup()

  expect_true(Matrix::isSymmetric(s$A_sys, tol = 1e-10))

  # All eigenvalues positive (PD check on small system)
  evs <- eigen(as.matrix(s$A_sys), symmetric = TRUE,
               only.values = TRUE)$values
  expect_true(all(evs > 0))
})

# ---- run: output dimensions & type -------------------------------

test_that("run() produces X_hat with correct dimensions", {
  p <- make_problem()
  g <- p$gen
  s <- ADMMSolver$new(g$Y, g$M, p$W, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05)
  s$run()
  expect_equal(dim(s$X_hat), c(g$params$n_cells, g$params$T_years))
  expect_true(is.matrix(s$X_hat))
})

test_that("run() populates history with iteration counts", {
  p <- make_problem()
  s <- ADMMSolver$new(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05)
  suppressWarnings(s$run(max_iter = 200L))
  expect_true(is.data.frame(s$history))
  expect_true(all(c("iter", "r_primal", "r_dual") %in% names(s$history)))
  expect_true(nrow(s$history) >= 1L)
  expect_true(nrow(s$history) <= 200L)
})

# ---- solve_deconv wrapper ----------------------------------------

test_that("solve_deconv returns a matrix by default", {
  p <- make_problem()
  g <- p$gen
  out <- solve_deconv(g$Y, g$M, p$W, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(g$params$n_cells, g$params$T_years))
})

test_that("solve_deconv returns ADMMSolver when return_solver = TRUE", {
  p   <- make_problem()
  g   <- p$gen
  sol <- solve_deconv(g$Y, g$M, p$W, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05,
                      return_solver = TRUE)
  expect_true(inherits(sol, "ADMMSolver"))
  expect_true(is.matrix(sol$X_hat))
})

# ---- data fidelity -----------------------------------------------

test_that("With tiny lambda, Y ≈ X_hat %*% t(M) (data fit)", {
  # No noise, tiny regularisation → near-perfect reconstruction of Y
  gen <- SyntheticDataGenerator$new(n_row = 3, n_col = 3,
                                    T_years = 10, noise_sd = 0,
                                    seed = 99L)
  gen$simulate()
  D  <- build_temporal_penalty(10, k = 1)
  Ls <- build_spatial_laplacian(gen$grid, queen = FALSE)

  # rho scaled down to match small lambda; warnings about non-convergence
  # are suppressed — the data-fit test is valid regardless
  X_hat <- suppressWarnings(
    solve_deconv(gen$Y, gen$M, W = NULL, D = D, Ls = Ls,
                 lambda_t = 1e-4, lambda_s = 1e-4,
                 rho = 1e-3, max_iter = 3000L, tol_abs = 1e-5)
  )

  Y_hat <- X_hat %*% t(gen$M)
  rmse  <- sqrt(mean((Y_hat - gen$Y)^2))
  expect_lt(rmse, 0.05)
})

# ---- regularisation effect ---------------------------------------

test_that("Higher lambda_t produces smoother temporal signals", {
  p <- make_problem(seed = 11L)
  g <- p$gen

  # low-lambda case may hit max_iter; suppress the non-convergence warning
  # as the roughness comparison holds even at a partial solution
  X_low  <- suppressWarnings(
    solve_deconv(g$Y, g$M, p$W, p$D, p$Ls,
                 lambda_t = 0.01, lambda_s = 0.01)
  )
  X_high <- solve_deconv(g$Y, g$M, p$W, p$D, p$Ls,
                         lambda_t = 5.0, lambda_s = 0.01)

  # Second-order temporal roughness of each row
  roughness <- function(X) {
    mean(apply(X, 1, function(xi) sum(diff(diff(xi))^2)))
  }
  expect_lt(roughness(X_high), roughness(X_low))
})

# ---- compute_objective -------------------------------------------

test_that("compute_objective returns finite named vector", {
  p   <- make_problem()
  sol <- solve_deconv(p$gen$Y, p$gen$M, p$W, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05,
                      return_solver = TRUE)
  obj <- sol$compute_objective()
  expect_true(all(is.finite(obj)))
  expect_true(all(c("data", "temporal", "spatial", "total") %in% names(obj)))
  expect_equal(unname(obj["total"]),
               unname(obj["data"] + obj["temporal"] + obj["spatial"]),
               tolerance = 1e-10)
})

# ---- NULL W (unweighted LS) --------------------------------------

test_that("NULL weight matrix is accepted (unweighted LS)", {
  p   <- make_problem()
  g   <- p$gen
  out <- solve_deconv(g$Y, g$M, W = NULL, p$D, p$Ls,
                      lambda_t = 0.1, lambda_s = 0.05)
  expect_equal(dim(out), c(g$params$n_cells, g$params$T_years))
})
