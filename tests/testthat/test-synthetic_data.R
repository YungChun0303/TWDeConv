## tests/testthat/test-synthetic_data.R

test_that("SyntheticDataGenerator produces correct matrix dimensions", {
  gen <- SyntheticDataGenerator$new(
    n_row = 5, n_col = 4, T_years = 12, window = 5, seed = 1L
  )
  gen$simulate()

  n       <- 5 * 4          # 20 cells
  T       <- 12
  w       <- 5
  n_p     <- T - w + 1L     # 8 periods

  expect_equal(dim(gen$X), c(n,   T  ))
  expect_equal(dim(gen$M), c(n_p, T  ))
  expect_equal(dim(gen$Y), c(n,   n_p))
})

test_that("Convolution matrix rows sum to 1", {
  gen <- SyntheticDataGenerator$new(T_years = 10, window = 5)
  gen$build_convolution_matrix()
  expect_equal(unname(rowSums(gen$M)), rep(1, nrow(gen$M)), tolerance = 1e-12)
})

test_that("Y equals X %*% t(M) plus bounded noise", {
  gen <- SyntheticDataGenerator$new(
    n_row = 3, n_col = 3, T_years = 8, window = 5,
    noise_sd = 0, seed = 99L             # zero noise → exact equality
  )
  gen$simulate()

  Y_expected <- gen$X %*% t(gen$M)
  expect_equal(gen$Y, Y_expected, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("as_sf returns an sf object with the right number of rows", {
  gen <- SyntheticDataGenerator$new(n_row = 4, n_col = 4, T_years = 8,
                                    seed = 7L)
  gen$simulate()
  out <- gen$as_sf()
  expect_s3_class(out, "sf")
  expect_equal(nrow(out), 4 * 4)
})

test_that("simulate() is reproducible with same seed", {
  g1 <- SyntheticDataGenerator$new(seed = 42L)$simulate()
  g2 <- SyntheticDataGenerator$new(seed = 42L)$simulate()
  expect_equal(g1$X, g2$X)
  expect_equal(g1$Y, g2$Y)
})
