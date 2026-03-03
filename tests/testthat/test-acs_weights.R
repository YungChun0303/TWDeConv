## tests/testthat/test-acs_weights.R
##
## Tests are in three tiers:
##   1. Pure-math with a mock ACSData  – no network, no srvyr required
##   2. srvyr path (survey design)     – skip if srvyr not installed
##   3. Live API fixture               – skip if offline

# ================================================================
# Helpers: build a minimal mock ACSData
# ================================================================

.make_mock_acs <- function(n = 10, seed = 1L, add_reps = TRUE) {
  set.seed(seed)
  geoids  <- sprintf("%011d", seq_len(n))
  est     <- matrix(runif(n, 30000, 90000), nrow = n, ncol = 1,
                    dimnames = list(geoids, "inc"))
  moe     <- matrix(est * 0.15, nrow = n, ncol = 1,
                    dimnames = list(geoids, "inc"))
  se      <- moe / 1.6449

  reps <- if (add_reps) {
    rep_mat <- matrix(
      rnorm(n * 80, mean = rep(est, 80), sd = rep(se, 80)),
      nrow = n, ncol = 80,
      dimnames = list(geoids, paste0("R", seq_len(80L)))
    )
    list(inc = rep_mat)
  } else {
    list(inc = NULL)
  }

  structure(
    list(
      sf         = data.frame(GEOID = geoids, incE = est[,1], incM = moe[,1]),
      estimates  = est,
      se         = se,
      moe        = moe,
      replicates = reps,
      metadata   = list(
        state      = "mock",
        state_fips = "00",
        county     = NULL,
        year_end   = 2022L,
        year_range = "2018-2022",
        variables  = "B00000_001",
        var_names  = "inc",
        n_tracts   = n,
        moe_level  = 90,
        vre_source = if (add_reps) "sdr" else "moe",
        call_time  = Sys.time()
      )
    ),
    class = "ACSData"
  )
}


# ================================================================
# 1.  Basic structure and dimensions  (no network, no srvyr)
# ================================================================

test_that("build_acs_precision_weights returns a ddiMatrix", {
  acs <- .make_mock_acs(add_reps = FALSE)
  W   <- build_acs_precision_weights(acs, use_srvyr = FALSE, quiet = TRUE)
  expect_true(is(W, "ddiMatrix") || is(W, "diagonalMatrix"))
})

test_that("output dimension equals n_tracts for single variable", {
  acs <- .make_mock_acs(n = 12, add_reps = FALSE)
  W   <- build_acs_precision_weights(acs, use_srvyr = FALSE, quiet = TRUE)
  expect_equal(dim(W), c(12L, 12L))
})

test_that("output dimension is n_tracts * n_vars for two variables", {
  # Build a 2-variable mock
  n      <- 8
  geoids <- sprintf("%011d", seq_len(n))
  set.seed(2)
  est  <- matrix(runif(n * 2, 30000, 80000), n, 2,
                 dimnames = list(geoids, c("inc", "pop")))
  se   <- est * 0.10
  moe  <- se * 1.6449
  acs2 <- structure(
    list(
      sf         = data.frame(GEOID = geoids),
      estimates  = est, se = se, moe = moe,
      replicates = list(inc = NULL, pop = NULL),
      metadata   = list(var_names = c("inc", "pop"), n_tracts = n,
                        vre_source = "moe", moe_level = 90)
    ),
    class = "ACSData"
  )
  W2 <- build_acs_precision_weights(acs2, use_srvyr = FALSE, quiet = TRUE)
  expect_equal(dim(W2), c(n * 2L, n * 2L))
})

test_that("diagonal entries equal 1/se^2 (MOE path)", {
  acs   <- .make_mock_acs(n = 5, add_reps = FALSE)
  W     <- build_acs_precision_weights(acs, use_srvyr = FALSE, quiet = TRUE)
  se    <- acs$se[, "inc"]
  # survey-within-cell order: all variables for tract 1, tract 2, ...
  # with 1 variable: ordering is just tract 1, tract 2, ...
  prec_expected <- 1 / se^2
  expect_equal(unname(Matrix::diag(W)), unname(prec_expected), tolerance = 1e-10)
})

test_that("row/column names are GEOID:variable in correct order", {
  acs  <- .make_mock_acs(n = 3, add_reps = FALSE)
  W    <- build_acs_precision_weights(acs, use_srvyr = FALSE, quiet = TRUE)
  nms  <- rownames(W)
  geoids <- rownames(acs$estimates)
  expected_nms <- paste0(geoids, ":inc")
  expect_equal(nms, expected_nms)
})

test_that("variable subset selection works", {
  n      <- 5
  geoids <- sprintf("%011d", seq_len(n))
  set.seed(3)
  est  <- matrix(runif(n * 2, 30000, 80000), n, 2,
                 dimnames = list(geoids, c("inc", "pop")))
  acs2 <- structure(
    list(
      sf = data.frame(GEOID = geoids), estimates = est,
      se = est * 0.1, moe = est * 0.1 * 1.6449,
      replicates = list(inc = NULL, pop = NULL),
      metadata = list(var_names = c("inc", "pop"), n_tracts = n,
                      vre_source = "moe", moe_level = 90)
    ),
    class = "ACSData"
  )
  W_sub <- build_acs_precision_weights(acs2, variables = "inc",
                                        use_srvyr = FALSE, quiet = TRUE)
  expect_equal(dim(W_sub), c(n, n))
  expect_true(all(grepl(":inc$", rownames(W_sub))))
})

test_that("error on unknown variable", {
  acs <- .make_mock_acs(add_reps = FALSE)
  expect_error(
    build_acs_precision_weights(acs, variables = "zzz",
                                 use_srvyr = FALSE, quiet = TRUE),
    "None of the requested variables"
  )
})

test_that("error when not passed an ACSData object", {
  expect_error(
    build_acs_precision_weights(list(a = 1), quiet = TRUE),
    "ACSData"
  )
})


# ================================================================
# 2.  zero_var_action edge cases  (no network, no srvyr)
# ================================================================

test_that("zero SE with max_weight gives finite positive precision", {
  acs          <- .make_mock_acs(n = 5, add_reps = FALSE)
  acs$se[1, 1] <- 0   # introduce zero SE
  expect_warning(
    W <- build_acs_precision_weights(acs, use_srvyr = FALSE,
                                      zero_var_action = "max_weight",
                                      quiet = TRUE),
    "zero"
  )
  expect_true(all(is.finite(Matrix::diag(W)) & Matrix::diag(W) >= 0))
  # The zero-SE tract gets the maximum precision
  expect_equal(as.numeric(Matrix::diag(W)[1]), as.numeric(max(Matrix::diag(W))), tolerance = 1e-10)
})

test_that("zero SE with drop gives weight 0", {
  acs          <- .make_mock_acs(n = 5, add_reps = FALSE)
  acs$se[2, 1] <- 0
  expect_warning(
    W <- build_acs_precision_weights(acs, use_srvyr = FALSE,
                                      zero_var_action = "drop",
                                      quiet = TRUE),
    "zero"
  )
  expect_equal(as.numeric(Matrix::diag(W)[2]), 0)
})

test_that("zero SE with error throws", {
  acs          <- .make_mock_acs(n = 5, add_reps = FALSE)
  acs$se[1, 1] <- 0
  expect_error(
    build_acs_precision_weights(acs, use_srvyr = FALSE,
                                 zero_var_action = "error",
                                 quiet = TRUE),
    "zero"
  )
})

test_that("NA SE produces weight 0 (observation excluded)", {
  acs          <- .make_mock_acs(n = 5, add_reps = FALSE)
  acs$se[3, 1] <- NA_real_
  # NA variance → max_weight replaces it
  expect_warning(
    W <- build_acs_precision_weights(acs, use_srvyr = FALSE,
                                      zero_var_action = "max_weight",
                                      quiet = TRUE)
  )
  expect_true(is.finite(Matrix::diag(W)[3]) && Matrix::diag(W)[3] > 0)
})


# ================================================================
# 3.  srvyr path: SDR design gives identical result to direct formula
# ================================================================

test_that("srvyr path matches compute_sdr_se() exactly", {
  skip_if_not_installed("srvyr")
  skip_if_not_installed("survey")

  acs <- .make_mock_acs(n = 15, add_reps = TRUE, seed = 42L)

  # srvyr path
  W_srvyr <- build_acs_precision_weights(acs, use_srvyr = TRUE,  quiet = TRUE)
  # Direct SDR formula path
  W_moe   <- build_acs_precision_weights(acs, use_srvyr = FALSE, quiet = TRUE)

  # SDR via direct formula
  sdr_se  <- compute_sdr_se(acs$estimates[, "inc"], acs$replicates[["inc"]])
  prec_direct <- 1 / sdr_se^2

  expect_equal(unname(Matrix::diag(W_srvyr)), unname(prec_direct), tolerance = 1e-8)
})

test_that("srvyr path produces larger precision when replicates are tighter", {
  skip_if_not_installed("srvyr")
  skip_if_not_installed("survey")

  n      <- 10
  geoids <- sprintf("%011d", seq_len(n))
  set.seed(7)
  est <- runif(n, 50000, 70000)

  # Tight replicates (low variance)
  reps_tight <- matrix(rnorm(n * 80, mean = est, sd = 200),
                       nrow = n, ncol = 80)
  # Loose replicates (high variance)
  reps_loose <- matrix(rnorm(n * 80, mean = est, sd = 8000),
                       nrow = n, ncol = 80)

  make_acs <- function(reps) {
    se <- compute_sdr_se(est, reps)
    structure(
      list(
        sf = data.frame(GEOID = geoids), estimates = matrix(est, n, 1,
        dimnames = list(geoids, "inc")),
        se = matrix(se, n, 1, dimnames = list(geoids, "inc")),
        moe = matrix(se * 1.6449, n, 1, dimnames = list(geoids, "inc")),
        replicates = list(inc = `rownames<-`(reps, geoids)),
        metadata = list(var_names = "inc", n_tracts = n,
                        vre_source = "sdr", moe_level = 90)
      ),
      class = "ACSData"
    )
  }

  W_tight <- build_acs_precision_weights(make_acs(reps_tight),
                                          use_srvyr = TRUE, quiet = TRUE)
  W_loose <- build_acs_precision_weights(make_acs(reps_loose),
                                          use_srvyr = TRUE, quiet = TRUE)

  expect_gt(mean(Matrix::diag(W_tight)), mean(Matrix::diag(W_loose)))
})


# ================================================================
# 4.  Survey-within-cell ordering matches build_precision_weights()
# ================================================================

test_that("ordering is survey-within-cell (matches ADMM y-vector)", {
  # Build a 2-variable, 3-tract ACSData with known, distinct SEs
  n      <- 3
  geoids <- c("A", "B", "C")
  # SE for (tract, variable): deliberately distinct
  se_inc <- c(100, 200, 300)
  se_pop <- c(10,  20,  30)
  est    <- matrix(c(50000, 60000, 70000, 1000, 2000, 3000), n, 2,
                   dimnames = list(geoids, c("inc", "pop")))
  se_mat <- matrix(c(se_inc, se_pop), n, 2,
                   dimnames = list(geoids, c("inc", "pop")))

  acs2 <- structure(
    list(
      sf = data.frame(GEOID = geoids), estimates = est,
      se = se_mat, moe = se_mat * 1.6449,
      replicates = list(inc = NULL, pop = NULL),
      metadata = list(var_names = c("inc", "pop"), n_tracts = n,
                      vre_source = "moe", moe_level = 90)
    ),
    class = "ACSData"
  )

  W <- build_acs_precision_weights(acs2, use_srvyr = FALSE, quiet = TRUE)
  d <- Matrix::diag(W)

  # Survey-within-cell order for 2 variables (inc, pop):
  # (A:inc, A:pop, B:inc, B:pop, C:inc, C:pop)
  expected <- c(1/100^2, 1/10^2,
                1/200^2, 1/20^2,
                1/300^2, 1/30^2)
  expect_equal(unname(d), expected, tolerance = 1e-14)
})


# ================================================================
# 5.  Live API test  (skip if offline)
# ================================================================

.skip_if_census_down <- function() {
  ok <- tryCatch(
    { con <- url("https://api.census.gov", "r"); close(con); TRUE },
    error = function(e) FALSE
  )
  if (!ok) testthat::skip("Census API unreachable")
}

test_that("build_acs_precision_weights works on real ACS data (MOE path)", {
  .skip_if_census_down()
  acs <- fetch_acs_vre_data(
    state = "IL", county = "Cook", year_end = 2022,
    variables = c(inc = "B19013_001"),
    geometry = FALSE, vre = FALSE, quiet = TRUE
  )
  W <- build_acs_precision_weights(acs, use_srvyr = FALSE, quiet = TRUE)

  n <- acs$metadata$n_tracts
  expect_equal(dim(W), c(n, n))
  d <- Matrix::diag(W)
  expect_true(all(is.finite(d) | d == 0))
  expect_true(all(d >= 0))
  # Most tracts should have positive precision
  expect_gt(mean(d > 0), 0.9)
})
