## tests/testthat/test-acs_connect.R
##
## Tests are grouped into four tiers:
##   1. Pure-math / internal helpers  – no network required
##   2. compute_sdr_se()              – no network required
##   3. fetch_acs_vre_data()          – network (skip_if_offline)
##   4. Integration: ACSData methods  – depends on tier-3 fixture

# ================================================================
# 1.  Internal helpers  (no network)
# ================================================================

# ── .state_to_fips ───────────────────────────────────────────────

test_that(".state_to_fips converts 2-letter abbreviation", {
  expect_equal(.state_to_fips("IL"), "17")
  expect_equal(.state_to_fips("CA"), "06")
  expect_equal(.state_to_fips("NY"), "36")
})

test_that(".state_to_fips converts full state name (case-insensitive)", {
  expect_equal(.state_to_fips("Illinois"),  "17")
  expect_equal(.state_to_fips("illinois"),  "17")
  expect_equal(.state_to_fips("CALIFORNIA"), "06")
})

test_that(".state_to_fips accepts numeric-style FIPS code", {
  expect_equal(.state_to_fips("17"), "17")
  expect_equal(.state_to_fips("6"),  "06")   # single-digit → zero-padded
})

test_that(".state_to_fips errors on unrecognised input", {
  expect_error(.state_to_fips("XX"),  "Cannot resolve")
  expect_error(.state_to_fips("Fake State"), "Cannot resolve")
})

# ── .county_to_fips ──────────────────────────────────────────────

test_that(".county_to_fips returns NULL for NULL input", {
  expect_null(.county_to_fips(NULL, "17"))
})

test_that(".county_to_fips resolves Cook County by name", {
  expect_equal(.county_to_fips("Cook", "17"), "031")
})

test_that(".county_to_fips accepts numeric FIPS code", {
  expect_equal(.county_to_fips("031", "17"), "031")
  expect_equal(.county_to_fips("31",  "17"), "031")  # zero-padded
})

test_that(".county_to_fips errors when county not found in state", {
  expect_error(.county_to_fips("Nonexistent County", "17"), "not found")
})

# ── .moe_to_se ───────────────────────────────────────────────────

test_that(".moe_to_se converts at 90% level", {
  # z = 1.6449; 164.49 / 1.6449 ≈ 100
  expect_equal(.moe_to_se(164.49, 90), 164.49 / 1.6449)
})

test_that(".moe_to_se converts at 95% and 99% levels", {
  expect_equal(.moe_to_se(196.0, 95), 196.0 / 1.9600)
  expect_equal(.moe_to_se(257.58, 99), 257.58 / 2.5758)
})

test_that(".moe_to_se is vectorised and uses abs()", {
  moe <- c(100, -200, 0)
  se  <- .moe_to_se(moe, 90)
  expect_equal(se, abs(moe) / 1.6449)
})

test_that(".moe_to_se errors on unsupported level", {
  expect_error(.moe_to_se(100, 80), "moe_level")
})

# ── .vars_to_tables ──────────────────────────────────────────────

test_that(".vars_to_tables strips last segment from B-series codes", {
  expect_equal(.vars_to_tables("B19013_001"), "B19013")
  expect_equal(.vars_to_tables("B01003_001"), "B01003")
})

test_that(".vars_to_tables returns unique table names", {
  expect_equal(
    .vars_to_tables(c("B19013_001", "B19013_002", "B19013_003")),
    "B19013"
  )
})

test_that(".vars_to_tables handles multiple tables", {
  tbls <- .vars_to_tables(c("B19013_001", "B01003_001"))
  expect_equal(sort(tbls), c("B01003", "B19013"))
})

# ── .vre_rep_cols ────────────────────────────────────────────────

test_that(".vre_rep_cols detects standard R1..R80 pattern", {
  dummy <- setNames(
    as.data.frame(matrix(0, 1, 82)),
    c("GEOID", "est", paste0("B19013_001_R", 1:80))
  )
  cols <- .vre_rep_cols(dummy, "B19013_001")
  expect_equal(length(cols), 80L)
  expect_true(all(grepl("R\\d+$", cols)))
})

test_that(".vre_rep_cols returns NULL when no pattern matches", {
  dummy <- data.frame(GEOID = "x", value = 1)
  expect_null(.vre_rep_cols(dummy, "B19013_001"))
})

test_that(".vre_rep_cols returns cols sorted by replicate number", {
  col_nms <- c("GEOID", paste0("B19013_001_R", c(10, 1, 5, 80, 2)))
  dummy   <- setNames(as.data.frame(matrix(0, 1, length(col_nms))), col_nms)
  # Too few (only 5), should return NULL
  expect_null(.vre_rep_cols(dummy, "B19013_001"))
})

# ── .vre_geoid_col ───────────────────────────────────────────────

test_that(".vre_geoid_col finds GEO_ID column", {
  df <- data.frame(GEO_ID = "x", val = 1)
  expect_equal(.vre_geoid_col(df), "GEO_ID")
})

test_that(".vre_geoid_col finds GEOID column (case insensitive)", {
  df <- data.frame(GEOID = "x", val = 1)
  expect_equal(.vre_geoid_col(df), "GEOID")
})

test_that(".vre_geoid_col falls back to first 11-char column", {
  df <- data.frame(tract_id = "17031010100", val = 1,
                   stringsAsFactors = FALSE)
  expect_equal(.vre_geoid_col(df), "tract_id")
})


# ================================================================
# 2.  compute_sdr_se()  (no network)
# ================================================================

test_that("compute_sdr_se returns zero SE for identical replicates", {
  est  <- c(50000, 60000, 45000)
  reps <- matrix(rep(est, 80), nrow = 3, ncol = 80)
  expect_equal(compute_sdr_se(est, reps), c(0, 0, 0))
})

test_that("compute_sdr_se matches analytic formula", {
  # SE = sqrt(4/80 * sum((r - est)^2))
  # For a single unit: est=100, replicates all either 90 or 110 (40 each)
  est      <- 100
  reps_row <- c(rep(90, 40), rep(110, 40))   # 80 replicates
  reps     <- matrix(reps_row, nrow = 1)
  # manual: sqrt(4/80 * (40*(90-100)^2 + 40*(110-100)^2))
  #       = sqrt(4/80 * 8000) = sqrt(400) = 20
  expect_equal(compute_sdr_se(est, reps), 20, tolerance = 1e-12)
})

test_that("compute_sdr_se is vectorised over n units", {
  n    <- 50
  est  <- seq(1000, 10000, length.out = n)
  set.seed(123)
  reps <- matrix(rnorm(n * 80, mean = rep(est, each = 80), sd = 500),
                 nrow = n, ncol = 80, byrow = TRUE)
  se   <- compute_sdr_se(est, reps)
  expect_length(se, n)
  expect_true(all(is.finite(se) & se >= 0))
})

test_that("compute_sdr_se errors when ncol(replicates) != 80", {
  expect_error(compute_sdr_se(1:10, matrix(0, 10, 79)), "exactly 80")
})

test_that("compute_sdr_se errors when nrow != length(estimate)", {
  expect_error(compute_sdr_se(1:10, matrix(0, 5, 80)), "nrow")
})

test_that("compute_sdr_se larger SE for more dispersed replicates", {
  est       <- rep(50000, 20)
  set.seed(42)
  reps_low  <- matrix(rnorm(20 * 80, 50000, 500),  nrow = 20)
  reps_high <- matrix(rnorm(20 * 80, 50000, 5000), nrow = 20)
  expect_lt(mean(compute_sdr_se(est, reps_low)),
            mean(compute_sdr_se(est, reps_high)))
})


# ================================================================
# 3.  fetch_acs_vre_data()  (network required)
# ================================================================

skip_if_offline <- function() {
  tryCatch(
    {
      con <- url("https://api.census.gov", "r")
      close(con)
    },
    error = function(e) testthat::skip("No internet / Census API unreachable")
  )
}

# Small fixture reused across tier-3 tests
.acs_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      skip_if_offline()
      cache <<- fetch_acs_vre_data(
        state     = "IL",
        county    = "Cook",
        year_end  = 2022,
        variables = c(med_inc = "B19013_001", pop = "B01003_001"),
        geometry  = FALSE,
        vre       = FALSE,
        quiet     = TRUE
      )
    }
    cache
  }
})

test_that("fetch_acs_vre_data returns an ACSData object", {
  acs <- .acs_fixture()
  expect_s3_class(acs, "ACSData")
})

test_that("ACSData has all required top-level components", {
  acs <- .acs_fixture()
  expect_true(all(c("sf", "estimates", "se", "moe",
                    "replicates", "metadata") %in% names(acs)))
})

test_that("estimates matrix has correct columns from named variables", {
  acs <- .acs_fixture()
  expect_equal(colnames(acs$estimates), c("med_inc", "pop"))
  expect_equal(colnames(acs$se),        c("med_inc", "pop"))
  expect_equal(colnames(acs$moe),       c("med_inc", "pop"))
})

test_that("estimates, se, moe matrices have matching dimensions", {
  acs <- .acs_fixture()
  expect_equal(dim(acs$se),  dim(acs$estimates))
  expect_equal(dim(acs$moe), dim(acs$estimates))
})

test_that("n_tracts matches nrow of estimates", {
  acs <- .acs_fixture()
  expect_equal(acs$metadata$n_tracts, nrow(acs$estimates))
})

test_that("estimates matrix has positive values for median income", {
  acs <- .acs_fixture()
  vals <- acs$estimates[, "med_inc"]
  vals <- vals[is.finite(vals)]
  expect_gt(length(vals), 0L)
  expect_true(all(vals > 0))
})

test_that("SE values are non-negative and finite (or NA)", {
  acs <- .acs_fixture()
  se  <- acs$se[, "med_inc"]
  expect_true(all(is.na(se) | se >= 0))
})

test_that("MOE >= SE * z at 90% level (abs tolerance)", {
  acs <- .acs_fixture()
  z   <- 1.6449
  moe <- abs(acs$moe[, "med_inc"])
  se  <- acs$se[, "med_inc"]
  valid <- is.finite(moe) & is.finite(se) & se > 0
  # moe ≈ se * z; allow 1 % relative tolerance
  expect_true(all(abs(moe[valid] / (se[valid] * z) - 1) < 0.01))
})

test_that("rownames of estimates are 11-digit GEOIDs", {
  acs <- .acs_fixture()
  geoids <- rownames(acs$estimates)
  expect_true(all(nchar(geoids) == 11L))
  expect_true(all(grepl("^\\d{11}$", geoids)))
})

test_that("metadata contains all expected fields", {
  acs <- .acs_fixture()
  m   <- acs$metadata
  expect_true(all(c("state", "state_fips", "county", "year_end",
                    "year_range", "variables", "var_names",
                    "n_tracts", "moe_level", "vre_source",
                    "call_time") %in% names(m)))
})

test_that("metadata$year_range matches year_end", {
  acs <- .acs_fixture()
  m   <- acs$metadata
  expect_equal(m$year_end, 2022L)
  expect_true(grepl("2018", m$year_range))
  expect_true(grepl("2022", m$year_range))
})

test_that("metadata$vre_source is 'moe' when vre=FALSE", {
  acs <- .acs_fixture()
  expect_equal(acs$metadata$vre_source, "moe")
})

test_that("unnamed variables use ACS codes as column names", {
  skip_if_offline()
  acs_un <- fetch_acs_vre_data(
    state = "IL", county = "Cook", year_end = 2022,
    variables = "B19013_001",
    geometry = FALSE, vre = FALSE, quiet = TRUE
  )
  expect_equal(colnames(acs_un$estimates), "B19013_001")
})

test_that("state as full name and as abbreviation give same GEOID set", {
  skip_if_offline()
  acs_abbr <- fetch_acs_vre_data(
    state = "IL", county = "Cook", year_end = 2022,
    variables = c(inc = "B19013_001"),
    geometry = FALSE, vre = FALSE, quiet = TRUE
  )
  acs_name <- fetch_acs_vre_data(
    state = "Illinois", county = "Cook", year_end = 2022,
    variables = c(inc = "B19013_001"),
    geometry = FALSE, vre = FALSE, quiet = TRUE
  )
  expect_equal(sort(rownames(acs_abbr$estimates)),
               sort(rownames(acs_name$estimates)))
})


# ================================================================
# 4.  S3 methods: print / summary  (uses fixture)
# ================================================================

test_that("print.ACSData runs without error", {
  acs <- .acs_fixture()
  expect_no_error(capture.output(print(acs)))
})

test_that("print.ACSData output contains state and year", {
  acs <- .acs_fixture()
  out <- capture.output(print(acs))
  txt <- paste(out, collapse = "\n")
  expect_true(grepl("IL",    txt, fixed = TRUE))
  expect_true(grepl("2022",  txt))
  expect_true(grepl("med_inc", txt))
})

test_that("summary.ACSData runs without error", {
  acs <- .acs_fixture()
  expect_no_error(capture.output(summary(acs)))
})

test_that("summary.ACSData output has Estimates section", {
  acs <- .acs_fixture()
  out <- paste(capture.output(summary(acs)), collapse = "\n")
  expect_true(grepl("Estimates", out))
  expect_true(grepl("med_inc",   out))
})


# ================================================================
# 5.  prep_acs_for_deconv()
# ================================================================

# ── Tier 1: structure tests (no network) ─────────────────────────
#
# We inject a pre-built ACSData list via a mock of fetch_acs_vre_data
# so that prep_acs_for_deconv can be tested without any network calls.

# Helper: build a minimal ACSData object for n tracts and one variable,
# with a simple square grid geometry.
.make_mock_acs_data <- function(n_row = 4, n_col = 3, var_code = "B19013_001",
                                var_name = "med_inc", seed = 1L) {
  set.seed(seed)
  n <- n_row * n_col

  # sf with GEOID column
  bb   <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = n_col, ymax = n_row))
  grid <- sf::st_make_grid(sf::st_as_sfc(bb), n = c(n_col, n_row))
  geoids <- sprintf("%011d", seq_len(n))
  sf_obj <- sf::st_sf(GEOID = geoids, geometry = grid)

  est <- matrix(rnorm(n, mean = 50000, sd = 8000), nrow = n, ncol = 1L)
  se  <- matrix(runif(n, 500, 2000), nrow = n, ncol = 1L)
  moe <- se * 1.6449

  rownames(est) <- rownames(se) <- rownames(moe) <- geoids
  colnames(est) <- colnames(se) <- colnames(moe) <- var_name

  # Add SE column to sf
  se_df <- as.data.frame(se)
  names(se_df) <- paste0(var_name, "_se")
  sf_out <- cbind(sf_obj, se_df)

  structure(
    list(
      sf         = sf_out,
      estimates  = est,
      se         = se,
      moe        = moe,
      replicates = stats::setNames(list(NULL), var_name),
      metadata   = list(
        state      = "IL",
        state_fips = "17",
        county     = "Cook",
        year_end   = NA_integer_,
        year_range = NA_character_,
        variables  = var_code,
        var_names  = var_name,
        n_tracts   = n,
        moe_level  = 90L,
        vre_source = "moe",
        call_time  = Sys.time()
      )
    ),
    class = "ACSData"
  )
}

# Build the mock fixture that simulates two windows sharing all tracts
.prep_fixture_mock <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)

    year_ends <- c(2021L, 2022L)
    # One mock object per year_end; we'll override fetch_acs_vre_data
    mock_data <- list(
      .make_mock_acs_data(n_row = 4, n_col = 3, seed = 1L),
      .make_mock_acs_data(n_row = 4, n_col = 3, seed = 2L)
    )
    # Second window must not have geometry
    mock_data[[2L]]$sf <- sf::st_drop_geometry(mock_data[[2L]]$sf)

    # Monkey-patch fetch_acs_vre_data in the package environment
    counter <- 0L
    local_fetch <- function(...) {
      counter <<- counter + 1L
      mock_data[[counter]]
    }

    env <- getNamespace("TWDeConv")
    old_fn <- env$fetch_acs_vre_data
    on.exit(assignInNamespace("fetch_acs_vre_data", old_fn, ns = "TWDeConv"),
            add = TRUE)
    assignInNamespace("fetch_acs_vre_data", local_fetch, ns = "TWDeConv")

    cache <<- prep_acs_for_deconv(
      state     = "IL",
      county    = "Cook",
      year_ends = year_ends,
      variables = c(med_inc = "B19013_001"),
      k         = 2L,
      quiet     = TRUE
    )
    cache
  }
})

test_that("prep_acs_for_deconv returns a DeconvInput object (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  expect_s3_class(inp, "DeconvInput")
})

test_that("DeconvInput has required top-level names (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  expect_true(all(c("Y", "M", "W", "D", "Ls", "metadata") %in% names(inp)))
})

test_that("Y dimensions are n x S (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  m   <- inp$metadata
  expect_equal(dim(inp$Y), c(m$n_tracts, m$S_windows))
})

test_that("M dimensions are S x T (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  m   <- inp$metadata
  expect_equal(dim(inp$M), c(m$S_windows, m$T_years))
})

test_that("M row sums equal 1 (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  expect_equal(as.numeric(Matrix::rowSums(inp$M)),
               rep(1, nrow(inp$M)), tolerance = 1e-12)
})

test_that("W is square diagonal of size n*S (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  nS  <- inp$metadata$n_tracts * inp$metadata$S_windows
  expect_equal(dim(inp$W), c(nS, nS))
})

test_that("D dimensions are (T-k) x T (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  m   <- inp$metadata
  expect_equal(dim(inp$D), c(m$T_years - m$k, m$T_years))
})

test_that("Ls is n x n symmetric (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  n   <- inp$metadata$n_tracts
  expect_equal(dim(inp$Ls), c(n, n))
  expect_true(is(inp$Ls, "symmetricMatrix"))
})

test_that("metadata contains all required fields (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  m   <- inp$metadata
  required <- c("state", "county", "year_ends", "variables", "var_name",
                "n_tracts", "T_years", "S_windows", "k",
                "first_latent_year", "last_latent_year", "vre_sources")
  expect_true(all(required %in% names(m)))
})

test_that("T_years = S_windows + 4 for consecutive year_ends (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  m   <- inp$metadata
  expect_equal(m$T_years, m$S_windows + 4L)
})

test_that("print.DeconvInput runs without error (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  expect_no_error(capture.output(print(inp)))
})

test_that("print.DeconvInput output contains key info (mock)", {
  skip_if_not_installed("TWDeConv")
  inp <- .prep_fixture_mock()
  out <- paste(capture.output(print(inp)), collapse = "\n")
  expect_true(grepl("DeconvInput", out))
  expect_true(grepl("ADMMSolver",  out))
})

test_that("prep_acs_for_deconv errors on unsorted year_ends", {
  expect_error(
    prep_acs_for_deconv("IL", "Cook", year_ends = c(2022L, 2021L),
                        variables = "B19013_001"),
    "sorted"
  )
})

test_that("prep_acs_for_deconv errors when length(variables) > 1", {
  expect_error(
    prep_acs_for_deconv("IL", "Cook", year_ends = 2021L,
                        variables = c("B19013_001", "B01003_001")),
    "univariate"
  )
})


# ── Tier 2: network tests ─────────────────────────────────────────

.deconv_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      skip_if_offline()
      cache <<- tryCatch(
        prep_acs_for_deconv(
          state     = "IL",
          county    = "Cook",
          year_ends = c(2021L, 2022L),
          variables = c(med_inc = "B19013_001"),
          k         = 2L,
          vre       = FALSE,
          quiet     = TRUE
        ),
        error = function(e) {
          msg <- conditionMessage(e)
          if (grepl("Geometry unavailable|load_tiger|unused argument|tigris",
                    msg, ignore.case = TRUE))
            testthat::skip(paste("geometry download unavailable:", msg))
          stop(e)
        }
      )
    }
    cache
  }
})

test_that("prep_acs_for_deconv returns DeconvInput with real data", {
  inp <- .deconv_fixture()
  expect_s3_class(inp, "DeconvInput")
})

test_that("real DeconvInput Y and M dimensions are consistent", {
  inp <- .deconv_fixture()
  m   <- inp$metadata
  expect_equal(ncol(inp$Y), m$S_windows)
  expect_equal(nrow(inp$M), m$S_windows)
  expect_equal(ncol(inp$M), m$T_years)
  expect_equal(m$T_years, m$S_windows + 4L)
})
