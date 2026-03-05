## ============================================================
## R/acs_connect.R
## ACS 5-year data ingestion with Variance Replicate Estimates
##
## Data sources
##   Estimates + MOE : Census API via tidycensus::get_acs()
##   Geometries      : Census TIGER/Line via tidycensus (geometry=TRUE)
##   VRE replicates  : Census Bureau file server (www2.census.gov)
##                     Falls back to MOE-derived SE on download failure.
##
## VRE methodology
##   The ACS uses Successive Difference Replication (SDR) to quantify
##   sampling variance.  For each summary statistic the Census Bureau
##   publishes 80 replicate estimates theta^_r (r = 1...80) so that
##     SE(theta^) = sqrt(4/80 . Sigma_r (theta^_r - theta^)^2)
##   See: Census Bureau (2018) ACS Design & Methodology, Ch. 12.
## ============================================================

# ================================================================
# 1.  Internal helpers
# ================================================================

#' Convert state name / abbreviation / FIPS to zero-padded 2-digit FIPS
#' @noRd
.state_to_fips <- function(state) {
  fc     <- tidycensus::fips_codes
  state  <- trimws(as.character(state))

  # Already a 2-digit FIPS code
  if (grepl("^\\d{1,2}$", state))
    return(sprintf("%02d", as.integer(state)))

  fc_uniq <- fc[!duplicated(fc$state_code), ]

  # 2-letter abbreviation
  idx <- match(toupper(state), toupper(fc_uniq$state))
  if (!is.na(idx)) return(fc_uniq$state_code[idx])

  # Full state name
  idx <- match(toupper(state), toupper(fc_uniq$state_name))
  if (!is.na(idx)) return(fc_uniq$state_code[idx])

  stop("Cannot resolve state identifier '", state,
       "'. Supply a 2-letter abbreviation, full name, or FIPS code.")
}

#' Convert county name / partial name / FIPS to 3-digit FIPS
#' Returns NULL if county is NULL (= all counties).
#' @noRd
.county_to_fips <- function(county, state_fips) {
  if (is.null(county)) return(NULL)
  county <- trimws(as.character(county))
  if (grepl("^\\d{1,3}$", county)) return(sprintf("%03d", as.integer(county)))
  # Partial name match within state
  fc  <- tidycensus::fips_codes
  sub <- fc[fc$state_code == state_fips, ]
  idx <- grep(county, sub$county, ignore.case = TRUE)
  if (length(idx) == 0)
    stop("County '", county, "' not found in state FIPS ", state_fips, ".")
  if (length(idx) > 1)
    warning("Multiple county matches for '", county,
            "'; using first: ", sub$county[idx[1]], call. = FALSE)
  sub$county_code[idx[1]]
}

#' Extract unique base table names from ACS variable codes
#' "B19013_001" -> "B19013";  "S1701_C01_001" -> "S1701"
#' @noRd
.vars_to_tables <- function(variables) {
  unique(sub("_[^_]+$", "", variables))
}

#' MOE -> SE conversion (ACS publishes 90 % CI by default)
#' @noRd
.moe_to_se <- function(moe, level = 90) {
  z_tbl <- c("90" = 1.6449, "95" = 1.9600, "99" = 2.5758)
  z     <- unname(z_tbl[as.character(level)])  # [ ] returns NA on miss; unname strips label
  if (is.na(z)) stop("`moe_level` must be 90, 95, or 99.")
  abs(moe) / z
}

#' Build a prioritised list of Census Bureau VRE download URLs
#'
#' The Census Bureau's file server layout has changed slightly across years.
#' We try several URL patterns and accept the first that responds.
#' @noRd
.vre_url_candidates <- function(year, state_fips, table_name) {
  base <- "https://www2.census.gov/programs-surveys/acs/replicate_estimates"
  c(
    # Pattern 1 -- primary (used 2016 - present)
    sprintf("%s/%d/data/5-year/%s/ACS_%d_5YR_%s_%s.zip",
            base, year, state_fips, year, table_name, state_fips),
    # Pattern 2 -- alternate underscore separator
    sprintf("%s/%d/data/5_year/%s/ACS_%d_5YR_%s_%s.zip",
            base, year, state_fips, year, table_name, state_fips),
    # Pattern 3 -- no year in filename (older releases)
    sprintf("%s/%d/data/5-year/%s/%s_%s.zip",
            base, year, state_fips, table_name, state_fips)
  )
}

#' Try to download a VRE zip from the Census Bureau; return parsed data.frame or NULL
#' @noRd
.download_vre <- function(year, state_fips, table_name, cache_dir, quiet,
                           timeout_sec = 60) {
  candidates <- .vre_url_candidates(year, state_fips, table_name)
  zip_path   <- file.path(cache_dir,
                           sprintf("VRE_%d_%s_%s.zip",
                                   year, state_fips, table_name))

  # -- Try each URL candidate -----------------------------------------
  if (!file.exists(zip_path)) {
    downloaded <- FALSE
    for (url in candidates) {
      if (!quiet) message("  [VRE] Trying: ", basename(url))
      resp <- tryCatch(
        httr::GET(url,
                  httr::write_disk(zip_path, overwrite = TRUE),
                  httr::timeout(timeout_sec)),
        error = function(e) NULL
      )
      if (!is.null(resp) && !httr::http_error(resp)) {
        downloaded <- TRUE
        if (!quiet) message("  [VRE] Downloaded ", table_name, " for state ", state_fips)
        break
      }
      if (file.exists(zip_path)) unlink(zip_path)
    }
    if (!downloaded) return(NULL)
  } else {
    if (!quiet) message("  [VRE] Cache hit: ", table_name, " for state ", state_fips)
  }

  # -- Unzip and read CSV --------------------------------------------
  csv_files <- tryCatch(
    utils::unzip(zip_path, exdir = file.path(cache_dir,
                  sprintf("VRE_%d_%s_%s", year, state_fips, table_name))),
    error = function(e) NULL
  )
  if (is.null(csv_files) || length(csv_files) == 0) {
    if (!quiet) message("  [VRE] Unzip failed for ", table_name)
    unlink(zip_path)
    return(NULL)
  }

  csv_path <- csv_files[grepl("\\.csv$", csv_files, ignore.case = TRUE)][1L]
  if (is.na(csv_path)) return(NULL)

  df <- tryCatch(
    utils::read.csv(csv_path, stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA", "N")),
    error = function(e) NULL
  )
  df
}

#' Fetch geometry via tigris for a given geography
#' Returns an sf object or NULL on failure.
#' @noRd
.fetch_tigris_geometry <- function(geography, state_fips, county_fips, year) {
  if (!requireNamespace("tigris", quietly = TRUE)) return(NULL)
  switch(geography,
    "tract"       = tigris::tracts(
                      state = state_fips, county = county_fips, year = year),
    "block group" = tigris::block_groups(
                      state = state_fips, county = county_fips, year = year),
    "county"      = tigris::counties(state = state_fips, year = year),
    "zcta"        = tigris::zctas(year = year),
    stop("tigris geometry fallback not implemented for geography = '",
         geography, "'")
  )
}

#' Find the GEOID column in a VRE data.frame
#' Census VRE files use various names: GEO_ID, GEOID, geo_id, etc.
#' @noRd
.vre_geoid_col <- function(vre_df) {
  candidates <- grep("^(GEO_?ID|GEOID|geo_?id)$",
                     names(vre_df), ignore.case = TRUE, value = TRUE)
  if (length(candidates) == 0L) {
    # Fallback: first column whose values look like Census GEOIDs
    # (all-digit strings of consistent length, stripped of any "1400000US" prefix)
    for (nm in names(vre_df)) {
      vals <- sub("^\\d+US", "", as.character(vre_df[[nm]]))
      vals <- na.omit(vals)
      if (length(vals) > 0L && all(grepl("^\\d+$", vals)) &&
          length(unique(nchar(vals))) == 1L)
        return(nm)
    }
    return(names(vre_df)[1L])  # last resort: first column
  }
  candidates[1L]
}

#' Extract the 80 replicate estimate columns for one variable from a VRE data.frame
#'
#' VRE files use several naming conventions across Census years:
#'   B19013_001_R1  ... B19013_001_R80
#'   B19013_001_R001 ... B19013_001_R080
#'   B19013_001E_R1  ... B19013_001E_R80
#' Returns a character vector of column names (length 80) sorted by
#' replicate number, or NULL if no pattern matches.
#' @noRd
.vre_rep_cols <- function(vre_df, variable) {
  col_nms <- names(vre_df)
  patterns <- c(
    sprintf("^%s_R(\\d+)$",  variable),        # B19013_001_R1
    sprintf("^%sE_R(\\d+)$", variable),        # B19013_001E_R1
    sprintf("^%s_r(\\d+)$",  variable),        # lowercase
    sprintf("^%sE_r(\\d+)$", variable)
  )
  for (pat in patterns) {
    m <- grep(pat, col_nms, value = TRUE)
    if (length(m) >= 80L) {
      nums    <- as.integer(regmatches(m, regexpr("\\d+$", m)))
      ordered <- m[order(nums)]
      return(ordered[seq_len(80L)])
    }
  }
  NULL
}


# ================================================================
# 2.  compute_sdr_se()  (exported)
# ================================================================

#' Compute SDR Standard Errors from 80 ACS Replicate Estimates
#'
#' @description
#' Applies the Successive Difference Replication (SDR) variance formula
#' used by the U.S. Census Bureau for the American Community Survey:
#'
#' \deqn{\widehat{\text{SE}}(\hat{\theta}) =
#'   \sqrt{\frac{4}{80} \sum_{r=1}^{80} (\hat{\theta}_r - \hat{\theta})^2}}
#'
#' where \eqn{\hat{\theta}} is the point estimate and
#' \eqn{\hat{\theta}_1, \ldots, \hat{\theta}_{80}} are the 80 replicate
#' estimates from the ACS Variance Replicate Estimate (VRE) tables.
#'
#' @param estimate Numeric vector of length \eqn{n} (one value per geographic
#'   unit, e.g., Census tract).
#' @param replicates Numeric matrix of dimensions \eqn{n \times 80} whose
#'   row \eqn{i} contains the 80 replicate estimates for unit \eqn{i}.
#'
#' @return Numeric vector of length \eqn{n} of SDR standard errors.
#'
#' @references
#' U.S. Census Bureau (2018). \emph{Understanding and Using American
#' Community Survey Data: What All Data Users Need to Know.}
#' U.S. Government Publishing Office, Washington, DC.
#' Chapter 8, Appendix A.
#'
#' @examples
#' set.seed(42)
#' n    <- 50
#' est  <- rnorm(n, mean = 55000, sd = 8000)
#' reps <- matrix(rnorm(n * 80, mean = est, sd = 3000), nrow = n, ncol = 80)
#' se   <- compute_sdr_se(est, reps)
#' summary(se)
#'
#' @export
compute_sdr_se <- function(estimate, replicates) {
  estimate   <- as.numeric(estimate)
  replicates <- as.matrix(replicates)
  if (ncol(replicates) != 80L)
    stop("`replicates` must have exactly 80 columns.")
  if (nrow(replicates) != length(estimate))
    stop("`nrow(replicates)` must equal `length(estimate)`.")
  diffs <- sweep(replicates, 1L, estimate, FUN = "-")  # n \u00D7 80
  sqrt((4 / 80) * rowSums(diffs^2))
}


# ================================================================
# 3.  fetch_acs_vre_data()  (exported)
# ================================================================

#' Fetch ACS 5-Year Estimates, Geometries, and Variance Replicate Estimates
#'
#' @description
#' Retrieves ACS 5-year data for a user-specified geographic unit via
#' [tidycensus::get_acs()], augmented with Variance Replicate Estimate (VRE)
#' tables from the Census Bureau file server.  The result is a structured
#' list containing the point-estimate [`sf`][sf::sf] object, standard error
#' matrices, and (when available) the full \eqn{n \times 80} replicate
#' matrices required for Successive Difference Replication (SDR) variance
#' estimation.
#'
#' @details
#' ## Standard errors
#'
#' Two SE sources are supported, selected automatically:
#'
#' \describe{
#'   \item{`"moe"` (fallback)}{`SE = MOE / z`, where `z = 1.6449` for the
#'     default 90 % confidence intervals published by the ACS.  Always
#'     available.}
#'   \item{`"sdr"` (preferred)}{SE computed from the 80 VRE replicate
#'     estimates via [compute_sdr_se()].  Requires a successful download of
#'     the Census Bureau's VRE CSV files.}
#' }
#'
#' ## VRE file download
#'
#' VRE tables are pre-computed replicate-based estimates published by the
#' Census Bureau at
#' `https://www2.census.gov/programs-surveys/acs/replicate_estimates/`.
#' They are available for B- and C-series tables at the **tract** level from
#' 2013 to the most recent 5-year release.  For all other geographies (block
#' group, ZCTA, county, etc.) VRE is not published and the function
#' automatically falls back to MOE-derived SE.  The function tries several
#' URL patterns and caches downloaded zip files in `cache_dir` to avoid
#' redundant downloads.
#'
#' ## Geometry
#'
#' Geometries are fetched via the Census TIGER/Line server through
#' `tidycensus`.  Set `geometry = FALSE` if the TIGER server is unavailable
#' or the geometry is not needed.
#'
#' ## Geographic unit
#'
#' The `geography` argument is passed directly to [tidycensus::get_acs()].
#' Common values include `"tract"` (default), `"block group"`, `"county"`,
#' and `"zcta"`.  Note that `county` filtering is ignored by `tidycensus`
#' for `"zcta"` — a warning is issued in that case.
#'
#' @param state Character. State name (e.g., `"Illinois"`), 2-letter
#'   abbreviation (`"IL"`), or 2-digit FIPS code (`"17"`).
#' @param county Character or `NULL`. County name, partial name, or 3-digit
#'   FIPS code.  `NULL` (default) retrieves all units in the state.
#'   Ignored for `geography = "zcta"`.
#' @param year_end Integer. End year of the ACS 5-year period
#'   (e.g., `2022` for ACS 2018-2022).  Must be >= 2013.
#' @param variables Named or unnamed character vector of ACS variable codes,
#'   e.g. `c(med_inc = "B19013_001", pop = "B01003_001")`.  Names become
#'   column prefixes in the output; defaults to the variable codes.
#' @param geography Character. Geographic unit passed to
#'   [tidycensus::get_acs()].  Common values: `"tract"` (default),
#'   `"block group"`, `"county"`, `"zcta"`.  VRE-based SDR standard errors
#'   are only available for `"tract"`; all other geographies fall back to
#'   MOE-derived SE.
#' @param geometry Logical. Attach polygon geometries via TIGER/Line
#'   (default `TRUE`).  Set to `FALSE` if the geometry server is unavailable.
#' @param vre Logical. Attempt to download VRE replicate tables and compute
#'   SDR standard errors (default `TRUE`).  Automatically disabled for
#'   geographies other than `"tract"`.  Falls back to MOE on failure.
#' @param moe_level Numeric. ACS margin-of-error confidence level: `90`
#'   (default, ACS standard), `95`, or `99`.
#' @param cache_dir Character. Directory for caching downloaded VRE files
#'   (default `tempdir()`).  Use a persistent path across sessions to avoid
#'   re-downloading large state files.
#' @param timeout_sec Numeric. HTTP timeout in seconds for VRE downloads
#'   (default `90`).
#' @param census_api_key Character or `NULL`. Census API key.  If `NULL`,
#'   uses the key previously installed via [tidycensus::census_api_key()].
#' @param quiet Logical. Suppress progress messages (default `FALSE`).
#'
#' @return A list of class `"ACSData"` with components:
#' \describe{
#'   \item{`$sf`}{[`sf`][sf::sf] data frame with one row per geographic unit:
#'     `GEOID`, `NAME`, `{var}E` (estimate), `{var}M` (MOE),
#'     `{var}_se` (standard error), and (when `vre = TRUE` and download
#'     succeeds) `{var}_R1` ... `{var}_R80` replicate columns.
#'     `geometry` is `NA` when `geometry = FALSE`.}
#'   \item{`$estimates`}{n  x  p numeric matrix of point estimates.
#'     Rows named by `GEOID`; columns named by `var_names`.}
#'   \item{`$se`}{n  x  p numeric matrix of standard errors (SDR or MOE-derived).}
#'   \item{`$moe`}{n  x  p numeric matrix of raw margins of error at `moe_level` %.}
#'   \item{`$replicates`}{Named list of p elements, each an n  x  80 numeric
#'     matrix of replicate estimates.  `NULL` entries indicate variables for
#'     which VRE data were unavailable.}
#'   \item{`$metadata`}{Named list: `state`, `state_fips`, `county`,
#'     `geography`, `year_end`, `year_range`, `variables`, `var_names`,
#'     `n_units`, `moe_level`, `vre_source` (`"sdr"` or `"moe"`),
#'     `call_time`.}
#' }
#'
#' @examples
#' \dontrun{
#' # Tract-level (default) -- Cook County, IL, ACS 2018-2022
#' acs <- fetch_acs_vre_data(
#'   state     = "IL",
#'   county    = "Cook",
#'   year_end  = 2022,
#'   variables = c(med_inc = "B19013_001"),
#'   geography = "tract",
#'   cache_dir = "~/acs_cache"
#' )
#'
#' # Block-group level (VRE unavailable; falls back to MOE-derived SE)
#' acs_bg <- fetch_acs_vre_data(
#'   state     = "IL",
#'   county    = "Cook",
#'   year_end  = 2022,
#'   variables = c(med_inc = "B19013_001"),
#'   geography = "block group"
#' )
#'
#' # ZCTA level (county filter is ignored)
#' acs_zip <- fetch_acs_vre_data(
#'   state     = "IL",
#'   year_end  = 2022,
#'   variables = c(med_inc = "B19013_001"),
#'   geography = "zcta"
#' )
#' }
#'
#' @seealso [compute_sdr_se()], [build_precision_weights()],
#'   [build_spatial_laplacian()], [tidycensus::get_acs()]
#' @export
fetch_acs_vre_data <- function(state,
                                county        = NULL,
                                year_end,
                                variables,
                                geography     = "tract",
                                geometry      = TRUE,
                                vre           = TRUE,
                                moe_level     = 90,
                                cache_dir     = tempdir(),
                                timeout_sec   = 90,
                                census_api_key = NULL,
                                quiet         = FALSE) {

  # -- 0. Validate & resolve identifiers --------------------------
  if (!requireNamespace("tidycensus", quietly = TRUE))
    stop("Package 'tidycensus' is required.")

  geography <- tolower(trimws(geography))

  # VRE replicate tables are only published for Census tracts
  vre_geographies <- "tract"
  if (vre && !geography %in% vre_geographies) {
    if (!quiet)
      message(sprintf(
        "  VRE replicate tables are not available for geography = '%s'; ",
        geography), "falling back to MOE-derived SE.")
    vre <- FALSE
  }

  if (vre && !requireNamespace("httr", quietly = TRUE)) {
    message("Package 'httr' needed for VRE download; falling back to MOE.")
    vre <- FALSE
  }

  year_end <- as.integer(year_end)
  if (year_end < 2013L)
    stop("`year_end` must be >= 2013 (earliest available ACS 5-year).")

  if (!is.null(census_api_key))
    tidycensus::census_api_key(census_api_key, install = FALSE)

  state_fips  <- .state_to_fips(state)

  # ZCTAs do not support county-level filtering in tidycensus
  if (geography == "zcta" && !is.null(county)) {
    warning("county filtering is not supported for geography = 'zcta'; ",
            "the `county` argument will be ignored.", call. = FALSE)
    county <- NULL
  }
  county_fips <- if (geography == "zcta") NULL else
                   .county_to_fips(county, state_fips)

  # -- 1. Canonicalise variable names ------------------------------
  # Strip trailing "E" suffix users sometimes include (B19013_001E -> B19013_001)
  # Save names BEFORE as.character() since as.character() drops the names attribute
  user_names  <- names(variables)
  variables   <- sub("E$", "", as.character(variables))
  var_names   <- if (!is.null(user_names) && length(user_names) == length(variables))
                   user_names else variables
  var_names   <- make.names(var_names, unique = TRUE)
  named_vars  <- stats::setNames(variables, var_names)

  year_range  <- sprintf("%d\u2013%d", year_end - 4L, year_end)
  if (!quiet)
    message(sprintf("Fetching ACS 5-year %s | geography=%s | state=%s | %d variable(s)...",
                    year_range, geography, state, length(variables)))

  # -- 2. Fetch point estimates (always without geometry for reliability) ------
  acs_raw <- tryCatch(
    tidycensus::get_acs(
      geography = geography,
      variables = named_vars,
      state     = state_fips,
      county    = county_fips,
      year      = year_end,
      survey    = "acs5",
      output    = "wide",
      geometry  = FALSE,
      moe_level = moe_level,
      quiet     = quiet
    ),
    error = function(e) stop("get_acs() failed: ", conditionMessage(e))
  )

  n_units <- nrow(acs_raw)
  if (!quiet) message(sprintf("  Retrieved %d %s(s).", n_units, geography))

  # -- 3. Extract estimate and MOE matrices ------------------------
  # tidycensus wide output: columns named {var_name}E and {var_name}M
  est_cols <- paste0(var_names, "E")
  moe_cols <- paste0(var_names, "M")

  # Detect actual column names (some tidycensus versions suffix differently)
  actual_cols <- if (inherits(acs_raw, "sf"))
    names(sf::st_drop_geometry(acs_raw)) else names(acs_raw)

  missing_e <- setdiff(est_cols, actual_cols)
  if (length(missing_e) > 0)
    stop("Expected estimate columns not found: ",
         paste(missing_e, collapse = ", "),
         ". Check that variable codes are valid for year ", year_end, ".")

  drop_geom <- if (inherits(acs_raw, "sf")) sf::st_drop_geometry else identity

  est_mat <- as.matrix(drop_geom(acs_raw)[, est_cols, drop = FALSE])
  moe_mat <- as.matrix(drop_geom(acs_raw)[, moe_cols, drop = FALSE])
  se_mat  <- .moe_to_se(moe_mat, moe_level)

  colnames(est_mat) <- var_names
  colnames(moe_mat) <- var_names
  colnames(se_mat)  <- var_names
  rownames(est_mat) <- rownames(moe_mat) <- rownames(se_mat) <- acs_raw$GEOID

  # -- 4. VRE download ---------------------------------------------
  replicates_list <- vector("list", length(var_names))
  names(replicates_list) <- var_names
  vre_source <- "moe"

  if (vre) {
    if (!quiet) message("Fetching VRE replicate tables...")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    tables_needed <- .vars_to_tables(variables)
    vre_cache     <- list()   # keyed by table name

    for (tbl in tables_needed) {
      vre_df <- .download_vre(year_end, state_fips, tbl, cache_dir, quiet,
                               timeout_sec)
      if (!is.null(vre_df)) vre_cache[[tbl]] <- vre_df
    }

    if (length(vre_cache) > 0L) {
      n_sdr <- 0L

      for (j in seq_along(variables)) {
        var_j  <- variables[j]
        tbl_j  <- .vars_to_tables(var_j)
        vre_df <- vre_cache[[tbl_j]]
        if (is.null(vre_df)) next

        # Align rows to our tract GEOIDs
        geoid_col  <- .vre_geoid_col(vre_df)
        vre_df[[geoid_col]] <- as.character(vre_df[[geoid_col]])

        # VRE GEOIDs may have "1400000US" prefix -- strip to 11-digit FIPS
        vre_df[[geoid_col]] <- sub("^\\d+US", "", vre_df[[geoid_col]])
        row_idx <- match(acs_raw$GEOID, vre_df[[geoid_col]])

        # Identify replicate columns
        rep_cols <- .vre_rep_cols(vre_df, var_j)

        if (!is.null(rep_cols) && length(rep_cols) == 80L &&
            sum(!is.na(row_idx)) > 0L) {

          rep_mat <- matrix(NA_real_, nrow = n_units, ncol = 80L)
          valid   <- !is.na(row_idx)
          rep_raw <- as.matrix(vre_df[row_idx[valid], rep_cols, drop = FALSE])
          storage.mode(rep_raw) <- "double"
          rep_mat[valid, ]   <- rep_raw
          rownames(rep_mat)  <- acs_raw$GEOID
          colnames(rep_mat)  <- paste0("R", seq_len(80L))
          replicates_list[[j]] <- rep_mat

          # Update SE with SDR formula where replicates are non-NA
          valid_rows <- !is.na(est_mat[, j]) & rowSums(!is.na(rep_mat)) == 80L
          if (sum(valid_rows) > 0L) {
            se_mat[valid_rows, j] <- compute_sdr_se(
              est_mat[valid_rows, j],
              rep_mat[valid_rows, , drop = FALSE]
            )
            n_sdr <- n_sdr + 1L
          }
        }
      }

      if (n_sdr > 0L) {
        vre_source <- "sdr"
        if (!quiet)
          message(sprintf("  SDR SE computed for %d/%d variable(s).",
                          n_sdr, length(variables)))
      } else {
        if (!quiet)
          message("  VRE columns not matched; using MOE-derived SE.")
      }
    } else {
      if (!quiet)
        message("  VRE download failed for all tables; using MOE-derived SE.")
    }
  }

  # -- 4b. Fetch geometry separately (when requested) ----------------------
  # Tries tidycensus first, then tigris (which uses a local cache).
  geom_sf <- NULL
  if (geometry) {
    # Attempt 1: tidycensus
    geom_sf <- tryCatch(
      tidycensus::get_acs(
        geography = geography,
        variables = named_vars[1L],
        state     = state_fips,
        county    = county_fips,
        year      = year_end,
        survey    = "acs5",
        geometry  = TRUE,
        quiet     = quiet
      ),
      error = function(e) {
        if (!quiet)
          message("  tidycensus geometry failed; trying tigris cache. ",
                  "Error: ", conditionMessage(e))
        NULL
      }
    )
    # Attempt 2: tigris (has its own local cache via options(tigris_use_cache))
    if (is.null(geom_sf)) {
      geom_sf <- tryCatch(
        .fetch_tigris_geometry(geography, state_fips, county_fips, year_end),
        error = function(e) {
          if (!quiet)
            message("  tigris geometry also failed. Error: ", conditionMessage(e))
          NULL
        }
      )
    }
    if (is.null(geom_sf) && !quiet)
      message("  Geometry unavailable; $sf will be a plain data.frame.")
  }

  # -- 5. Augment sf with SE (and replicate columns) ---------------
  se_df           <- as.data.frame(se_mat)
  names(se_df)    <- paste0(var_names, "_se")
  acs_out         <- cbind(acs_raw, se_df)

  for (j in seq_along(var_names)) {
    if (!is.null(replicates_list[[j]])) {
      rep_df        <- as.data.frame(replicates_list[[j]])
      names(rep_df) <- paste0(var_names[j], "_R", seq_len(80L))
      acs_out       <- cbind(acs_out, rep_df)
    }
  }

  # Attach geometry if available (makes acs_out an sf object)
  if (!is.null(geom_sf) && inherits(geom_sf, "sf")) {
    idx     <- match(acs_out$GEOID, geom_sf$GEOID)
    acs_out <- sf::st_sf(
      acs_out,
      geometry = sf::st_geometry(geom_sf)[idx],
      crs      = sf::st_crs(geom_sf)
    )
  }

  # -- 6. Assemble output ------------------------------------------
  structure(
    list(
      sf         = acs_out,
      estimates  = est_mat,
      se         = se_mat,
      moe        = moe_mat,
      replicates = replicates_list,
      metadata   = list(
        state      = state,
        state_fips = state_fips,
        county     = county,
        geography  = geography,
        year_end   = year_end,
        year_range = year_range,
        variables  = variables,
        var_names  = var_names,
        n_units    = n_units,
        moe_level  = moe_level,
        vre_source = vre_source,
        call_time  = Sys.time()
      )
    ),
    class = "ACSData"
  )
}


# ================================================================
# 4.  S3 methods for ACSData
# ================================================================

#' Print method for ACSData
#' @param x An `ACSData` object.
#' @param ... Ignored.
#' @export
print.ACSData <- function(x, ...) {
  m <- x$metadata
  has_geom <- inherits(x$sf, "sf") &&
              !all(sf::st_is_empty(sf::st_geometry(x$sf)))
  has_rep  <- !all(vapply(x$replicates, is.null, logical(1L)))

  cat("\u2500\u2500 ACSData (TWDeConv) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  Geography : %s  |  %s%s\n",
              m$geography,
              if (is.null(m$county)) "all counties, " else paste0(m$county, " county, "),
              m$state))
  cat(sprintf("  Period    : ACS 5-year %s  (end year %d)\n",
              m$year_range, m$year_end))
  cat(sprintf("  Units     : %d\n", m$n_units))
  cat(sprintf("  Variables : %d\n", length(m$var_names)))
  for (j in seq_along(m$var_names)) {
    has_r <- !is.null(x$replicates[[j]])
    cat(sprintf("    [%d] %-20s  (code: %s%s)\n",
                j, m$var_names[j], m$variables[j],
                if (has_r) ", VRE \u2713" else ""))
  }
  cat(sprintf("  SE source : %s  |  MOE level: %d%%\n",
              m$vre_source, m$moe_level))
  cat(sprintf("  Geometry  : %s\n",
              if (has_geom) "attached (sf)" else "not available"))
  cat(sprintf("  Replicates: %s\n",
              if (has_rep) "80 per variable (SDR)" else "none (MOE fallback)"))
  cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  invisible(x)
}

# ================================================================
# 5.  prep_acs_for_deconv()  (exported)
# ================================================================

#' Prepare ACS Data for Spatio-Temporal Deconvolution
#'
#' @description
#' High-level pipeline wrapper that fetches multiple years of ACS 5-year
#' estimates, aligns tracts across windows, and assembles all matrices
#' required by `ADMMSolver$new()`.
#'
#' @param state Character. State name, 2-letter abbreviation, or FIPS code.
#' @param county Character or `NULL`. County name, partial name, or 3-digit
#'   FIPS. `NULL` retrieves all units in the state. Ignored for
#'   `geography = "zcta"`.
#' @param year_ends Integer vector of ACS window end-years (must be sorted
#'   ascending). Each entry corresponds to one ACS 5-year window
#'   (e.g., `2022` for ACS 2018-2022). Must be >= 2013.
#' @param variables Single ACS variable code (deconvolution is univariate per
#'   call), e.g. `"B19013_001"`. May be named.
#' @param geography Character. Geographic unit passed to
#'   [tidycensus::get_acs()].  Common values: `"tract"` (default),
#'   `"block group"`, `"county"`, `"zcta"`.  VRE-based SDR standard errors
#'   are only available for `"tract"`; all others fall back to MOE-derived SE.
#' @param k Integer. Temporal difference order for the penalty matrix D.
#'   Default `2L`.
#' @param queen Logical. Use Queen contiguity for the spatial Laplacian?
#'   Default `TRUE`.
#' @param vre Logical. Attempt VRE download for SDR standard errors?
#'   Default `TRUE`. Automatically disabled for non-tract geographies.
#' @param use_srvyr Reserved parameter (unused). Default `TRUE`.
#' @param cache_dir Character. Cache directory for VRE downloads.
#'   Default `tempdir()`.
#' @param census_api_key Character or `NULL`. Census API key.
#' @param quiet Logical. Suppress progress messages. Default `FALSE`.
#'
#' @return A named list of class `"DeconvInput"` with components:
#' \describe{
#'   \item{`Y`}{n  x  S numeric matrix of ACS estimates (units  x  windows).}
#'   \item{`M`}{S  x  T sparse convolution matrix from
#'     [build_acs_convolution_matrix()].}
#'   \item{`W`}{(nS)  x  (nS) sparse diagonal precision weight matrix.}
#'   \item{`D`}{(T-k)  x  T sparse temporal difference matrix.}
#'   \item{`Ls`}{n  x  n sparse symmetric spatial Laplacian.}
#'   \item{`metadata`}{Named list: `state`, `county`, `geography`,
#'     `year_ends`, `variables`, `var_name`, `n_units`, `T_years`,
#'     `S_windows`, `k`, `first_latent_year`, `last_latent_year`,
#'     `vre_sources`.}
#' }
#'
#' @examples
#' \dontrun{
#' inp <- prep_acs_for_deconv(
#'   state     = "IL",
#'   county    = "Cook",
#'   year_ends = 2019:2022,
#'   variables = c(med_inc = "B19013_001"),
#'   k         = 2L,
#'   cache_dir = "~/acs_cache"
#' )
#' print(inp)
#'
#' # Pass directly to the ADMM solver:
#' solver <- ADMMSolver$new(
#'   Y = inp$Y, M = inp$M, W = inp$W, D = inp$D, Ls = inp$Ls,
#'   lambda_t = 10, lambda_s = 1, rho = 1, ridge = 1e-4
#' )
#' }
#'
#' @seealso [fetch_acs_vre_data()], [build_acs_convolution_matrix()],
#'   [build_precision_weights()], [build_empirical_laplacian()],
#'   [build_temporal_penalty()]
#' @export
prep_acs_for_deconv <- function(state,
                                county         = NULL,
                                year_ends,
                                variables,
                                geography      = "tract",
                                k              = 2L,
                                queen          = TRUE,
                                vre            = TRUE,
                                use_srvyr      = TRUE,
                                cache_dir      = tempdir(),
                                census_api_key = NULL,
                                quiet          = FALSE) {

  # -- 1. Validate inputs ----------------------------------------
  year_ends <- as.integer(year_ends)
  if (length(year_ends) < 1L)
    stop("`year_ends` must have at least one element.")
  if (is.unsorted(year_ends))
    stop("`year_ends` must be sorted in ascending order.")
  if (length(variables) != 1L)
    stop("`variables` must be a single ACS variable code ",
         "(deconvolution is univariate per call).")

  k          <- as.integer(k)
  start_year <- year_ends[1L]
  end_year   <- year_ends[length(year_ends)]
  S          <- length(year_ends)

  # -- 2. Fetch ACS data for each window ------------------------
  if (!quiet) message(sprintf("Fetching ACS data for %d window(s) ...", S))

  acs_list    <- vector("list", S)
  vre_sources <- character(S)
  names(vre_sources) <- as.character(year_ends)

  for (i in seq_len(S)) {
    yr <- year_ends[i]
    if (!quiet) message(sprintf("  Window %d/%d: year_end = %d", i, S, yr))
    acs_list[[i]] <- fetch_acs_vre_data(
      state          = state,
      county         = county,
      year_end       = yr,
      variables      = variables,
      geography      = geography,
      geometry       = (i == 1L),   # geometry only for first call
      vre            = vre,
      cache_dir      = cache_dir,
      census_api_key = census_api_key,
      quiet          = quiet
    )
    vre_sources[i] <- acs_list[[i]]$metadata$vre_source
  }

  # -- 3. Align GEOIDs across windows ---------------------------
  geoid_sets    <- lapply(acs_list, function(a) rownames(a$estimates))
  common_geoids <- Reduce(intersect, geoid_sets)

  if (length(common_geoids) == 0L)
    stop("No common tracts found across all ACS windows.")

  n_dropped <- length(geoid_sets[[1L]]) - length(common_geoids)
  if (n_dropped > 0L)
    warning(n_dropped,
            " tract(s) dropped due to missing data in one or more windows.",
            call. = FALSE)

  # -- 4 & 5. Build Y and SE_mat --------------------------------
  var_name <- acs_list[[1L]]$metadata$var_names[1L]

  Y      <- matrix(NA_real_, nrow = length(common_geoids), ncol = S)
  SE_mat <- matrix(NA_real_, nrow = length(common_geoids), ncol = S)
  rownames(Y) <- rownames(SE_mat) <- common_geoids
  colnames(Y) <- colnames(SE_mat) <- paste0("window_", year_ends)

  for (i in seq_len(S)) {
    Y[, i]      <- acs_list[[i]]$estimates[common_geoids, var_name]
    SE_mat[, i] <- acs_list[[i]]$se[common_geoids, var_name]
  }

  # -- 6. Build W ------------------------------------------------
  W <- build_precision_weights(SE_mat)

  # -- 7. Build Ls -----------------------------------------------
  sf_aligned <- acs_list[[1L]]$sf

  # If geometry wasn't attached during ACS fetch (e.g. TIGER server was down),
  # try tigris directly as a second independent attempt.
  if (!inherits(sf_aligned, "sf")) {
    if (!quiet)
      message("  Geometry not in ACS result; trying tigris directly...")
    state_fips_ls  <- acs_list[[1L]]$metadata$state_fips
    geom_raw <- tryCatch(
      .fetch_tigris_geometry(geography, state_fips_ls, county, year_ends[1L]),
      error = function(e) {
        if (!quiet)
          message("  tigris also failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(geom_raw) && inherits(geom_raw, "sf"))
      sf_aligned <- geom_raw
  }

  if (!inherits(sf_aligned, "sf"))
    stop("Geometry unavailable: re-run with geometry=TRUE (the default) ",
         "or ensure the Census TIGER server is reachable.")

  # Subset to common GEOIDs preserving order
  sf_aligned <- sf_aligned[match(common_geoids, sf_aligned$GEOID), ]
  Ls <- build_empirical_laplacian(sf_aligned, queen = queen)

  # -- 8. Build M ------------------------------------------------
  if (all(diff(year_ends) == 1L)) {
    # Consecutive windows: use the dedicated function
    M <- build_acs_convolution_matrix(start_year, end_year)
  } else {
    # Non-consecutive: build sparse M manually with arbitrary column positions
    first_latent <- min(year_ends) - 4L
    last_latent  <- max(year_ends)
    T_nc         <- last_latent - first_latent + 1L
    latent_years <- first_latent:last_latent

    i_idx <- rep(seq_len(S), each = 5L)
    j_idx <- unlist(lapply(seq_len(S), function(s) {
      window_years <- (year_ends[s] - 4L):year_ends[s]
      match(window_years, latent_years)
    }))

    M <- Matrix::sparseMatrix(
      i    = i_idx,
      j    = j_idx,
      x    = 1 / 5,
      dims = c(S, T_nc)
    )
    rownames(M) <- paste0("window_", year_ends)
    colnames(M) <- paste0("year_",   latent_years)
  }

  T_years         <- ncol(M)
  first_latent_yr <- as.integer(sub("year_", "", colnames(M)[1L]))
  last_latent_yr  <- as.integer(sub("year_", "", colnames(M)[T_years]))

  # -- 9. Build D ------------------------------------------------
  D <- build_temporal_penalty(T_years, k)

  # -- 10. Assemble and return -----------------------------------
  structure(
    list(
      Y  = Y,
      M  = M,
      W  = W,
      D  = D,
      Ls = Ls,
      metadata = list(
        state             = state,
        county            = county,
        geography         = geography,
        year_ends         = year_ends,
        variables         = as.character(variables),
        var_name          = var_name,
        n_units           = length(common_geoids),
        T_years           = T_years,
        S_windows         = S,
        k                 = k,
        first_latent_year = first_latent_yr,
        last_latent_year  = last_latent_yr,
        vre_sources       = vre_sources
      )
    ),
    class = "DeconvInput"
  )
}


#' Print method for DeconvInput
#' @param x A `DeconvInput` object.
#' @param ... Ignored.
#' @export
print.DeconvInput <- function(x, ...) {
  m <- x$metadata
  cat("\u2500\u2500 DeconvInput (TWDeConv) \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  Geography : %s  |  %s%s\n",
              m$geography,
              if (is.null(m$county)) "all counties, "
              else paste0(m$county, " county, "),
              m$state))
  cat(sprintf("  Variable  : %s  (code: %s)\n",
              m$var_name, m$variables))
  cat(sprintf("  Windows   : %d  (ACS end years %d \u2013 %d)\n",
              m$S_windows, m$year_ends[1L],
              m$year_ends[length(m$year_ends)]))
  cat(sprintf("  Latent yrs: %d  (%d \u2013 %d)\n",
              m$T_years, m$first_latent_year, m$last_latent_year))
  cat(sprintf("  Units (n) : %d\n", m$n_units))
  cat(sprintf("  Y dims    : %d \u00d7 %d  (units \u00d7 windows)\n",
              nrow(x$Y), ncol(x$Y)))
  cat(sprintf("  M dims    : %d \u00d7 %d  (windows \u00d7 latent years)\n",
              nrow(x$M), ncol(x$M)))
  cat(sprintf("  W dims    : %d \u00d7 %d  (diagonal precision)\n",
              nrow(x$W), ncol(x$W)))
  cat(sprintf("  D dims    : %d \u00d7 %d  (order-%d temporal diff)\n",
              nrow(x$D), ncol(x$D), m$k))
  cat(sprintf("  Ls dims   : %d \u00d7 %d  (spatial Laplacian)\n",
              nrow(x$Ls), ncol(x$Ls)))
  cat(sprintf("  SE source : %s\n",
              paste(unique(m$vre_sources), collapse = ", ")))
  cat("\n  Pass to solver:\n")
  cat("    solver <- ADMMSolver$new(\n")
  cat("      Y=inp$Y, M=inp$M, W=inp$W, D=inp$D, Ls=inp$Ls,\n")
  cat("      lambda_t=..., lambda_s=..., rho=1, ridge=1e-4\n")
  cat("    )\n")
  cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  invisible(x)
}


# ================================================================
# 6.  S3 methods for ACSData
# ================================================================

#' Summary method for ACSData
#' @param object An `ACSData` object.
#' @param ... Ignored.
#' @export
summary.ACSData <- function(object, ...) {
  m    <- object$metadata
  est  <- object$estimates
  se   <- object$se

  cat("\u2500\u2500 ACSData Summary \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  cat(sprintf("  %d %s(s)  |  %s  |  %s\n",
              m$n_units, m$geography, m$year_range, m$state))
  cat("\n  Estimates:\n")
  for (j in seq_along(m$var_names)) {
    vals <- est[, j]
    vals <- vals[is.finite(vals)]
    cat(sprintf("    %-20s  min=%9.1f  median=%9.1f  max=%9.1f  NA=%d\n",
                m$var_names[j],
                min(vals), stats::median(vals), max(vals),
                sum(!is.finite(est[, j]))))
  }
  cat("\n  Standard Errors [", m$vre_source, "]:\n", sep = "")
  for (j in seq_along(m$var_names)) {
    vals <- se[, j]
    vals <- vals[is.finite(vals)]
    cat(sprintf("    %-20s  median SE=%8.1f  CV(median)=%.3f\n",
                m$var_names[j],
                stats::median(vals),
                stats::median(vals) / abs(stats::median(object$estimates[, j],
                                                         na.rm = TRUE) + 1e-9)))
  }
  cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")
  invisible(object)
}
