# TWDeConv

<!-- badges: start -->
[![R-CMD-check](https://github.com/YungChun0303/TWDeConv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YungChun0303/TWDeConv/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**TWDeConv** is an R package for spatio-temporal deconvolution of multi-year
moving-average estimates. Its primary application is recovering latent *annual*
signals from published **ACS 5-year** estimates (which are 5-year moving
averages over Census tracts), but the engine is general and works for any
regularly or irregularly gridded spatial domain with overlapping temporal
windows.

## How it works

The ACS 5-year estimate for tract *i* and window ending in year *s* is
(approximately):

$$Y_{is} = \frac{1}{5}\sum_{t=s-4}^{s} X_{it} + \varepsilon_{is}$$

where **X** is the unobserved latent annual signal. **TWDeConv** recovers
**X** by solving a Kronecker-structured spatio-temporal generalised lasso via
ADMM:

$$\min_{X} \;\tfrac{1}{2}\| W^{1/2}(\text{vec}(Y) - (I_n \otimes M)\text{vec}(X^\top))\|_2^2 + \lambda_t \| (I_n \otimes D)\text{vec}(X^\top)\|_1 + \lambda_s \| (L_s \otimes I_T)\text{vec}(X^\top)\|_1$$

- **M** — S × T convolution matrix encoding the 5-year averaging
- **W** — diagonal precision matrix (inverse ACS sampling variances)
- **D** — *k*-th order temporal difference matrix (trend-filtering penalty)
- **L_s** — spatial graph Laplacian (contiguity-based smoothness penalty)

## Installation

```r
# Install from GitHub (requires remotes)
remotes::install_github("YungChun0303/TWDeConv")
```

The package requires R ≥ 4.1. A Census API key is needed for ACS data
fetching; register for free at <https://api.census.gov/data/key_signup.html>
and install with `tidycensus::census_api_key("YOUR_KEY", install = TRUE)`.

## Quick start — synthetic benchmark

```r
library(TWDeConv)

# 1. Generate synthetic data on a 10×10 grid, T = 20 years
gen <- SyntheticDataGenerator$new(
  n_row = 10, n_col = 10, T_years = 20, window = 5,
  spatial_range = 3, temporal_rho = 0.7, noise_sd = 0.05, seed = 42
)
gen$simulate()          # grid → latent signal → convolution matrix → observed Y

# 2. Build penalty matrices
D  <- build_temporal_penalty(gen$params$T_years, k = 2)
Ls <- build_spatial_laplacian(gen$grid, queen = TRUE)
W  <- build_precision_weights(matrix(0.05, nrow(gen$Y), ncol(gen$Y)))

# 3. Recover the latent annual signal
X_hat <- solve_deconv(
  Y = gen$Y, M = gen$M, W = W, D = D, Ls = Ls,
  lambda_t = 0.5, lambda_s = 0.1, rho = 1.0
)
dim(X_hat)   # 100 × 20

# 4. Evaluate over a (λ_t, λ_s) grid
suite <- run_simulation_suite(
  gen,
  lambda_t_grid = 10^seq(-2, 1, length.out = 6),
  lambda_s_grid = 10^seq(-2, 1, length.out = 6),
  D = D, Ls = Ls
)
print(suite)            # best RMSE / CCC

# 5. Visualise one cell
plot_cell_signals(gen, suite$best$X_hat)
plot_lambda_heatmap(suite)
```

## Quick start — ACS 5-year data (Cook County, IL)

```r
library(TWDeConv)

# Fetch + prepare matrices for two consecutive ACS windows (2021, 2022)
inp <- prep_acs_for_deconv(
  state     = "IL",
  county    = "Cook",
  year_ends = 2019:2022,                     # four 5-year windows
  variables = c(med_inc = "B19013_001"),      # median household income
  k         = 2L,
  cache_dir = "~/acs_cache"                  # persistent VRE cache
)
print(inp)   # dimensions, latent year range, SE source

# Solve
X_hat <- solve_deconv(
  Y = inp$Y, M = inp$M, W = inp$W, D = inp$D, Ls = inp$Ls,
  lambda_t = 5, lambda_s = 0.5, rho = 1.0
)
dim(X_hat)   # n_tracts × T_years
```

## Package modules

| Module | File | Key exports |
|--------|------|-------------|
| Synthetic data | `R/synthetic_data.R` | `SyntheticDataGenerator` |
| Penalty matrices | `R/penalty_matrices.R` | `build_temporal_penalty`, `build_spatial_laplacian`, `build_empirical_laplacian`, `build_precision_weights`, `build_acs_convolution_matrix` |
| ADMM solver | `R/admm_solver.R` | `ADMMSolver`, `solve_deconv` |
| Simulation suite | `R/simulation_suite.R` | `run_simulation_suite`, `plot_cell_signals`, `plot_lambda_heatmap` |
| ACS ingestion | `R/acs_connect.R` | `fetch_acs_vre_data`, `compute_sdr_se`, `prep_acs_for_deconv` |
| ACS weights | `R/acs_weights.R` | `build_acs_precision_weights` |

## Citation

If you use **TWDeConv** in published research, please cite:

```
Li, Y. (2025). TWDeConv: Spatio-temporal deconvolution of moving-average
estimates. R package version 0.1.0.
https://github.com/YungChun0303/TWDeConv
```

## License

MIT © 2025 Yungchun Li
