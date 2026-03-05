# TWDeConv

<!-- badges: start -->
[![R-CMD-check](https://github.com/YungChun0303/TWDeConv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YungChun0303/TWDeConv/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**TWDeConv** is an R package for spatio-temporal deconvolution of multi-year moving-average estimates. 
Its primary application is recovering latent *annual* signals from published **ACS 5-year** estimates (e.g., 5-year moving averages over Census tracts), 
but the engine is general and works for any regularly or irregularly gridded spatial domain with overlapping temporal windows. 
It formalizes the conceptual three-phase approach from the sources into a unified, Kronecker-structured generalized lasso engine

## Core Model: The Spatio-Temporal Engine
The package recovers the unobserved latent annual signal ($X$) by solving a single convex optimization problem. This objective function balances three primary goals:

1. Data Fidelity: It ensures the recovered signal, when averaged over five years (via convolution matrix $\mathbf{M}$, matches the observed ACS estimates, weighted by their precision ($\mathbf{W}$).
2. Temporal Smoothness: It applies a *k*-th order trend-filtering penalty ($\mathbf{D}$) to capture structured paths and identify sudden "level shifts" or changing trends.
3. Spatial Continuity: It utilizes a spatial graph Laplacian ($L_s$) to penalize divergence between adjacent Census tracts, allowing the model to "borrow strength" across space.

The model is solved using the **Alternating Direction Method of Multipliers (ADMM)**, which is the industry standard for handling non-differentiable ($L_1$) penalties while maintaining computational efficiency for large-scale national data


## How it works

The ACS 5-year estimate for tract *i* and window ending in year *s* is (approximately):

$$Y_{is} = \frac{1}{5}\sum_{t=s-4}^{s} X_{it} + \varepsilon_{is}$$

where **X** is the unobserved latent annual signal. **TWDeConv** recovers
**X** by solving a Kronecker-structured spatio-temporal generalised lasso via
ADMM:

$$\min_{X} \;\tfrac{1}{2}\| W^{1/2}(\vec{Y}) - (I_n \otimes M)\vec{X^\top})\|_2^2 + \lambda_t \| (I_n \otimes D)\text{vec}(X^\top)\|_1 + \lambda_s \| (L_s \otimes I_T)\text{vec}(X^\top)\|_1$$

- **$M$** — S × T convolution matrix encoding the 5-year averaging
- **$W$** — diagonal precision matrix (inverse ACS sampling variances)
- **$D$** — *k*-th order temporal difference matrix (trend-filtering penalty)
- **$L_s$** — spatial graph Laplacian (contiguity-based smoothness penalty)

## Key Package Modules and Commands
The package is organized into specialized modules that handle everything from data ingestion to advanced optimization.
1. ACS Data Ingestion and Precision Modeling
These commands implement the project's requirement for "exact" margins of error and precision-weighted least squares.
- `fetch_acs_vre_data()`: Pulls the 80 separate variance replicates from **Variance Replicate Estimate (VRE)** tables.
- `compute_sdr_se()`: Uses **Successive Difference Replication (SDR)** to calculate official-standard standard errors that account for complex survey designs.
- `build_acs_precision_weights()`: Constructs the diagonal precision matrix (**$W$**) used in the objective function to ensure higher-quality estimates have more influence.
- `prep_acs_for_deconv()`: A preprocessing tool that includes a **tigris fallback** to ensure spatial geometry is correctly assigned to the data even if metadata is missing.
2. Matrix Construction and Penalties
These functions build the operators required for the Kronecker-structured formulation.
- `build_temporal_penalty()`: Generates the k-th order discrete difference matrix (D) for trend filtering.
- `build_spatial_laplacian()`: Constructs a contiguity-based smoothness penalty based on the spatial network of Census units.
- `build_acs_convolution_matrix()`: Encodes the 5-year moving average process (matrix M) that links latent annual values to observed data.
3. The ADMM Solver
This is the "engine room" of the package, fulfilling the strategic roadmap’s call for a scalable solver.
- `ADMMSolver`: A class-based implementation of the ADMM algorithm designed to decouple data-fidelity and penalty terms.
- `solve_deconv()`: The primary user-facing command that executes the deconvolution and returns the recovered latent annual signal.
4. Simulation and Validation Suite
Because tract-level annual "ground truth" is often unavailable, these commands allow researchers to test the model on synthetic data.
- `SyntheticDataGenerator`: Creates known annual signals for benchmarking.
- `run_simulation_suite()`: Automates the testing of the model against synthetic datasets using metrics like RMSE and Lin’s CCC.
- `plot_lambda_heatmap()`: Assists in hyperparameter tuning by visualizing how different values of $\lambda_t$(temporal) and $\lambda_s$(spatial) affect recovery accuracy.

This structured implementation allows the TWDeConv package to scale to the massive geographic resolution of the entire United States while maintaining the inferential rigor required for social policy research

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
Chun, Y. and Mei, X. (2026). TWDeConv: Spatio-temporal deconvolution of moving-average
estimates. R package version 0.1.0.
https://github.com/YungChun0303/TWDeConv
```

## License

MIT © 2025 Chun, Yung
