# Robust and Adaptive Optimization — A Quantitative Framework for Bond Futures Basis Strategies

**Imperial College London — Electrical and Electronic Engineering — Final Year Project (2026)**  
**Author:** Lucas Lasanta | **Supervisor:** Prof. Tony Constantinides

---

## Overview

This repository contains the full MATLAB implementation of a two-layer Lagrange Programming Neural Network (LPNN) framework for bond futures basis trading across four U.S. Treasury maturities (2-year, 5-year, 10-year and 30-year).

The strategy exploits mean-reversion in the CTD-adjusted futures-cash basis. Layer 1 computes a robust yield-curve hedge ratio for each maturity. Layer 2 allocates capital across the four hedged basis trades using a robust min-max objective with ellipsoidal factor uncertainty, a Tikhonov smoothness penalty and a transaction cost term, implemented as a warm-started primal-dual LPNN dynamical system. A CVX interior-point solver solves the identical problem at each date as a benchmark.

The framework achieves gross Sharpe ratios of 1.77 (LPNN) and 1.79 (CVX) over the 2019-2025 out-of-sample test period, falling to 1.76 and 1.78 after realistic institutional transaction costs of 0.1 basis points, with a maximum drawdown of 1.42%. The LPNN-CVX return correlation is 0.9998 at the production configuration, rising to 1.0000 under the adaptive step recommended by the discretisation analysis.

---

## Requirements

| Requirement | Version tested |
|---|---|
| MATLAB | R2021a or later |
| Statistics and Machine Learning Toolbox | Included with most MATLAB licences |
| CVX (with SeDuMi solver) | Free academic licence — cvxr.com/cvx |

CVX is only required for the `benchmarks.m` module (CVX robust and mean-variance benchmarks). The LPNN itself runs entirely on native MATLAB with no external dependencies.

---

## Data Requirements

The strategy requires Bloomberg data that is **not included** in this repository due to data licensing restrictions. The following CSV files must be placed in the same folder as the `.m` files before running:

| File | Contents |
|---|---|
| `model_data2.csv` | Weekly basis returns, repo rate, yield curve changes (10 tenors) |
| `futures.csv` | Daily front-month futures prices for TU, FV, TY, US |
| `CTD_selected.csv` | Daily CTD bond identification, conversion factors, DV01 |
| `basis.csv` | Daily gross CTD-adjusted basis levels for all four maturities |

All data was sourced from Bloomberg terminals. TU (2-year) basis data is only consistently available from July 2021 onwards; the framework handles this automatically using a factor-model fallback signal for earlier dates.

---

## File Structure

| File | Purpose |
|---|---|
| `main.m` | Master script — runs the full pipeline in order |
| `load_data.m` | Reads and aligns the four CSV data files |
| `run_pca.m` | Estimates PCA loadings on training period only (no look-ahead) |
| `calibrate_params.m` | Grid search for gamma and kappa on the validation period (2017-2018) |
| `backtest.m` | Rolling weekly backtest over the test period (2019-2025) |
| `layer1_hedge.m` | Layer 1: robust hedge ratio via LPNN (one per maturity) |
| `layer2_portfolio.m` | Layer 2: robust portfolio allocation via LPNN primal-dual system |
| `benchmarks.m` | CVX robust and mean-variance benchmarks at each rebalancing date |
| `lw_shrink.m` | Ledoit-Wolf analytical covariance shrinkage |
| `evaluate.m` | Performance metrics, signal diagnostics and output figures |
| `tikhonov_frontier.m` | Tikhonov frontier sweep across lambda_T values |
| `sensitivity_analysis.m` | Sensitivity analysis across key parameter configurations |
| `transaction_cost_analysis.m` | Ex-post transaction cost stress test |
| `discretisation_analysis.m` | Discretisation comparison: Fixed Euler, RK4, Adaptive Euler, Euler above smooth bound |
| `plot_signal_flow_graph.m` | Signal flow graph of the Layer-2 LPNN primal-dual system |

---

## How to Run

1. Clone or download this repository
2. Place the four Bloomberg CSV files in the same folder as the `.m` files
3. Open MATLAB and navigate to that folder
4. Open `main.m` and press Run (or type `main` in the command window)

The pipeline runs automatically: data loading, PCA, parameter calibration, backtest, evaluation, discretisation analysis, the adaptive-step verification backtest and the iteration-count convergence table. Runtime is approximately one and a half hours in this full configuration, dominated by the CVX benchmark solves inside each backtest. Commenting out the verification backtests at the bottom of `main.m` reduces this to roughly 20 minutes for the core backtest and evaluation.

---

## Parameter Reference

All parameters are set in `main.m` inside the `p` struct. The key ones are:

| Parameter | Value used | Meaning |
|---|---|---|
| `p.window` | 104 | Rolling estimation window (weeks) |
| `p.K` | 4 | Number of factors (3 PCA + repo) |
| `p.tau` | 0.1 | LPNN step size (20% of the stability bound tau_c = 0.500) |
| `p.c` | 1 | Augmented Lagrangian penalty coefficient |
| `p.n_iter` | 2500 | LPNN iterations per rebalancing date |
| `p.lambda_T` | 0.01 | Tikhonov smoothness penalty strength |
| `p.w_min / p.w_max` | -1 / +1 | Portfolio weight bounds |
| `p.chi2_alpha` | 0.95 | Confidence level for ellipsoidal uncertainty set |
| `p.ridge_lam` | 1e-4 | Ridge regularization for factor loading estimation |
| `p.basis_window` | 52 | Basis z-score rolling window (weeks) |
| `p.gamma` | Calibrated | Risk aversion (calibrated on validation period) |
| `p.kappa` | Calibrated | Transaction cost penalty (calibrated on validation period) |

`p.gamma` and `p.kappa` are calibrated automatically by `calibrate_params.m` on the 2017-2018 validation period. The calibrated values are gamma = 0.001 and kappa = 0.0001, with a best validation Sharpe of 1.671.

Two notes on the solver configuration. The critical step size tau_c = 2/(c*M) = 0.500 is an exact structural constant of the model parameters rather than a market quantity, since the Tikhonov Laplacian has zero row sums and contributes nothing along the budget direction. The discretisation analysis recommends the adaptive step tau = 0.45 (90% of tau_c) for any redeployment; the production results retain tau = 0.10, with the adaptive configuration verified separately in `main.m`.

---

## Output Files

Running `main.m` produces the following figures saved in the working directory:

| File | Content |
|---|---|
| `fig_cumulative_returns` | Cumulative return paths for all four strategies |
| `fig_drawdowns` | Drawdown series over the test period |
| `fig_lpnn_weights` | LPNN portfolio weights by maturity |
| `fig_active_weights` | Active weights relative to equal weight |
| `fig_activity` | Gross exposure and L1 distance from equal weight |
| `fig_lyapunov_energy` | LPNN Lyapunov energy convergence diagnostic |
| `fig_hedge_ratios` | Layer 1 hedge ratios (confirms duration-neutrality) |
| `fig_epsilon` | Ellipsoidal uncertainty radius over time |
| `fig_mu_hat` | Basis z-score expected return signal by maturity |
| `fig_tikhonov_frontier` | Sharpe, turnover, exposure and drawdown vs lambda_T |
| `fig_transaction_costs` | Net Sharpe and net return vs transaction cost level |
| `fig_discretisation` | Best-so-far Lyapunov energy for four discretisation schemes |

---

## Reproducing Specific Analyses

The analytical extensions are controlled by commenting or uncommenting lines near the bottom of `main.m`:

```matlab
% Analytical extensions (each runs multiple full backtests):
% sensitivity_analysis(returns_w, F, dates_w, dv01, basis_levels, p);  % ~2 hours
% tikhonov_frontier(returns_w, F, dates_w, dv01, basis_levels, p);     % ~3 hours
% transaction_cost_analysis(results);                                  % seconds
% plot_signal_flow_graph();                                            % seconds
```

Each analysis is self-contained and can be run after the main backtest has completed, with the `results` struct in the workspace. The discretisation analysis and the adaptive-step verification backtest run by default as part of `main.m`.

---

## Notes

- All reported results are strictly out-of-sample. The test period (January 2019 to December 2025) was not used for any parameter selection.
- PCA loadings are estimated on the training period only (up to end of 2016) and held fixed throughout, avoiding look-ahead bias in the factor construction.
- The LPNN warm-starts from the previous week's state at each rebalancing date. The `results` struct from `backtest.m` must be in the workspace if individual analyses are run separately after the main backtest.
- Results do not constitute investment advice.

---

## Reference

Zhang, S. and Constantinides, A.G. (1992). Lagrange Programming Neural Networks. IEEE Transactions on Circuits and Systems II, 39(7), 441-452.
