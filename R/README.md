# R/ — Core Simulation Library

This folder contains the shared R functions used throughout the project. These are sourced (not packaged) by the simulation scripts and analysis reports.

## Files

| File | Purpose |
|------|---------|
| `data_generation_code.R` | Data-generating process (DGP): generates cluster sizes with variable ICC, treatment effect heterogeneity, and correlations between cluster size, treatment proportion, and impact. Core building block for all simulations. |
| `simulation_support_code_crt.R` | Functions for running a single CRT simulation trial in the super-population framework (regenerates blocks each run). Includes `single.CRT()` and scenario-list generation. |
| `simulation_support_code_fixed_block.R` | Companion to `simulation_support_code_crt.R` for the fixed-blocks framework, where block characteristics are generated once and held constant across rerandomizations. Used in the main cloud simulations. |
| `calibrated_simulation_code.R` | Functions for fingerprint simulations and rerandomization simulations — calibrated runs that target a specific real dataset by matching its structural fingerprint. |
| `generate_data_from_fingerprint.R` | Generates synthetic data that mirrors the structure (block/cluster IDs, sample sizes) of a real dataset, combined with user-supplied ICC and other parameters. |
| `aggregate_simulation_results.R` | Aggregates raw simulation output into performance metrics (bias, RMSE, coverage, etc.) using superpopulation and finite-sample frameworks. |
| `group_estimators.R` | Groups estimators that produce identical point estimates or variance estimates (within tolerance) across simulation runs, useful for identifying redundancies. |
| `MCSE_functions.R` | Monte Carlo standard error (MCSE) utilities for assessing simulation precision. |
| `generate_nonnegative_integers.R` | Helper to sample non-negative integers with a given mean and SD, used for generating variable cluster sizes. |
| `util.R` | Miscellaneous utility functions (e.g., `scat()` for formatted console output). |

## Seeding convention

- **Grid generation** (`simulation/generate_grids.R`, `generate_grids_lowdf.R`): a single global `set.seed(40444)` at the top of the script, so the scenario grid itself is reproducible.
- **Per-scenario simulation runs** (`simulation/estimator_performance_simulation.R`, `estimator_performance_simulation_fixed_district.R`, `run_single_simulation.R`): each cloud/array job seeds with `set.seed(index^2 + index + 1)`, where `index` (or `SCENARIO_ID`) is the job's scenario index. This gives each scenario a distinct, reproducible seed derived deterministically from its index, rather than reusing one seed across all scenarios or seeding from wall-clock time.
- The only other `set.seed()` calls in `R/` (in `util.R`, `group_estimators.R`, `generate_nonnegative_integers.R`) are fixed constants (`123`, `1039`) inside `if (FALSE) { ... }` demo/testing blocks — they don't affect any simulation run.
