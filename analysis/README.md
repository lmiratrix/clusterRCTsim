# analysis/ — Main Simulation Analysis

This folder contains the reports and analysis scripts for the large multi-factor simulation study. The primary output is the simulation supplement accompanying the paper.

## Reports & Scripts

| File | Output | Description |
|------|--------|-------------|
| `Simulation_WriteUp.Rmd` | `Simulation_WriteUp.pdf` | **Main simulation supplement.** Full write-up of the multi-factor simulation study comparing ATE estimators across blocked, cluster-randomized designs. This is the document submitted/distributed as a paper supplement. |
| `load_and_prep_data.R` |  | Loads aggregated simulation results from `results_agg/` and prepares data frames for analysis. Sourced by most report Rmds. |

## 

## Shared Utility Functions (`R/`)

Shared plotting/tree-fitting helpers sourced by the reports in `analysis/`, `empirical/`, and `simulationII/`. Unlike `R/` at the repo root (DGP and simulation-runner code), these are analysis/reporting utilities only.

## 

| File | Purpose |
|------|---------|
| `create_analysis_tree.R` | Fits an rpart tree with a `num_leaves` pruning target, a legacy column `name_map`, and a signed color scale. Used by `Simulation_WriteUp.Rmd` and `empirical/empirical_results_report.Rmd`. |
| `create_analysis_tree_minmax.R` | A separate, non-interchangeable rpart tree fitter with a `min_leaves`/`max_leaves` pruning API. Used by `simulationII/fingerprint_simulation_report.Rmd`. See the header comments in both files for why these are kept distinct rather than merged. |
| `create_custom_boxplot.R` | Shared boxplot styling helper used across the simulation supplement. |

**Note on the two tree functions:** `create_analysis_tree()` and `create_analysis_tree_minmax()` look like duplicates of each other (same purpose, same underlying `rpart` engine) but have genuinely different pruning APIs and are each load-bearing for a different report. They are deliberately kept as separate files/functions rather than merged — see the header comment in each file.
