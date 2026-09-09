# analysis/R/ — Analysis Utility Functions

Shared plotting/tree-fitting helpers sourced by the reports in `analysis/`, `empirical/`, and `simulationII/`. Unlike `R/` at the repo root (DGP and simulation-runner code), these are analysis/reporting utilities only.

## Files

| File | Purpose |
|------|---------|
| `create_analysis_tree.R` | Fits an rpart regression tree to simulation/empirical output with a `num_leaves` pruning target, a legacy column `name_map`, and a signed color scale. Also defines `create_analysis_tree_v2()`, a `partykit::ctree`-based variant. Used by `analysis/Simulation_WriteUp.Rmd`, `empirical/empirical_results_report.Rmd`, `analysis/small_block_simulation_analysis.Rmd`, `analysis/factor_significance_exploration.R`. |
| `create_analysis_tree_minmax.R` | A separate rpart tree fitter with a `min_leaves`/`max_leaves` pruning API (`create_analysis_tree_minmax()`). Used only by `simulationII/fingerprint_simulation_report.Rmd`. |
| `create_custom_boxplot.R` | Shared boxplot styling helper used across the simulation supplement. |

## Note on the two tree functions

`create_analysis_tree()` and `create_analysis_tree_minmax()` look like duplicates of each other (same purpose, same underlying `rpart` engine) but have genuinely different pruning APIs and are each load-bearing for a different report. They are deliberately kept as separate files/functions rather than merged — see the header comment in each file. Don't consolidate them without checking both reports' output first.
