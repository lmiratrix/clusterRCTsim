# analysis/ — Main Simulation Analysis

This folder contains the reports and analysis scripts for the large multi-factor simulation study. The primary output is the simulation supplement accompanying the paper.

## Reports

| File | Output | Description |
|------|--------|-------------|
| `Simulation_WriteUp.Rmd` | `Simulation_WriteUp.pdf` | **Main simulation supplement.** Full write-up of the multi-factor simulation study comparing ATE estimators across blocked, cluster-randomized designs. This is the document submitted/distributed as a paper supplement. |
| `small_block_simulation_analysis.Rmd` | `small_block_simulation_analysis.pdf` / `.html` | Analysis focused on scenarios with small numbers of clusters per block (low-df regime). |
| `calculate_paper_results.Rmd` | `calculate_paper_results.pdf` | Earlier version of the simulation analysis; calculates the numerical results cited in the paper. |
| `calculate_paper_results2.Rmd` | — | Extended version of `calculate_paper_results.Rmd`. |
| `simulation_analysis.Rmd` | — | Older simulation analysis (circa 2019), predecessor to the current write-up. |

## Analysis Scripts

| File | Purpose |
|------|---------|
| `load_and_prep_data.R` | Loads aggregated simulation results from `results_agg/` and prepares data frames for analysis. Sourced by most report Rmds. |
| `evaluate_single_scenario.R` | Runs and inspects a single simulation scenario in detail; useful for debugging or targeted exploration. |
| `explore_simulation_scenario_grids.R` | Examines the structure and coverage of the simulation scenario grid. |
| `factor_significance_exploration.R` | Explores which simulation factors most strongly predict estimator performance differences. |
| `linear_regression_on_sim_results.R` | Fits linear models to simulation performance metrics to identify key predictors. |
| `coverage_investigation.R` | Focused investigation of CI coverage rates across scenarios and estimators. |
| `compare_DBSite_to_FIRC.R` | Compares two specific estimator families: design-based site-level and fully interacted random coefficients. |
| `calculate_paper_results.R` | Script version (non-Rmd) of the paper results calculation. |
| `make_tree_of_sim_results.R` | Builds regression trees (via `R/create_analysis_tree.R`) explaining which simulation factors drive bias in specific method subsets. Relocated here from the repo root. |

## Shared Utility Functions (`R/`)

| File | Purpose |
|------|---------|
| `create_analysis_tree.R` | Fits an rpart tree with a `num_leaves` pruning target, a legacy column `name_map`, and a signed color scale. Used by `Simulation_WriteUp.Rmd` and `empirical/empirical_results_report.Rmd`. |
| `create_analysis_tree_minmax.R` | A separate, non-interchangeable rpart tree fitter with a `min_leaves`/`max_leaves` pruning API. Used by `simulationII/fingerprint_simulation_report.Rmd`. See the header comments in both files for why these are kept distinct rather than merged. |
| `create_custom_boxplot.R` | Shared boxplot styling helper used across the simulation supplement. |
