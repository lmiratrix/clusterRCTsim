# empirical/ — Empirical Analysis on Real CRT Datasets

This folder contains all code for applying the clusterRCT estimators to a collection of real cluster-randomized trial datasets. The goal is to show how different estimators compare on actual data and to characterize the range of estimates they produce.

The underlying data lives on various secure servers and is not stored in this repository. Clean-room datasets are processed there; results are exported to `empirical/results/`.

**Reproducibility without server access:** the empirical report (`empirical_results_report.Rmd`, via `load_and_prep_empirical_results.R`) reads from the committed, aggregated `empirical/results/*.RDS` files, not from raw data — so the main text's empirical figures/tables ARE reproducible without secure-server access. Regenerating those aggregated `.RDS` files themselves (via `master_analyze_all_data.R`) requires server access.



## Reports

| File | Output | Description |
|------------------|---------------------|----------------------------------|
| `empirical_results_report.Rmd` | `empirical_results_report.pdf` | **Main empirical report.** Documents results reported in the paper: summary statistics on studies, range of estimates across methods, tables and figures. |
| `nonblocked_experiment_report.Rmd` | `nonblocked_experiment_report.pdf` | Parallel analysis restricted to non-blocked (completely randomized) CRT studies. |
| `make_study_flashcards.Rmd` | `make_study_flashcards.pdf` | One-page summaries ("flashcards") for each study in the dataset. |
|                                    |                                    |                                                              |

## Core Analysis Scripts

| File | Purpose |
|-----------------------------|-------------------------------------------|
| `master_analyze_all_data.R` | **Master script to be run on raw data.** Applies `compare_methods()` and `describe_clusterRCT()` to all non-clean-room datasets, combines with clean-room output, and saves  to `results/master241204_combined.RDS`. _This cannot be run without the original data_. |
| `load_empirical_results.R` | Loads raw empirical results from the `results/` subdirectory (these are the sets of ATE estimates and standard errors for all the methods on all the datasets). |
| `load_and_prep_empirical_results.R` | Loads and pre-processes empirical results, calculates range statistics and flags. Sourced by the report Rmds. |
| `calculate_all_ranges.R` | Calculates the range of estimates (spread across methods) for each study outcome. |
| `data_processing_functions.R` | Utility functions for reading and cleaning individual study datasets for doing analysis. |

## Supporting Scripts

| File | Purpose |
|-----------------------------|-------------------------------------------|
| `final_paper_figures.R` | Generates the publication-quality figures for the paper (blocked studies). |
| `final_paper_figures_nonblocked.R` | Generates figures for the non-blocked studies. |
| `make_estimator_feature_table.R` | Builds the table summarizing estimator properties. |
| `R/makeDotplot.R` | Creates dot-plot visualizations of estimator comparisons. |
| `check_quality_control.R` | Quality control checks on the empirical data and results. |

## Results Subdirectory (`results/`)

Contains exported analysis results:

\- `master241204_combined.RDS` — Combined results from all datasets (Dec 2024 run)

\- `master240528_combined.RDS` — Earlier combined results (May 2024 run)

\- `jg_results.csv` — Results for Josh Gilbert's study collection.

\- `study_descriptions.csv` — Summary descriptions of all included studies, produced by `load_and_prep_empirical_results.R`. Read by `simulationII/load_master_rerandomization_results.R` (Simulation II's data-loading step, needed to render `analysis/Simulation_WriteUp.Rmd`) - committed for that reason, unlike most other intermediate empirical outputs.

\- Excel files with estimator comparison tables
