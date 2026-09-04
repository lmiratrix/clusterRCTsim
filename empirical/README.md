# empirical/ — Empirical Analysis on Real CRT Datasets

This folder contains all code for applying the clusterRCT estimators to a collection of real cluster-randomized trial datasets. The goal is to show how different estimators compare on actual data and to characterize the range of estimates they produce.

The underlying data lives on a secure server (`/data/share/pse/IES CRT/`) and is not stored in this repository. Clean-room datasets are processed there; results are exported to `empirical/results/`. Scripts that touch this data read the root path from `IES_CRT_ROOT <- Sys.getenv("IES_CRT_DATA_ROOT", "/data/share/pse/IES CRT")` — set the `IES_CRT_DATA_ROOT` environment variable if your mount point differs.

**Reproducibility without server access:** the empirical report (`empirical_results_report.Rmd`, via `load_and_prep_empirical_results.R`) reads from the committed, aggregated `empirical/results/*.RDS` files, not from raw data — so the main text's empirical figures/tables ARE reproducible without secure-server access. Only regenerating those aggregated `.RDS` files themselves (via `master_analyze_all_data.R`) requires server access.

> **Note:** Scripts that run fingerprint simulations on these datasets (`master_run_fingerprints.R`, `master_run_fingerprints_testing.R`, `testing_fingerprint_simR.R`) live in [`fingerprint/`](../fingerprint/README.md), since they are part of the fingerprint simulation workflow.

## Reports

| File | Output | Description |
|------------------|---------------------|----------------------------------|
| `empirical_results_report.Rmd` | `empirical_results_report.pdf` | **Main empirical report.** Documents results reported in the paper: summary statistics on studies, range of estimates across methods, tables and figures. |
| `nonblocked_experiment_report.Rmd` | `nonblocked_experiment_report.pdf` / `.docx` | Parallel analysis restricted to non-blocked (completely randomized) CRT studies. |
| `make_study_flashcards.Rmd` | `make_study_flashcards.pdf` | One-page summaries ("flashcards") for each study in the dataset. |
| `standard_error_exploration.Rmd` | — | Exploratory look at standard error behavior across the empirical studies. |

## Core Analysis Scripts

| File | Purpose |
|-----------------------------|-------------------------------------------|
| `master_analyze_all_data.R` | **Master script.** Applies `compare_methods()` and `describe_clusterRCT()` to all non-clean-room datasets, combines with clean-room output, and saves `results/master241204_combined.RDS`. |
| `load_empirical_results.R` | Loads raw empirical results from the `results/` subdirectory. |
| `load_and_prep_empirical_results.R` | Loads and pre-processes results, calculates range statistics and flags. Sourced by the report Rmds. |
| `calculate_all_ranges.R` | Calculates the range of estimates (spread across methods) for each study outcome. |
| `data_processing_functions.R` | Utility functions for reading and cleaning individual study datasets. |

## Supporting Scripts

| File | Purpose |
|-----------------------------|-------------------------------------------|
| `final_paper_figures.R` | Generates the publication-quality figures for the paper (blocked studies). |
| `final_paper_figures_nonblocked.R` | Generates figures for the non-blocked studies. |
| `make_estimator_feature_table.R` | Builds the table summarizing estimator properties. |
| `makeDotplot.R` | Creates dot-plot visualizations of estimator comparisons. |
| `excel_supp.R` | Produces Excel supplementary tables. |
| `josh_study_exploration.R` | Exploration of datasets from Josh Gilbert's study collection. |
| `check_discrepancy_DB_v_Agg.R` | Checks for discrepancies between design-based and aggregate estimators. |
| `check_quality_control.R` | Quality control checks on the empirical data and results. |
| `random_result_exploration.R` | Ad hoc exploration of selected results. |
| `beaconmentor_test.R` | Test script for the Beacon Mentor study. |

## Results Subdirectory (`results/`)

Contains exported analysis results:

\- `master241204_combined.RDS` — Combined results from all datasets (Dec 2024 run)

\- `master240528_combined.RDS` — Earlier combined results (May 2024 run)

\- `jg_results.csv`, `jg_describe.csv` — Results for Josh Gottlieb's study collection

\- `study_descriptions.csv` — Summary descriptions of all included studies - Excel files with estimator comparison tables
