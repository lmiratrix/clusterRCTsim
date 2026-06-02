# fingerprint/ — Fingerprint and Rerandomization Simulations

This folder contains code and reports for two related simulation approaches:

1. **Fingerprint simulations**: Synthetic simulations calibrated to match the structural characteristics of a specific real dataset (its "fingerprint" — block structure, cluster sizes, ICCs, etc.). Used to understand how estimators would behave on data *like* that study.

2. **Rerandomization simulations**: Simulations that fix a generated dataset and repeatedly rerandomize treatment assignment, comparing how variance estimators perform under this finite-population perspective.

## Reports

| File | Output | Description |
|------|--------|-------------|
| `fingerprint_simulation_report.Rmd` | `fingerprint_simulation_report.pdf` | **Main fingerprint simulation report.** Summarizes fingerprint simulation results across the real-data studies. |
| `fingerprint_simulation_report_rerand_shadow.Rmd` | `fingerprint_simulation_report_rerand_shadow.pdf` | Shadow report that presents rerandomization simulation results in the same format as the fingerprint report, enabling direct comparison. |
| `rerandomization_simulation_report.Rmd` | `rerandomization_simulation_report.pdf` | Main rerandomization simulation report: compares SE estimators under repeated rerandomization of fixed data. |
| `rerandomization_sim_report_within.Rmd` | `rerandomization_sim_report_within.pdf` | Rerandomization analysis for within-study (individual-level outcome) designs. |

## Key Scripts

| File | Purpose |
|------|---------|
| `master_run_fingerprints.R` | **Master runner.** Loads all real datasets from the secure server and runs a fingerprint simulation on each, saving per-study `.rds` files to `results_fingerprint/`. |
| `master_run_fingerprints_testing.R` | Test/development version of `master_run_fingerprints.R` for iterating on the pipeline. |
| `full_fingerprint_simulation.R` | Runs the full set of fingerprint simulations across all datasets (older/alternate entry point). |
| `fingerprint_simulation_demo.R` | Short demo of how fingerprint simulations work; useful for understanding the approach. |
| `testing_fingerprint_simR.R` | Tests and validates the fingerprint simulation code on individual empirical datasets. |
| `load_master_fingerprint_results.R` | Loads and combines fingerprint simulation results from `results_fingerprint/` for use in reports. |
| `load_master_permutation_results.R` | Loads and combines rerandomization simulation results from `results_rerand*/` for use in reports. |

## Results Subdirectories

| Subdirectory | Contents |
|---|---|
| `results_fingerprint/` | Per-study `.rds` files from fingerprint simulations (one file per study/outcome, named by study ID). |
| `results_rerand/`, `results_rerand_master/`, `results_rerand_sim_public/` | Per-study `.rds` files from rerandomization simulations under various configurations. |
| `figs/` | Generated figure files (PDF) from the reports. |
| `attic/` | Older/superseded versions of analysis scripts. |
| `scratch/` | Scratch scripts used during development. |
| `exploratory_sims/` | Exploratory simulation scripts predating the current reports. |
