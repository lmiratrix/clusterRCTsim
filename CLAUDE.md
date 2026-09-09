# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

This is a research repository (not an R package) containing the simulation study and empirical analysis code for comparing ATE estimators in blocked, cluster-randomized trials (CRTs). It is the simulation/analysis companion to the [`clusterRCT`](https://github.com/lmiratrix/clusterRCT) R package, which implements the estimators themselves. Scripts throughout the repo `source()` shared functions from `R/` rather than loading a package — there is no `DESCRIPTION`/`NAMESPACE`, so this is not `devtools::load_all()`-able.

Start with `README.md` at the repo root — it has an up-to-date table of every folder's purpose and the canonical local/cluster workflows. Most subfolders (`R/`, `analysis/`, `empirical/`, `simulationII/`, `reports/`, `covariate_simulations/`, `data_prep/`, `cloud/`, `results/`, `results_agg/`, `demo_scripts/`) have their own `README.md` — read the relevant one before working in that folder rather than re-deriving structure from file listing.

## Setup

```r
source("install_needed_packages.R")
```

Installs `tidyverse`, `randomizr`, `estimatr`, `lme4`, `lmerTest`, `multiwayvcov`, and the two GitHub-only dependencies `lmiratrix/clusterRCT` and `lmiratrix/blkvar` (via `devtools::install_github(..., force = TRUE)`).

## Tests

Tests use testthat but are run standalone (no package to `check()`), against `tests/testthat/`:

```r
library(testthat)
test_dir(here::here("tests/testthat/"))
```

or simply `source("tests/testthat.R")`. To run a single test file, `source()` the relevant `R/` file(s) it depends on first (each test file sources its own dependencies at the top), then `testthat::test_file("tests/testthat/test-DGP.R")`. All path resolution in this repo uses `here::here()`, so R's working directory should be the repo root (or an RStudio project session opened via `clusterRCTsim.Rproj`).

Coverage: `test-DGP.R` (data generation), `test-simulation.R` (superpopulation single-CRT runner), `test-simulation_fixed_block.R` (fixed-blocks runner), `test-generate_data_from_fingerprint.R` (fingerprint-based data generation).

## Architecture

**Two DGP/simulation frameworks, used for different purposes:**
- *Superpopulation framework* (`R/simulation_support_code_crt.R`, `single.CRT()`): regenerates block structure fresh on every simulation run.
- *Fixed-blocks framework* (`R/simulation_support_code_fixed_block.R`): generates block characteristics once, then holds them constant across repeated rerandomizations. This is what the main cloud/grid simulation (`simulation/estimator_performance_simulation_fixed_district.R`) uses.

Both build on `R/data_generation_code.R`, the core DGP: generates cluster sizes with variable ICC, treatment-effect heterogeneity, and correlations between cluster size/treatment proportion/impact.

A third mode, **fingerprint simulation** (`R/calibrated_simulation_code.R`, `R/generate_data_from_fingerprint.R`), generates synthetic data that mirrors the structure (block/cluster IDs, sample sizes) of a specific real dataset, combined with user-supplied ICC/R² parameters — used to ask "how would estimators behave on data like this real study?" A related but distinct mode, rerandomization simulation (in `simulationII/`), fixes one generated dataset and repeatedly rerandomizes treatment assignment.

**End-to-end pipeline for the main grid simulation:**
1. `simulation/generate_grids.R` (or `generate_grids_lowdf.R`) builds the scenario grid → `scenarios.rds`. Prior grids are snapshotted to `archive/` (dated) before regenerating, since `scenarios.rds` is overwritten in place — see `archive/README` note in the root README for the convention.
2. Code + `scenarios.rds` pushed to the Harvard FASRC (Cannon) cluster; `cloud/run_it.slurm` runs a SLURM array job, one task per grid chunk, calling `simulation/estimator_performance_simulation_fixed_district.R`. `cloud/execute_runfullsim.sh` orchestrates submit → wait → post-process.
3. `data_prep/` is a strict 3-step pipeline run in order (or via `attic/fully_prepare_simulation_results.R`): `bundle_odyssey_results.R` (combine per-job `.rds` → `complete_odyssey_results.rds`) → `clean_simulation_results.R` (drop failed/malformed runs, unnest, produce `scenarios`/`simdata`/`sruns` in `processed_simulation_results.RData`) → `Combine_Results_Across_Simulations.R` (merge batches into `results_agg/`). Note: some scripts in this folder have "conflicted copy" variants from a Dropbox sync issue — always use the non-conflicted version.
4. Aggregation math (bias, RMSE, coverage, MCSE) lives in `R/aggregate_simulation_results.R` and `R/MCSE_functions.R`; `R/group_estimators.R` collapses estimators that produce numerically identical estimates/SEs, useful before reporting since several named estimators coincide algebraically in some regimes.
5. `analysis/Simulation_WriteUp.Rmd` is the primary report/paper supplement consuming `results_agg/`.

**Empirical side** (`empirical/`) applies the same `clusterRCT` estimators to real CRT datasets whose underlying individual-level data lives only on a secure server (`/data/share/pse/IES CRT/`) — not in this repo. `empirical/master_analyze_all_data.R` is the master entry point; fingerprint simulations *on* these empirical datasets live in `simulationII/` (not `empirical/`) since they're part of that separate workflow, despite depending on empirical data-loading code.

**Large/derived files are git-partitioned deliberately**: raw per-job cluster output (`results/simulation_results_*.rds`) is NOT committed (regenerate via cluster reruns); processed/aggregated outputs in `results_agg/` ARE committed so reports are reproducible without rerunning the full simulation.

## Conventions

- All file paths are resolved with `here::here()` from the repo root — don't hardcode paths relative to a script's own location.
- Simulation scripts typically start with `rm(list = ls())` and end by writing a `.rds`/`.RData` output — treat them as standalone entry points, not sourceable libraries (unlike everything in `R/`, which is meant to be sourced).
- When modifying `simulation/generate_grids.R`, remember the archive convention: copy the current `generate_grids.R` + `scenarios.rds` pair into `archive/` with a creation-date suffix before regenerating, so prior grids used for submitted results are never silently lost.
