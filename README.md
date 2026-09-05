# clusterRCTsim — Simulation Study for Cluster RCT Estimators

This repository contains the simulation study and empirical analysis code for comparing ATE estimators in blocked, cluster-randomized trials (CRTs). In particular, the code in this repository was used to generate the results presented in "An Applied Researcher's Guide to Estimating Effects from Multisite Cluster Randomized Trials: Estimands, Estimators, and Estimates," a paper accepted to JREE in 2026. Note that this repository is a collection of scripts, not an R package in any traditional sense.

The scripts in this repository do rely on a package `clusterRCT`, also on GitHub (at <https://github.com/lmiratrix/clusterRCT>). That package implements all the methods tested, and has various exploratory tools for people analysing blocked, cluster randomized experiments (or cluster randomized experiments) of their own.

## Repository Structure


### Main Analyses

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`analysis/`](analysis/README.md) | Reports and scripts for the large multi-factor simulation study. Primary output: `Simulation_WriteUp.pdf` (the paper's simulation supplement). Also contains scripts for targeted analyses (coverage, factor importance, estimator comparisons). These scripts primarily use the final, aggregated simulation results. |
| [`empirical/`](empirical/README.md) | Applies all estimators to a collection of real CRT datasets. Produces the empirical results reported in the paper and generates figures. The underlying individual-level data lives on a secure server; this folder stores processing scripts and exported results. |
| [`fingerprint/`](fingerprint/README.md) | Supplementary simulations, in particular the "Fingerprint simulations" that use DGPs calibrated to the real datasets and "rerandomization simulations" that use the real datasets but permute treatment assignment to assess estimator stability. |

### Simulation Infrastructure

| Folder/File | Description |
|------------------------------------|------------------------------------|
| [`R/`](R/README.md) | Shared R functions: data-generating processes, simulation runners, aggregation utilities, and MCSE helpers. Sourced by scripts throughout the repo. |
| [`simulation/`](simulation/) | Main simulation entry-point scripts: `generate_grids.R`, `generate_grids_lowdf.R`, `estimator_performance_simulation.R`, `estimator_performance_simulation_fixed_district.R`, and `run_single_simulation.R`. |
| [`cloud/`](cloud/README.md) | SLURM batch scripts for running the large-scale simulation on the Harvard FASRC (Cannon) cluster under the `miratrix_lab` account. |
| [`data_prep/`](data_prep/README.md) | Pipeline (3 scripts) that bundles raw cluster job outputs into cleaned, aggregated `.RData` files ready for analysis. |
| `scenarios.rds` | Current scenario grid file (also mirrored in `results/`). |

### Results

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`results/`](results/README.md) | Scenario grid `.rds` files and raw simulation result files from cloud runs. Large per-job output files are not committed to git. |
| [`results_agg/`](results_agg/README.md) | Aggregated, processed simulation results (committed to git). Primary input for analysis reports. |

### Other

| Folder/File | Description |
|------------------------------------|------------------------------------|
| [`demo_scripts/`](demo_scripts/README.md) | Short demonstration scripts for key workflows (fingerprint sims, grid DGP, single-scenario runs). |
| [`figures/`](figures/README.md) | Generated PDF figures for the paper. |
| `install_needed_packages.R` | Installs all R packages required to run the simulations. |


------------------------------------------------------------------------

## Typical Workflows

### Running a small simulation locally

``` r
source("install_needed_packages.R")                                      # install dependencies
source("simulation/generate_grids.R")                                    # create scenarios.rds
source("simulation/estimator_performance_simulation_fixed_district.R")  # run sim
source("data_prep/clean_simulation_results.R")                           # clean output
```

### Running the full simulation on the cluster

1.  Push code and `scenarios.rds` to a computer cluster.
2.  Submit the simulation jobs via `cloud/execute_runfullsim.sh` or `sbatch cloud/run_it.slurm`.
3.  After jobs complete, run `cloud/results.slurm` to bundle outputs.
4.  Pull results to your local machine and run the `data_prep/` pipeline locally.
5.  Open `analysis/Simulation_WriteUp.Rmd` to regenerate the report.

### Reproducing empirical results

``` r
source("empirical/master_analyze_all_data.R")   # requires server data access
knit("empirical/empirical_results_report.Rmd")
```

## Acknowledgements

The current data generating code is adapted from two sources: the PUMP package simulation code and the data generation code found in the [Miratrix, Weiss & Henderson multisite estimands, estimators and estimates](https://www.tandfonline.com/doi/full/10.1080/19345747.2020.1831115) paper.


