# clusterRCTsim — Simulation Study for Cluster RCT Estimators

This repository contains the simulation study and empirical analysis code for comparing ATE estimators in blocked, cluster-randomized trials (CRTs). The estimators themselves are implemented in the companion [`clusterRCT`](https://github.com/lmiratrix/clusterRCT) R package.

## Repository Structure

### Core Functions

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`R/`](R/README.md) | Shared R functions: data-generating processes, simulation runners, aggregation utilities, and MCSE helpers. Sourced by scripts throughout the repo. |

### Main Analyses

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`analysis/`](analysis/README.md) | Reports and scripts for the large multi-factor simulation study. Primary output: `Simulation_WriteUp.pdf` (the paper's simulation supplement). Also contains scripts for targeted analyses (coverage, factor importance, estimator comparisons). |
| [`empirical/`](empirical/README.md) | Applies all estimators to a collection of real CRT datasets. Produces the empirical results reported in the paper and generates figures. The underlying individual-level data lives on a secure server; this folder stores processing scripts and exported results. |
| [`fingerprint/`](fingerprint/README.md) | Fingerprint simulations (calibrated to real datasets) and rerandomization simulations (fixed dataset, repeated treatment assignment). Reports compare SE estimators under these conditions. |
| [`reports/`](reports/README.md) | Standalone technical reports: superpopulation vs. finite-population SE comparisons, and a getting-started tutorial for new contributors. |
| [`covariate_simulations/`](covariate_simulations/README.md) | Simulations exploring covariate adjustment in blocked CRTs, including L1 vs. L2 covariates and the low-df regime. |

### Simulation Infrastructure

| Folder/File | Description |
|------------------------------------|------------------------------------|
| [`cloud/`](cloud/README.md) | SLURM batch scripts for running the large-scale simulation on the Harvard FASRC (Cannon) cluster under the `miratrix_lab` account. |
| [`data_prep/`](data_prep/README.md) | Pipeline (3 scripts) that bundles raw cluster job outputs into cleaned, aggregated `.RData` files ready for analysis. |
| `generate_grids.R` | Generates the full-factorial simulation scenario grid (`scenarios.rds`) used by the cloud run. |
| `generate_grids_lowdf.R` | Generates the low-df scenario grid variant. |
| `estimator_performance_simulation.R` | Script for running a simulation on a given scenario grid (super-population version). |
| `estimator_performance_simulation_fixed_district.R` | Fixed-district version of the performance simulation. |
| `scenarios.rds` | Current scenario grid file (also mirrored in `results/`). |

### Results

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`results/`](results/README.md) | Scenario grid `.rds` files and raw simulation result files from cloud runs. Large per-job output files are not committed to git. |
| [`results_agg/`](results_agg/README.md) | Aggregated, processed simulation results (committed to git). Primary input for analysis reports. |

### Other

| Folder/File | Description |
|------------------------------------|------------------------------------|
| [`small_simulations/`](small_simulations/README.md) | Small, self-contained simulations exploring specific questions locally. Many predate the current infrastructure. |
| [`demo_scripts/`](demo_scripts/README.md) | Short demonstration scripts for key workflows (fingerprint sims, grid DGP, single-scenario runs). |
| [`explorations/`](explorations/README.md) | Ad hoc scripts for targeted theoretical/methodological questions. |
| [`documentation/`](documentation/README.md) | Early project documentation. |
| [`figures/`](figures/README.md) | Generated PDF figures for the paper. |
| `install_needed_packages.R` | Installs all R packages required to run the simulations. |

------------------------------------------------------------------------

## Typical Workflows

### Running a small simulation locally

``` r
source("install_needed_packages.R")      # install dependencies
source("generate_grids.R")               # create scenarios.rds
source("estimator_performance_simulation_fixed_district.R")  # run sim
source("data_prep/clean_simulation_results.R")               # clean output
```

### Running the full simulation on the cluster

1.  Push code and `scenarios.rds` to FASRC.
2.  Submit via `cloud/execute_runfullsim.sh` or `sbatch cloud/run_it.slurm`.
3.  After jobs complete, run `cloud/results.slurm` to bundle outputs.
4.  Pull results and run the `data_prep/` pipeline locally.
5.  Open `analysis/Simulation_WriteUp.Rmd` to regenerate the report.

### Reproducing empirical results

``` r
source("empirical/master_analyze_all_data.R")   # requires server data access
knit("empirical/empirical_results_report.Rmd")
```

## Acknowledgements

The current data generating code is adapted from two sources: the PUMP package simulation code and the data generation code found in the Miratrix, Weiss & Henderson multisite estimands, estimators and estimates paper.

## Archive

The `archive/` folder holds dated snapshots of the grid-generation program and its output, saved before we change them so earlier versions are never lost.

The 2025-11-20 snapshot is the version used in the initial journal submission:

- `archive/generate_grids_2025-11-20.R`: the `generate_grids.R` that produced the 11/20/2025 grid (three block types: same, reasonable, messy; 1,296 scenarios).
- `archive/scenarios_2025-11-20.rds`: the grid that program generated, and the one the submitted results were run on.

When `generate_grids.R` is rerun it overwrites `scenarios.rds` in place, so we copy the current pair into `archive/` with the grid's creation date before regenerating.

# WARNING: Readme from this point on needs to be updated

## Running the simulation

We ran the simulation on the Odessey cluster at Harvard. The simulation takes several computers awhile to run, but it is easy to run smaller simulations targeting specific scenarios.

The Odyssey scripts will generate a large number of rds files, each file being the result of a fraction of the simulation runs. We then bundle these results into a single file, and then analyze that file.

## Overview of file structure

We have the following

- `simulation_support_code` - Code for generating data and running the simulation. Uses the blkvar package, which implements the various methods we explored.
- `data_prep` folder - Three files that form a pipeline that takes the raw output from the Odyssey cluster and makes some aggregate datafiles. These scripts are described below. There is also a bundling file (`fully_prepare_simulation_results.R`) that just sources these files in order.
- `analysis` folder - Collection of scripts that do the analysis of the simulation results presented in the paper.

## Instructions to get simulation results data ready

The steps are as follows:

1)  Gather results from the cluster into a single file.

```         
source( "prepare simulation results/bundle_odyssey_results.R" )
```

2)  Clean up single file, remove bad sim runs, add in metadata and summary data.

```         
source( "prepare simulation results/clean_simulation_results.R" )
```

3)  Calculate the performance metrics used for the report.

```         
source( "prepare simulation results/prepare_aggregate_results.R" )
```

This final file works with the following data frames:

- scenarios (list of scenario characteristics)
- simdata (the raw data; this can be a very big file)
- sruns (aggregated data by method and simulation run)

## Instructions to get simulation results

Source the `calculate_paper_results.R` file. This is the R code from an old R markdown script that was then converted to the simulation writeup tex document.

Data-generating code is adapted from the [PUMP](https://github.com/MDRCNY/PUMP) package simulation code and from the Miratrix, Weiss & Henderson (2021) multisite estimands paper. The original simulation infrastructure was built by Jonathan Che; the current multi-factor design was developed by Sophie Litschwartz and Luke Miratrix.
