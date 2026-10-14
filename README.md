# Code for replicating a comparison of methods for estimate average treatment effects in Blocked, Cluster Randomized Controlled Trials (RCTs)

This repository contains the code for running both a simulation study and empirical analysis for comparing ATE estimators in blocked, cluster-randomized trials (CRTs). In particular, the code in this repository was used to generate the results presented in "An Applied Researcher's Guide to Estimating Effects from Multisite Cluster Randomized Trials: Estimands, Estimators, and Estimates," a paper accepted to JREE in 2026. Note that this repository is a collection of scripts, not an R package in any traditional sense.

The scripts in this repository do rely on a package `clusterRCT`, also on GitHub (at <https://github.com/lmiratrix/clusterRCT>). That package implements all the methods tested, and has various exploratory tools for people analysing blocked, cluster randomized experiments (or cluster randomized experiments) of their own. That package might be useful, in other words, for everyday work.

This collection of scripts is useful for verifying the results of our paper, running simulations to see how the results hold up in different contexts, and for just generating synthetic multilevel data.

## Reproducing the Paper

The paper depends on a lot of complex analyses. We have used Claude to help generate some pointers and guidance for reproducing or understanding the analyses done. You can read through the structure below, but also see [REPLICATION.md](/REPLICATION.md) in this repository for more information about what code was used to generate what result in the main paper.

**A note on exact reproducibility.** Our goal in sharing this code is to let others understand our approach, adapt it, and run similar simulations of their own. Unfortunately, we were not able to rerun everything using the current code base, so the code will likely not reproduce the paper's reported numbers to the decimal. This is a research codebase that evolved over a multi-year project, and we unfortunately don't have a reliable record of exactly which version of everything (in particular, which commit of the companion `clusterRCT` package, itself under active development) produced the final numbers. Rerunning the pipeline today should give results that are qualitatively and substantively the same, but we can't guarantee an exact match to every table and figure. See [REPLICATION.md](/REPLICATION.md) for more detail on this.

## Repository Structure


### Main Analyses Code

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`analysis/`](analysis/README.md) | RMarkdown reports and scripts for analyzing the results of the large multi-factor simulation studies conducted. Primary output: `Simulation_WriteUp.pdf` (the paper's simulation supplement). All of these scripts primarily use the final, aggregated simulation results (see below). |
| [`empirical/`](empirical/README.md) | Applies all estimators to a collection of real CRT datasets. Produces the empirical results reported in the paper and generates figures. Some of the underlying individual-level data lives on a secure server; this folder stores processing scripts and exported results, and not the raw data. Other of the underlying data are publicly available in the [Item Response Warehouse](https://redivis.com/datasets/as2e-cv7jb41fd) (IRW, [Domingue et al., 2024](https://link.springer.com/article/10.3758/s13428-025-02796-y)) under the prefix “gilbert_meta”; these data could be analyzed after prepping them using code in  `prep_gilbert_data/` |
| [`prep_gilbert_data/`](prep_gilbert_data/README.md) | Processes that public IRW “gilbert_meta” collection end to end: item-level response data → per-outcome analytic datasets → `compare_methods()`/`describe_clusterRCT()` results and Appendix A descriptive tables. |
|                                                     |                                                              |
| [`simulation/`](simulation/) | The code to actually run the main multifactor simulation. Main simulation entry-point scripts: `generate_grids.R`, `generate_grids_lowdf.R`, `estimator_performance_simulation.R`, `estimator_performance_simulation_fixed_district.R`, and `run_single_simulation.R`. |
| [`simulationII/`](simulationII/README.md) | Supplementary simulation that uses a rerandomization apporoach to closely adhere to the existing data structure. |

### Simulation Infrastructure

| Folder/File | Description |
|------------------------------------|------------------------------------|
| [`R/`](R/README.md) | Shared R functions as a series of easily sourcable R scripts. These scripts include data-generating processes, simulation runners, aggregation utilities, and MCSE helpers. Sourced by scripts throughout the repo. |
|                                     |                                                              |
| [`cloud/`](cloud/README.md) | SLURM batch scripts for running the large-scale simulation on the Harvard FASRC (Cannon) cluster under the `miratrix_lab` account. |
| [`data_prep/`](data_prep/README.md) | Pipeline (3 scripts) that bundles raw cluster job outputs into cleaned, aggregated `.RData` files ready for analysis. |
| `scenarios.rds` | The used scenario grid file, that indexes all the simulations run. |

### Results

| Folder | Description |
|----------------------------|--------------------------------------------|
| [`results/`](results/README.md) | Folder to store scenario grid `.rds` files and raw simulation result files from cloud runs. Large per-job output files are not committed to git, so the git pull of this directory will be fairly empty. |
| [`results_agg/`](results_agg/README.md) | Aggregated, processed simulation results (these final results are committed to git). These files are the primary input for analysis reports. |

### Other

| Folder/File | Description |
|------------------------------------|------------------------------------|
| [`demo_scripts/`](demo_scripts/README.md) | Short demonstration scripts for key workflows. Useful for understanding the code base. |
| [`figures/`](figures/README.md) | Generated PDF figures for the paper are saved here. |
| `install_needed_packages.R` | Script to install all R packages required to run the simulations. |


------------------------------------------------------------------------

## Typical Workflows

### Running a small simulation locally

``` r
source("install_needed_packages.R")                                      # install dependencies
source("simulation/generate_grids.R")                                    # create scenarios.rds
source("simulation/estimator_performance_simulation_fixed_district.R")   # run sim
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

