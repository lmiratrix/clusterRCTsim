# prep_gilbert_data/ — Processing the Public (IRT Warehouse) CRT Collection

This folder processes the collection of public cluster-randomized trial datasets described in the paper (Gilbert, Himmelsbach, Soland, Joshi, & Domingue 2025; publicly available via the Item Response Warehouse, prefix `gilbert_meta`, at <https://redivis.com/datasets/as2e-cv7jb41fd>). For the purposes of this study, we use the version of the data publicly posted in the replication materials for Gilbert, et al., 2025: <https://doi.org/10.7910/DVN/C4TJCA>.

The script turns item-level response data into per-outcome analytic datasets and builds the descriptive study tables used in Appendix A. The rerandomization simulation (Simulation II) that depends on these data, and everything downstream of that, now lives in [`simulationII/`](../simulationII/README.md) alongside its restricted-data counterpart — see that folder's README. This is the *only* dataset source behind the paper's empirical results that isn't restricted-access — everything else in `empirical/` is not sharable.

These scripts were written jointly by Josh Gilbert and Luke Miratrix.

## Pipeline

Run in order:

| Step | Script | Input | Output |
|------------------|------------------|------------------|------------------|
| 1 | `01_prepare_item_level_data.R` | `datasets_list.Rdata` | `data/raw/crt_data_{id}.csv`, one row per person per outcome, Rasch-scored |
| 2 | `load_all_datasets.R` | `study_metadata.xlsx` + `data/raw/crt_data_{id}.csv` | (not a standalone step — sourced by steps 3-5; loads everything into memory and applies documented per-dataset fixes for known data issues) |
| 3 | `02_compare_and_describe.R` | (sources step 2) | `data/clean/jg_compare.csv`, `jg_describe.csv`, `jg_results.csv` |
| 4 | `03_make_study_tables.R` | `study_metadata.xlsx` + `data/raw/crt_data_{id}.csv` | `tabfig/study_summaries.xlsx`, `study_long.xlsx` (Appendix A descriptive tables) |
| — | `check_results_integrity.R` | `../empirical/results/jg_results.csv` | console-only sanity check, nothing saved |

`processing_functions.R` is the shared function library (`read_in()`, `grab_data()`, `data_processor()`, `cm()`/`dc()` wrappers around `compare_methods()`/`describe_clusterRCT()`) sourced by steps 2-3; it isn't a pipeline step itself.

**Running Simulation II on this data** (one of the possible next steps after this folder's pipeline) is `simulationII/04_run_rerandomization_simulation.R` — it is there since it runs the simulation rather than preparing data, alongside `simulationII/05_analyze_rerand_simulation.R` (a sanity check on its output) and `simulationII/master_run_rerandomization.R` (the equivalent runner for the restricted/curated studies). See `simulationII/README.md`.
