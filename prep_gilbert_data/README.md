# prep_gilbert_data/ — Processing the Public (IRT Warehouse) CRT Collection

This folder processes the collection of public cluster-randomized trial datasets described in the paper (Gilbert, Himmelsbach, Soland, Joshi, & Domingue 2025; publicly available via the Item Response Warehouse, prefix `gilbert_meta`, at <https://redivis.com/datasets/as2e-cv7jb41fd>). 

The script turns item-level response data into per-outcome analytic datasets and builds the descriptive study tables used in Appendix A. 

The rerandomization simulation (Simulation II) that depends on these data, and everything downstream of that, now lives in [`simulationII/`](../simulationII/README.md) alongside its restricted-data counterpart — see that folder's README. This is the *only* dataset source behind the paper's empirical results that isn't restricted-access — everything else in `empirical/` is not sharable.



## Pipeline

Run in order:

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `01_prepare_item_level_data.R` | Item-level `.Rdata` export (see below — **not included**) | `data/raw/crt_data_{id}.csv`, one row per person per outcome, Rasch-scored |
| 2 | `load_all_datasets.R` | `study_metadata.xlsx` + `data/raw/crt_data_{id}.csv` | (not a standalone step — sourced by steps 3-5; loads everything into memory and applies documented per-dataset fixes for known data issues) |
| 3 | `02_compare_and_describe.R` | (sources step 2) | `data/clean/jg_compare.csv`, `jg_describe.csv`, `jg_results.csv` |
| 4 | `03_make_study_tables.R` | `study_metadata.xlsx` + `data/raw/crt_data_{id}.csv` | `tabfig/study_summaries.xlsx`, `study_long.xlsx` (Appendix A descriptive tables) |
| — | `check_results_integrity.R` | `../empirical/results/jg_results.csv` | console-only sanity check, nothing saved |

`processing_functions.R` is the shared function library (`read_in()`, `grab_data()`, `data_processor()`, `cm()`/`dc()` wrappers around `compare_methods()`/`describe_clusterRCT()`) sourced by steps 2-3; it isn't a pipeline step itself.

**Running Simulation II on this data** (one of the possible next steps after this folder's pipeline) is `simulationII/04_run_rerandomization_simulation.R` — it is there since it runs the simulation rather than preparing data, alongside `simulationII/05_analyze_rerand_simulation.R` (a sanity check on its output) and `simulationII/master_run_rerandomization.R` (the equivalent runner for the restricted/curated studies). See `simulationII/README.md`.

### Step 1 is not independently runnable here

`01_prepare_item_level_data.R` starts from a `.Rdata` export of item-level response data that exists only on the original processor's machine — it isn't in this repo, isn't public, and isn't recoverable from anything else here. The script is included for **provenance/transparency** (this is exactly how `data/raw/crt_data_{id}.csv` was produced) but will `stop()` immediately unless you set `JG_ITEM_DATA_PATH` to your own copy of that export. In practice, treat `data/raw/crt_data_{id}.csv` as the actual starting point for this pipeline (see "What's not committed" below).



## What's not committed

Following this repo's convention (raw/large data isn't committed; small aggregated results are — see the root `README.md`), `data/raw/`, `data/clean/`, `data/new_datasets2025/`, and `tabfig/` are gitignored. The two files the rest of the repo actually depends on (`jg_results.csv`, `jg_describe.csv`) are committed at their consumption point, `empirical/results/`, not here — see the `if (FALSE)` block at the end of `02_compare_and_describe.R` for how they get there.

This means: **the scripts and `study_metadata.xlsx` here document and can regenerate the pipeline, but the actual `data/raw/crt_data_*.csv` files (37MB) and `data/new_datasets2025/gilbert_meta_*.csv` files (17MB) are not currently included.** Since this data is public and step 1 can't regenerate it independently (see above), a replicator without those files can't currently rerun this pipeline end-to-end — only reuse its already-committed output (`empirical/results/jg_results.csv` etc.) and `simulationII/`'s permutation results. Whether to commit those raw CSVs here so the pipeline is fully self-contained is an open call — flagging it rather than deciding it unilaterally.
