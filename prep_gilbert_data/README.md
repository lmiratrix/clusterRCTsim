# prep_gilbert_data/ — Processing the Public (IRT Warehouse) CRT Collection

This folder processes the collection of public cluster-randomized trial datasets described in the paper (Gilbert, Himmelsbach, Soland, Joshi, & Domingue 2025; publicly available via the Item Response Warehouse, prefix `gilbert_meta`, at <https://redivis.com/datasets/as2e-cv7jb41fd>). It turns item-level response data into per-outcome analytic datasets, runs every estimator on each one, builds the descriptive study tables used in Appendix A, and runs the rerandomization simulation (Simulation II) that feeds `fingerprint/`. This is the *only* dataset source behind the paper's empirical results that isn't restricted-access — everything else in `empirical/` depends on the secure IES CRT server (see `empirical/README.md`).

Adapted from code originally written and run by Josh Gilbert. A handful of purely exploratory/debugging scripts from the original source (`data_processor.R`, an earlier and now-unused superset-subset duplicate of `processing_functions.R`; `add_new.R`, `explore_single_dataset.R`, `explore_specific_datasets.R`, ad hoc one-off explorations; `individual_dataset_code.R`, explicitly marked defunct by its own header comment) were left out since nothing in the pipeline sources them.

## Pipeline

Run in order:

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `01_prepare_item_level_data.R` | Item-level `.Rdata` export (see below — **not included**) | `data/raw/crt_data_{id}.csv`, one row per person per outcome, Rasch-scored |
| 2 | `load_all_datasets.R` | `study_metadata.xlsx` + `data/raw/crt_data_{id}.csv` | (not a standalone step — sourced by steps 3-5; loads everything into memory and applies documented per-dataset fixes for known data issues) |
| 3 | `02_compare_and_describe.R` | (sources step 2) | `data/clean/jg_compare.csv`, `jg_describe.csv`, `jg_results.csv` |
| 4 | `03_make_study_tables.R` | `study_metadata.xlsx` + `data/raw/crt_data_{id}.csv` | `tabfig/study_summaries.xlsx`, `study_long.xlsx` (Appendix A descriptive tables) |
| 5 | `04_run_rerandomization_simulation.R` | (sources step 2) | `../fingerprint/results_rerand_sim_public/jg_crt_{id}.rds`, one file per blocked dataset |
| — | `05_analyze_rerand_simulation.R` | step 5's output | console-only sanity check, nothing saved |
| — | `check_results_integrity.R` | `../empirical/results/jg_results.csv` | console-only sanity check, nothing saved |

`processing_functions.R` is the shared function library (`read_in()`, `grab_data()`, `data_processor()`, `cm()`/`dc()` wrappers around `compare_methods()`/`describe_clusterRCT()`) sourced by steps 2-5; it isn't a pipeline step itself.

### Step 1 is not independently runnable here

`01_prepare_item_level_data.R` starts from a `.Rdata` export of item-level response data that exists only on the original processor's machine — it isn't in this repo, isn't public, and isn't recoverable from anything else here. The script is included for **provenance/transparency** (this is exactly how `data/raw/crt_data_{id}.csv` was produced) but will `stop()` immediately unless you set `JG_ITEM_DATA_PATH` to your own copy of that export. In practice, treat `data/raw/crt_data_{id}.csv` as the actual starting point for this pipeline (see "What's not committed" below).

### A known gap: 7 additional datasets aren't wired in yet

`study_metadata.xlsx` lists dataset IDs 76-82 (Thai et al.; Cabell et al.) as `cluster_rct == 1`, and raw item-level files for them exist at `data/new_datasets2025/gilbert_meta_{76-82}.csv` — but no script here turns them into `data/raw/crt_data_{76-82}.csv`, and they don't yet appear in `empirical/results/jg_results.csv` (which currently tops out at dataset ID 74, 60 distinct datasets total). Whoever picks this up next will need either a new version of step 1 that handles this newer, differently-shaped raw format, or a one-off equivalent of `crt_clean()` written for it.

### A bug fixed in the copy here

The original `02 an_crt.R` referenced an undefined variable `over` while computing `over_d` (`over_d <- setdiff( over, ... )`) — this would throw `object 'over' not found` if the describe/combine section were run fresh. Fixed to `setdiff( over_d, ... )` in `02_compare_and_describe.R` here.

## What's not committed

Following this repo's convention (raw/large data isn't committed; small aggregated results are — see the root `README.md`), `data/raw/`, `data/clean/`, `data/new_datasets2025/`, and `tabfig/` are gitignored. The two files the rest of the repo actually depends on (`jg_results.csv`, `jg_describe.csv`) are committed at their consumption point, `empirical/results/`, not here — see the `if (FALSE)` block at the end of `02_compare_and_describe.R` for how they get there. Likewise, `04_run_rerandomization_simulation.R`'s output writes directly into `fingerprint/results_rerand_sim_public/` (already gitignored, consistent with the rest of `fingerprint/`'s regenerable simulation outputs).

This means: **the scripts and `study_metadata.xlsx` here document and can regenerate the pipeline, but the actual `data/raw/crt_data_*.csv` files (37MB) and `data/new_datasets2025/gilbert_meta_*.csv` files (17MB) are not currently included.** Since this data is public and step 1 can't regenerate it independently (see above), a replicator without those files can't currently rerun this pipeline end-to-end — only reuse its already-committed output (`empirical/results/jg_results.csv` etc.) and `fingerprint/`'s permutation results. Whether to commit those raw CSVs here so the pipeline is fully self-contained is an open call — flagging it rather than deciding it unilaterally.
