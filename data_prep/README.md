# data_prep/ — Simulation Results Processing Pipeline

This folder contains a three-step pipeline that takes raw output files from the Harvard FASRC cluster and produces cleaned, aggregated data files used in analysis.

## Pipeline Steps

Run these scripts in order (or source `attic/fully_prepare_simulation_results.R` to run all three):

### Step 1: `bundle_odyssey_results.R`
Reads all individual `simulation_results_*.rds` files from the `results/` folder (one per cluster job) and combines them into a single `complete_odyssey_results.rds` file. Also reports run-time diagnostics.

### Step 2: `clean_simulation_results.R`
Loads `complete_odyssey_results.rds`, removes failed or malformed simulation runs, unnests the nested results, and adds metadata. Produces three data frames saved to `processed_simulation_results.RData`:
- `scenarios` — one row per simulation scenario
- `simdata` — one row per simulation run per scenario (can be very large)
- `sruns` — aggregated statistics by method and simulation run

### Step 3: `Combine_Results_Across_Simulations.R`
Combines results across multiple simulation batches (e.g., different scenario grids or time periods) into unified data files saved in `results_agg/`.

## Notes

- `bundle_odyssey_results.R` and `clean_simulation_results.R` each have a "conflicted copy" variant from a Dropbox sync issue — use the non-conflicted versions.
- The raw per-job result files are not committed to git (too large); they must be retrieved from the cluster or re-run.
