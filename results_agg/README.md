# results_agg/ — Aggregated Simulation Results

This folder stores the final aggregated simulation results used by the analysis reports. These files are produced by the `data_prep/` pipeline and are committed to git so analyses can be reproduced without re-running the full simulation.

## Files

| File | Description |
|------|-------------|
| `combined_simulation_results.RData` | Combined results across all simulation batches; the primary input for most analysis scripts. |
| `processed_simulation_results_agg_full_factor_2025-01-20.RData` | Aggregated results from the January 2025 full-factorial simulation run. |
| `processed_simulation_results_agg_full_factor_2025-04-03.RData` | Aggregated results from the April 2025 full-factorial simulation run. |
| `processed_simulation_results_agg_full_factor_2025-05-21.RData` | Aggregated results from the May 2025 full-factorial simulation run (most recent). |

Each `.RData` file contains aggregated performance metrics (bias, RMSE, coverage rates, SE ratios, etc.) organized by estimator, estimand, and simulation scenario factors.
