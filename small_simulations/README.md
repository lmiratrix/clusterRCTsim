# small_simulations/ — Small Exploratory Simulations

This folder contains small, self-contained simulation scripts used to explore specific questions or prototype ideas. Unlike the cloud simulations, these run locally in minutes and target narrow questions. Many predate the current simulation framework.

These scripts are not part of the reproducible pipeline for the paper (nothing else in the repo sources or depends on them) — they're standalone, ad hoc checks kept for reference.

## Scripts

| File | Purpose |
|------|---------|
| `mini_simulation.R` | Basic mini simulation comparing estimators in a simple blocked CRT. |
| `mini_simulation_unblocked.R` | Mini simulation for completely randomized (unblocked) CRT designs. |
| `mini_simulation_block_variation.R` | Explores how block-level variation in ATEs affects estimator comparisons. |
| `mini_simulation_massive_variation.R` | Stress test with very high treatment effect heterogeneity. |
| `mini_simulation_massive_variation_noblock.R` | Same stress test but for unblocked designs. |
| `mini_simulation_MLM_explorations.R` | Compares multilevel model (MLM) estimators to simpler alternatives. |
| `superpop_SE_simulation.R` | Early version of the superpopulation SE simulation (now in `reports/`). |
| `simulation_DGP.R` | Prototype DGP code used before `R/data_generation_code.R` was developed. |
| `simulation_testing.R` | Ad hoc testing of simulation infrastructure. |
| `explore_person_bias.R` | Investigates bias of person-weighted vs. cluster-weighted estimators. |
| `person_v_cluster_weighting_demo.R` | Demonstration of when person-weighted and cluster-weighted estimands diverge. |
| `analyze_fake_data_check.R` | Checks estimator behavior on synthetic data with known ground truth. |

## Results Files

| File | Contents |
|------|---------|
| `pointwise_similarity_results.rds` | Saved output from a simulation examining which estimators give identical point estimates. |
