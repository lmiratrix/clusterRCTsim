# covariate_simulations/ — Covariate Adjustment Simulations

This folder explores the behavior of estimators that adjust for covariates in blocked cluster-randomized trials, particularly in the low-df regime where there are few clusters per block.

## Report

| File | Output | Description |
|------|--------|-------------|
| `on_covariate_adjustment.rmd` | `on_covariate_adjustment.pdf` | **Main report.** Investigates how individual-level (L1) and cluster-level (L2) covariates interact with fully-interacted estimators, and when covariate adjustment helps versus hurts. Addresses collinearity problems when adding cluster-level covariates to a fully-interacted model. |

## Scripts

| File | Purpose |
|------|---------|
| `covariate_exploration_simulation.R` | Initial simulation exploring how covariates affect estimator performance. |
| `covariate_exploration_simulation_v3.R` | Third iteration; most current version. |
| `DGP_code_varying_cluster_sizes.R` | Data-generating process code with variable cluster sizes, used for the covariate sims. |
| `matched_pairs_simulation.R` | Simulation focused on matched-pairs designs (two clusters per block) with covariate adjustment. |
| `sanity_simulation.R` | Sanity-check simulation verifying basic properties of the covariate-adjusted estimators. |
| `sim_analysis_function.R` | Analysis functions shared across the simulation scripts in this folder. |

`covariate_exploration_simulation_v1.R` and `_v2.R` are superseded by `_v3.R` and have been moved to `attic/`.

## Scope note

Nothing in the JREE paper's main text or appendices (Estimands Main Text, Appendices A–D) references covariate-adjustment simulations, low-df matched-pairs simulations, or this folder by name — confirmed by a close read of all four appendices. This folder may belong to a separate manuscript/exploration rather than this paper's replication package; worth confirming before including/excluding it from any paper-facing replication documentation (see the open question in `REPLICATION_CHECKLIST.md`).
