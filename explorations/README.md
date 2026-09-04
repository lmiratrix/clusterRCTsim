# explorations/ — Ad Hoc Exploratory Analyses

Scripts for targeted explorations that do not fit into the main simulation or empirical folders. These typically address specific theoretical questions or referee/collaborator queries.

## Files

| File | Purpose |
|------|---------|
| `checking_cv_vs_size_ratio.R` | Examines the relationship between the coefficient of variation (CV) of cluster sizes and the ratio of standard errors across estimators. |
| `mike_weighing_versions_exploration.R` | Explores different weighting versions of the estimators following suggestions from Mike Weiss. |
| `explore_specific_theoretical_predictions.R` | Checks simulation output against specific theoretical predictions about estimator bias/behavior. Relocated here from the repo root. |
| `compare_methods_across_data_structures.R` | Runs all estimators on a handful of hand-constructed data structures (e.g. matched pairs) to sanity-check behavior on edge-case designs. Relocated here from the repo root. |

These scripts are not part of the reproducible pipeline for the paper (nothing else in the repo sources or depends on them) — they're standalone, ad hoc checks kept for reference.
