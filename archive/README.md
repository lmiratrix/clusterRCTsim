# archive/ — Dated Scenario-Grid Snapshots

`simulation/generate_grids.R` overwrites `scenarios.rds` in place every time it's rerun, so before regenerating the grid we copy the current `generate_grids.R` + `scenarios.rds` pair into this folder, named with the date the grid was created. This preserves exactly the code and grid that produced any past set of simulation results, even after the live grid-generation code moves on.

## Snapshots

| Date | Files | Notes |
|------|-------|-------|
| 2025-11-20 | `generate_grids_2025-11-20.R`, `scenarios_2025-11-20.rds` | The `generate_grids.R` version and resulting grid (three block types: same, reasonable, messy; 1,296 scenarios) used for the initial JREE journal submission. |

## Open item: recording the estimator-package version alongside each grid

This folder records *which scenario grid* produced a given set of results, but not *which version of the `clusterRCT`/`blkvar` packages* was installed when the corresponding cloud simulation ran. Since those packages implement the estimators themselves and are under active development, the grid alone isn't enough to exactly reproduce a past result — see the open questions in `REPLICATION_CHECKLIST.md` for pinning a specific commit/tag to each snapshot.
