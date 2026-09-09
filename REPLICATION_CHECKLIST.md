# clusterRCTsim replication-readiness checklist

Two rounds of work so far, both on the git branch `cleanup/replication-prep` off `main` (**main is untouched**):

1.  An earlier session did a code review of `R/` alongside a parallel review of the companion `clusterRCT` package (see its own `CLEANUP_CHECKLIST.md`), and made two small fixes to `R/util.R`.
2.  This session (2026-09-03/04) read the full paper (main text + Appendices A-D) plus the earlier session's notes, then worked through every unambiguous cleanup item that didn't require your judgment call. Anything that did need your input was left alone and is listed under "Open questions" below.

## The single most important open item — resolved as a documented limitation

**Nothing in this repo records which exact `clusterRCT`/`blkvar` commit or version produced the numbers in the submitted JREE paper**, nor which exact scenario-grid configuration was active at the time (see open question #6, below). `install_needed_packages.R` installs both via `devtools::install_github(..., force = TRUE)` with no ref/tag pinning, there's no `renv.lock` or captured `sessionInfo()` anywhere, and `archive/` (which has its own README documenting the *scenario grid* snapshot convention) has no record of the *package version* that ran alongside the 2025-11-20 grid.

**Decision (2026-09-09): rather than keep chasing an answer we don't have, we're documenting this honestly instead.** Luke: there are too many moving parts across a multi-year project to reliably reconstruct exactly which version generated the final results, and that's OK — the goal of sharing this code is to let others understand the approach and run similar simulations, not to guarantee bit-for-bit reproduction of every reported number. A polite, upfront note to that effect is now in both `README.md` ("Reproducing the Paper" section) and `REPLICATION.md` (right after its intro). No further action needed here unless you later do track down the actual commit/config, in which case the note can be updated or removed.

## Done this session (uncommitted — review with `git diff` / `git status`)

**Reconciled the duplicate `create_analysis_tree` function** (the single biggest landmine found): there were two files each defining a same-named, non-interchangeable `rpart`-based tree-fitting function — `create_analysis_tree.R` at the repo root (`min_leaves`/`max_leaves` API, used by `fingerprint/fingerprint_simulation_report.Rmd`) and `analysis/R/create_analysis_tree.R` (`num_leaves` + legacy column renaming + signed color scale, used by `analysis/Simulation_WriteUp.Rmd` and `empirical/empirical_results_report.Rmd`). Renamed the root one to `create_analysis_tree_minmax()` in `analysis/R/create_analysis_tree_minmax.R`, updated all 3 call sites in `fingerprint_simulation_report.Rmd`, added cross-reference header comments to both files so nobody "deduplicates" them again by accident, and documented the split in a new `analysis/R/README.md`. **Behavior of every existing report is unchanged** — this was a pure rename/relocation, not a logic change.

**Fixed two broken, pre-existing source() paths** left over from an old `analysis/Simulation_Write_Up/` → `analysis/` folder flattening that never got fully propagated: `explore_specific_theoretical_predictions.R` (both its `load_and_prep_data.R` and `create_analysis_tree.R` references) and `analysis/make_tree_of_sim_results.R`'s `load_and_prep_data.R` reference. Both scripts were silently non-functional before this fix.

**Relocated 4 loose top-level scripts** into themed, READMEd folders (all `git mv`, history preserved): `create_analysis_tree.R` → `analysis/R/create_analysis_tree_minmax.R` (above), `make_tree_of_sim_results.R` → `analysis/`, `compare_methods_across_data_structures.R` and `explore_specific_theoretical_predictions.R` → `explorations/`. Updated `analysis/README.md` and `explorations/README.md` accordingly.

**Fixed a real bug in `cloud/execute_runfullsim.sh`**: the stderr log path used an undefined `INDEX_VAR1` instead of `INDEX_VAR`, so every array job's stderr was silently written to the same overwritten file (`out_p_.stderr.txt`) instead of one per job. Now `logs/sbout_p${INDEX_VAR}.stderr.txt`, matching the stdout pattern.

**Labeled all Harvard-FASRC-specific HPC scripts** (`cloud/*.slurm`, `execute_runfullsim.sh`) with a one-line header noting they're `miratrix_lab`-account-specific and need adaptation for another institution.

**Deleted 3 untracked Dropbox "conflicted copy" junk files** in `fingerprint/attic/` (confirmed untracked/gitignored, pure sync duplicates) and removed `data_prep/README.md`'s now-stale warning about conflicted copies (the files it named no longer exist there).

**Centralized the hardcoded secure-server path** (`/data/share/pse/IES CRT/...`, previously a literal string repeated across 7 files) into a single `IES_CRT_ROOT <- Sys.getenv("IES_CRT_DATA_ROOT", "/data/share/pse/IES CRT")` constant in each of the 6 live files that use it (`fingerprint/master_run_fingerprints.R`, `empirical/{check_quality_control,beaconmentor_test,master_analyze_all_data}.R`; left `empirical/load_empirical_results.R`'s copy alone since it was dead code, see below). Documented the override variable in `empirical/README.md`. Left `empirical/attic/nonblocked_experiment_analysis.R` (historical, unused) and a commented-out example line in `empirical/excel_supp.R` untouched.

**Found and removed genuinely dead code while doing the above**: `empirical/load_empirical_results.R` set `PLOT_PATH` to the hardcoded secure path on one line, then unconditionally overwrote it with a repo-relative `here::here(...)` path on the very next line — the hardcoded value was never actually used. Deleted the dead line.

**Clarified where the empirical report's output actually goes**: `empirical/master_analyze_all_data.R` was writing its final `master241204_combined.RDS` only to the secure-server mirror of the repo (`/data/share/pse/IES CRT/clusterRCTsim/empirical/results/...`). Changed that specific line to `here::here("empirical/results/...")` so it's explicit that this is the file meant to land back in the actual git repo (this is the same file that's already committed and that `empirical_results_report.Rmd` reads from — see the reproducibility note below).

**Documented, in `empirical/README.md`, that the empirical report is reproducible without secure-server access**: traced `empirical_results_report.Rmd` → `load_and_prep_empirical_results.R` and confirmed it reads only from the committed, aggregated `empirical/results/*.RDS` files, not raw data. Only regenerating those aggregate files requires server access. This is a real strength worth a reviewer knowing about.

**Fixed a `here::here()` inconsistency**: `simulation/estimator_performance_simulation_fixed_district.R` (the main cloud/grid entry point) used plain relative `source("R/...")` instead of `here::here(...)`, contradicting the repo's own stated convention — fragile if sourced from any working directory other than repo root. Fixed.

**Documented the seeding convention** in a new section of `R/README.md`: global `set.seed(40444)` for grid generation, per-scenario `set.seed(index^2 + index + 1)` for cloud/array simulation runs, and confirmed the remaining fixed-constant seeds in `util.R`/`group_estimators.R`/`generate_nonnegative_integers.R` are all inside `if (FALSE) {...}` demo blocks and don't affect any real run.

**Wrote `archive/README.md`** documenting the dated scenario-grid-snapshot convention (previously only described in the root README), with a placeholder noting the still-missing package-version record (see top item above).

**Renamed** `analysis/attic/Archived/Simulation_WriteUp copy.Rmd` → `Simulation_WriteUp_copy.Rmd` (removed the space; unreferenced anywhere, confirmed before renaming). Note: several *other* files in that same `Archived/` folder also have spaces in their names (`Range Data Results.Rmd`, `Simulation Memo 2.Rmd`, etc.) — left untouched since they were outside the specific item flagged, but same fix would apply if you want it.

**Archived superseded covariate-simulation versions**: moved `covariate_exploration_simulation_v1.R`/`_v2.R` into `covariate_simulations/attic/` now that `_v3.R` is confirmed current; updated that folder's README, and added a scope note there (see open question below).

**Added "not part of the reproducible pipeline" disclaimers** to `explorations/README.md` and `small_simulations/README.md` (nothing else in the repo references either folder's contents).

**Added `CITATION.cff`** at the repo root, using the paper's actual title and author list (Miratrix, Weiss, Litschwartz, Hill, Warner) pulled directly from the submitted PDF's title page. Marked as `unpublished`/under review with a note to update once accepted (journal, volume, DOI, year).

**In progress (background):** a fork is building `REPLICATION.md`, a figure/table → producing-script map covering the main text and Appendices A-C (Appendix D is pure math, no code component). This is large and will land as a separate addition once done — check for it if it's not here yet.

## Open questions for you

These all need a judgment call only you can make — nothing was changed for any of these:

1.  ~~Which `clusterRCT`/`blkvar` commit produced the submitted results?~~ **Resolved as a documented limitation** — see the top of this file. Not pinned down, and we're saying so in `README.md`/`REPLICATION.md` rather than guessing.

2.  ~~`logs/` (609MB, \~2,000 files) is fully committed to git~~ **Resolved (2026-09-10).** Gitignored on `main`/`dev`/`cleanup/replication-prep` (the only 3 branches that ever had it), then purged from git history entirely via `git filter-repo` across all branches and force-pushed. `.git` went from 81MB → 54MB. Verified via a fresh clone from GitHub that no trace remains. As with the earlier Beacon Mentor purge, anyone with an existing clone needs to re-sync (re-clone, or `git fetch` + `git reset --hard origin/<branch>` per branch) — let your collaborators know if they have local clones of any of the 3 affected branches.

3.  ~~`util.R` reimplements `rerandomize()` locally~~ **Resolved (2026-09-10).** Confirmed byte-for-byte identical to `clusterRCT::rerandomize()`, then deleted the local copy. Callers are qualified as `clusterRCT::rerandomize()` rather than left bare — turned out most simulation entry points in this repo *don't* call `library(clusterRCT)` (they rely on `R/util.R`'s `get.estimates()` using qualified `clusterRCT::compare_methods()` instead), so a bare call would have silently broken `rerandomize.data()` in those scripts. Also qualified the one other bare call site, in `R/generate_data_from_fingerprint.R`, which had the added wrinkle of a same-named boolean parameter (`rerandomize = TRUE`) shadowing the function name for ordinary lookups. Verified with a functional test (no `library(clusterRCT)` loaded, matching the real entry points) that `rerandomize.data()` still works. Applied identically on `main` and `dev`.

4.  ~~Does `covariate_simulations/` belong to this paper at all?~~ **Resolved — moot.** The folder was removed from `main` (kept only on `dev`) as part of trimming `main` to a replication package; the scope question no longer matters for `main`, and `dev`'s copy of the folder's README still carries the original scope note for whoever picks it up later.

5.  ~~Several committed report PDFs look force-added~~ **Resolved (2026-09-10).** Removed `fingerprint/rerandomization_simulation_report.pdf`, `empirical/empirical_results_report.pdf`, and `analysis/attic/designing_simulation2.pdf` from git on both `main` and `dev` — they regenerate from source `.Rmd` (the first two) or are non-load-bearing historical scratch material (the third, already living in `attic/`).

6.  ~~Which scenario `group` produced the submitted results?~~ **Resolved as a documented limitation** — see the top of this file. `simulation/estimator_performance_simulation_fixed_district.R` still selects the active scenario grid by manually commenting/uncommenting lines (`group = "full_factor"` vs. several alternatives) rather than a parameter, and we're not attempting to reconstruct which one matches the paper after the fact.

7.  ~~`covariate_simulations/covariate_exploration_simulation.R` is untracked~~ **Resolved — moot on `main`** (folder no longer there), **and fixed on `dev`** (the file is now tracked/committed there, so it's no longer at risk of being silently lost).

## Not done — correctness findings needing your judgment (from the earlier session)

From that session's review of `R/data_generation_code.R`, `R/simulation_support_code_crt.R`, `R/simulation_support_code_fixed_block.R`, `R/calibrated_simulation_code.R`, `R/generate_data_from_fingerprint.R`, `R/aggregate_simulation_results.R`, `R/group_estimators.R`, `R/MCSE_functions.R`:

1.  `generate_data_from_fingerprint.R` uses `blkvar:::make_covariate_set` — an unexported internal of a dependency, via triple-colon. Fragile to `blkvar` internal changes, no test guards it. **Decision (2026-09-10): leave as is** — not ideal, but livable.
2.  `single.CRT()`'s own comment (`simulation_support_code_crt.R`, \~line 157) acknowledges design-based estimators can return infinite SEs with fewer than two treated/control clusters in a block; the size-generation logic mitigates but doesn't provably eliminate this across a large random grid sweep. **Decision (2026-09-10): ignore.**
3.  No bugs found in `MCSE_functions.R`, `group_estimators.R`, `generate_nonnegative_integers.R`, the core DGP math in `data_generation_code.R`, or `aggregate_simulation_results.R`. `single.CRT()`/`single.CRT.fixed.block()` computing ground-truth finite-sample ATEs directly from `Y1`/`Y0` before rerandomization is a sound design — a good sign for the simulation's core correctness.

## Good practices already in place (worth keeping)

- `results/` is properly gitignored (raw per-job cluster output, regenerable); `results_agg/` is intentionally tracked (aggregated outputs, so reports are reproducible without rerunning the full simulation) — exactly right.
- `attic/`, `luke_scratch/`, `explorations/`, `small_simulations/`, `demo_scripts/`, `documentation/` each have their own README self-describing as non-load-bearing scratch/historical/demo material.
- The `archive/` dated-snapshot convention (copy `generate_grids.R` + `scenarios.rds` before regenerating) is exactly the right instinct — see open question #1 for the one piece still missing.
- No stray `.Rhistory`/`.DS_Store`/`.RData` tracked outside the intentional `results_agg/` files.
- The empirical report is reproducible from committed aggregate data without secure-server access (confirmed and now documented this session).

## Test-mining: ready to port into `clusterRCT`'s test suite

(Already actioned the first one — see `clusterRCT`'s `tests/testthat/test-rerandomize.R`.)

- `tests/testthat/test-generate_data_from_fingerprint.R` (lines \~22-46) generates a fingerprint-matched synthetic dataset from `fakeCRT`'s structure and asserts `describe_clusterRCT()` reproduces the same structural stats on both real and synthetic data. Since `fakeCRT` ships in the `clusterRCT` package itself, the `describe_clusterRCT()`-on-`fakeCRT`/`fakeCRT2`/`fakeBrokeCRT` half of this is portable as a value-level regression test for `clusterRCT`.
- `tests/testthat/test-DGP.R`'s "overall variation is 1" / "capture ICC and R2" pattern (generate data with a *known* ICC/R², assert `lmer()` recovers it within tolerance) is a good template to copy for a `clusterRCT` regression test of `describe_clusterRCT()`'s unblocked-branch ICC calculation.
