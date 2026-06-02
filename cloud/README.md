# cloud/ — Harvard FASRC Cluster Scripts

This folder contains SLURM batch scripts for running the large-scale simulation on the Harvard FASRC (Cannon/Odyssey) cluster under the `miratrix_lab` account.

## Files

| File | Purpose |
|------|---------|
| `run_it.slurm` | Main job array script. Submits a SLURM array job (`-a 3-100`) where each task runs one chunk of the simulation. Requests 1 core, 30 min runtime, 1 GB memory on the `serial_requeue` partition. |
| `sbatch_runScript.slurm` | Wrapper script for submitting `run_it.slurm` via `sbatch`. |
| `results.slurm` | Script for post-processing: bundles and cleans simulation results on the cluster after the array job finishes. |
| `test_it.slurm` | Test script for a single trial run to verify the setup before submitting a full job array. |
| `execute_runfullsim.sh` | Shell script that orchestrates the full simulation workflow: submit jobs, wait, then run results processing. |

## Workflow

1. Push current code and `scenarios.rds` to the cluster.
2. Run `execute_runfullsim.sh` (or submit `run_it.slurm` manually with `sbatch`).
3. After jobs complete, run `results.slurm` to bundle outputs.
4. Pull the bundled results back locally and run the `data_prep/` pipeline.
