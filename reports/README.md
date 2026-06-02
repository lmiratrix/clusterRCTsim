# reports/ — Technical Background Reports

This folder contains standalone technical reports that inform the methodology used in the project but are not the primary deliverables. They explore theoretical properties of the standard errors and provide an orientation tutorial.

## Reports

| File | Output | Description |
|------|--------|-------------|
| `clusterRCT_SEs.Rmd` | `clusterRCT_SEs.pdf` | Examines when superpopulation and finite-population standard errors diverge substantially. Uses math derivations (under a simple model with site-mean/impact correlations) validated by simulation. Key finding: large divergence requires large treatment effect heterogeneity, but in that regime the finite-sample SEs are also large. |
| `superpop_SE_report.Rmd` | `superpop_SE_report.pdf` | Simulation study where block characteristics are fixed and clusters are repeatedly resampled. Compares how much point estimates vary across simulation runs (i.e., the superpopulation SE) for each estimator. Uses nested simulations: first generate a fixed district structure, then repeatedly randomize and estimate. |
| `tutorial.Rmd` | `tutorial.pdf` | Getting-started guide for running clusterRCTsim simulations, both locally and on the Harvard FASRC (Odyssey/Cannon) cluster. Written by Jonathan Che. Describes the file structure and workflow. |
