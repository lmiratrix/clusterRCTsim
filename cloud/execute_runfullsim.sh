
## MAIN SIMULATION LAUNCHING SCRIPT
## How Sophie runs this file:
## Copy and paste this code into the terminal directly, if desired.
##
## Institution-specific: written for the Harvard FASRC (Cannon) cluster under
## the miratrix_lab account. Adapt the module name, R_LIBS_USER path, and
## sbatch account/partition settings (see cloud/sbatch_runScript.slurm) for
## another institution's HPC environment.


## Load software modules
module load R/4.2.2-fasrc01
## Use user-downloaded libraries
unset R_LIBS_SITE
export R_LIBS_USER=$HOME/apps/R-4.2.2:$R_LIBS_USER

for INDEX_VAR in $(seq 2 1000); do

  #print out indexes
  echo "${INDEX_VAR}"

  #give indexes to R so it can find them.
  export INDEX_VAR 

  #Run R script, and produce output files
  sbatch -o logs/sbout_p${INDEX_VAR}.stdout.txt \
	 -e logs/sbout_p${INDEX_VAR}.stderr.txt \
         --job-name=runScr_p${INDEX_VAR} \
        cloud/sbatch_runScript.slurm
  
  sleep 1 # pause to be kind to the scheduler


done


