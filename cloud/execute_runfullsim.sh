
## Load software modules
module load R/4.2.2-fasrc01
## Use user-downloaded libraries
unset R_LIBS_SITE
export R_LIBS_USER=$HOME/apps/R-4.2.2:$R_LIBS_USER

for INDEX_VAR in $(seq 101 1000); do

  #print out indexes
  echo "${INDEX_VAR}"

  #give indexes to R so it can find them.
  export INDEX_VAR 

  #Run R script, and produce output files
  sbatch -o logs/sbout_p${INDEX_VAR}.stdout.txt \
	 -e out_p${INDEX_VAR1}_.stderr.txt \
         --job-name=runScr_p${INDEX_VAR} \
        sbatch_runScript.slurm
  
  sleep 1 # pause to be kind to the scheduler


done


