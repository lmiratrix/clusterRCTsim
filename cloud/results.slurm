#!/bin/bash

# NOTE: This is old and not used.
#SBATCH -A miratrix_lab # lab account
#SBATCH -n 1  #Number of cores
#SBATCH -t 60  #Runtime in minutes
#SBATCH -p test #Partition to submit to
#SBATCH --mem=4000 #Memory per node in MB (see also --mem-per-cpu)

## Load software modules
module load R/4.2.2-fasrc01
## Use user-downloaded libraries
unset R_LIBS_SITE
export R_LIBS_USER=$HOME/apps/R-4.2.2:$R_LIBS_USER

#R file to run, and txt files to produce for output and errors
R -f fully_prepare_simulation_results.R

