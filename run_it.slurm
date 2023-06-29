#!/bin/bash

#SBATCH -A miratrix_lab # lab account

#SBATCH -n 1  #Number of cores

#SBATCH -t 20  #Runtime in minutes

#SBATCH -p serial_requeue #Partition to submit to

#SBATCH --mem=1000 #Memory per node in MB (see also --mem-per-cpu)

#SBATCH -o logs/sim%a.out # Standard output
#SBATCH -e logs/sim%a.err # Standard error
#SBATCH -a 4-100

## Load software modules
module load R/4.2.2-fasrc01
## Use user-downloaded libraries
unset R_LIBS_SITE
export R_LIBS_USER=$HOME/apps/R-4.2.2:$R_LIBS_USER

#R file to run, and txt files to produce for output and errors
R -f estimator_performance_simulation.R

