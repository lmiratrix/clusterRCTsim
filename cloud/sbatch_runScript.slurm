#!/bin/bash

#SBATCH -A miratrix_lab # lab account

#SBATCH -n 1  #Number of cores

#SBATCH -t 180  #Runtime in minutes

#SBATCH -p serial_requeue #Partition to submit to

#SBATCH --mem=5000 #Memory per node in MB (see also --mem-per-cpu)

#R file to run, and txt files to produce for output and errors
R CMD BATCH estimator_performamce_simulation_fixed_district.R logs/R_output_${INDEX_VAR}.txt --no-save




