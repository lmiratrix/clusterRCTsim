
# Large simulation to look at stability of the different SE estimators and true
# performances of the ATE estimators.
#
# This script runs a bunch of scenarios and saves the massive amount of results
# to a file to be processed.


# clear workspace
rm(list=ls())


source( "simulation_support_code_crt.R" )

# Which set of simulation scenarios should we run?
#
# See the make_scenario_list function to figure out which table of
# simulations to pull from.
group = "sim_interceptcor"


library( tibble )
library( purrr )


#### Figure out simulation scenarios to run, and if we are testing or running ####

#SIZE_IMPACT_CORRELATE = FALSE
#VARIABLE_P = FALSE

if ( exists( "index" ) ) {
   scat( "Warning: Index already set!  Set to:\n" )
   print( index )
   rm( index )
}


# First try to get command line argument
args = commandArgs( trailingOnly = TRUE )
cat( "Command line arguments:\n" )
cat(args, sep = "\n")
if ( length( args ) > 0 ) {
  index = as.numeric( args[[1]] )
  scat( "Retrieved (and possibly overrode) index from command line: %d\n", index )
}

if ( !exists( "index" ) ) {
  # check if we're using a slurm array -- if so, use SLURM_ARRAY_TASK_ID
  if (Sys.getenv("SLURM_ARRAY_TASK_ID") != "") {
    index <- as.numeric(as.character(Sys.getenv("SLURM_ARRAY_TASK_ID")))
    scat( "Retrieved index from environment variable SLURM_ARRAY_TASK_ID: %d\n", index )
  } else {   # otherwise use INDEX_VAR, from execute_runfullsim.sh
    index <- as.numeric(as.character(Sys.getenv("INDEX_VAR")))
    scat( "Retrieved index from environment variable INDEX_VAR: %d\n", index )
  }
}


# Are we testing or no?  Set flags depending
if ( exists( "TESTING" ) || !exists( "index" ) || is.na( index ) ) {
  cat( "WARNING: In TEST mode for simulation script\n" )
  TESTING = TRUE
  index = 0 # c( 41, 42 )
  FILENAME = "results/testing_simulation_results_" 
  
  #  APPROX_NUM_TRIALS_PER_RUN = 20
  NUM_TRIALS_PER_SCENARIO = 2  
  
} else {
  cat( "Not testing, setting file names and approx num trials\n" )

  TESTING = FALSE
  FILENAME = "results/simulation_results_"

  NUM_TRIALS_PER_SCENARIO = 4
  #APPROX_NUM_TRIALS_PER_RUN = 700
  #NUM_TRIALS_PER_RUN = 4

}


# Set seed based on simulation index
set.seed( index^2 + index + 1 )



##### Start setup of simulation ####

make.file.name = function( index ) {
  if ( length( index ) > 1 ) {
    index = paste0( min(index), "_to_", max( index ) )
  }
  paste0( FILENAME, index, ".rds" )
}


##### Set up simulation factors and run simulation #####

scat( "Running simulation job %d\n", index )

ptm = proc.time()

# What simulation are we running?
#scenarios = get.scenario.by.id( index, group=group )
scenarios = make.scenario.list( group=group )

scenarios

if ( TESTING ) {
  cat( "Cutting scenario list down for testing\n" )
  scenarios = filter( scenarios, J.bar < 60 )
  scenarios = scenarios[1:5,]
}


#### Run the simulation ####

#NUM_TRIALS_PER_SCENARIO = round( APPROX_NUM_TRIALS_PER_RUN / nrow( scenarios ) )
scat( "Running %d trials per scenario\n", NUM_TRIALS_PER_SCENARIO )

scenarios$run = pmap( scenarios, run.scenario, R = NUM_TRIALS_PER_SCENARIO)

scat("**\n**\tTotal time elapsed:\n")
tot.time = proc.time() - ptm
print(tot.time)

sim.per.min = sim.per.min = nrow( scenarios ) * NUM_TRIALS_PER_SCENARIO / (tot.time["elapsed"] / 60)

scat( "Num trials / scenario = %d with %d scenarios = %d runs\n", NUM_TRIALS_PER_SCENARIO, nrow( scenarios ), NUM_TRIALS_PER_SCENARIO * nrow( scenarios ) )
scat("Realized simulations per minute = %.2f\n", sim.per.min )
# scenarios$sim.per.min = sim.per.min


scat( "**\n**")



#### Clean up the simulation results a wee bit ####

# Rearrange columns to be nice
scenarios = dplyr::select( scenarios, ID, everything() )
scenarios

# Add in (estimated) superpopulation person weighted ATE
scenarios = mutate( scenarios, 
                    ATE.super.student = map_dbl( run, ~ mean( .$ATE.finite.student ) ),
                    ATE.super.school = map_dbl( run, ~ mean(.$ATE.finite.school )),
                    ATE.super.dist.student = map_dbl( run, ~ mean(.$ATE.finite.dist.student ) ),
                    ATE.super.dist.school = map_dbl( run, ~ mean(.$ATE.finite.dist.school ) ))

# Rearrange and foreground relevant columns
scenarios = dplyr::select( scenarios, 
                           ID, K, J.bar, p.bar, sd.district, ICC.district, variable.J,
                           variable.p, n.bar, sd.school, ICC.school, variable.n,
                           size.impact.cor.district, size.impact.cor.school,
                           proptx.impact.cor.district,
                           
                           ATE.super.district = ATE.district, 
                           ATE.super.school, 
                           ATE.super.student, 
                           ATE.super.dist.student,
                           ATE.super.dist.school,
                           everything() )
scenarios

filename = make.file.name(index)
saveRDS( scenarios, file=filename )
scat( "Simulations saved to '%s'\n", filename )


# NOTE: ID refers to simulation scenario, NOT index variable!



