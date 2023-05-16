

##
## Process Odyssey results and make easier to manage summary results
##

# Collect all Odyssey results into a single file
source( "prepare simulation results/bundle_odyssey_results.R" )

# Clean up single file, remove bad sim runs, add in metadata and summary data
source( "prepare simulation results/clean_simulation_results.R" )

# Calculate performance metrics for the report
source( "prepare simulation results/prepare_aggregate_results.R" )

cat( "Finished meta aggregation script\n" )
