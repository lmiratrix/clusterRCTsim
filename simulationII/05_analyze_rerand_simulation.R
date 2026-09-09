

# Analyze the rerandomization simulations

# Read in all the rds files in simulationII/results_rerand_sim_public/

library( tidyverse )

rds_files = list.files( here::here( "simulationII/results_rerand_sim_public" ), full.names = TRUE )

data = map( rds_files, read_rds )
names(data) = basename(tools::file_path_sans_ext(rds_files))
data <- bind_rows( data, .id = "UniqueID" ) %>%
    as_tibble()

data


table( data$UniqueID, data$runID )
