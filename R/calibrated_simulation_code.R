
# Calibrated simulation code
#
# This file has code for running fingerprint and rerandomization
# simulations, which are simulations that target a specific dataset of
# set of estimates about a specific dataset.



library( tidyverse )
library( clusterRCT )
source( here::here( "R/util.R" ) )
source( here::here( "R/generate_data_from_fingerprint.R" ) )

# ****************************************************************************
# Fingerprint simulation code ----
# ****************************************************************************




#' Run a finite sample simulation by permuting treatment assignment
#' labels within each block.
run_fingerprint_simulation <- function( structure = NULL, parameters = NULL, 
                                        R = 100, parallel = FALSE,
                                        summarize_results = FALSE, ... ) {
  
  
  
  generate_and_analyze <- function( structure = NULL, parameters = NULL, ... ) {
    
    dat = generate_data_from_fingerprint( structure = structure,
                                          parameters = parameters )  
    
    compare_methods( Yobs ~ Z | clusterID | blockID, 
                     data = dat, ... )  
  }
  
  if ( parallel ) {
    library( furrr )
    plan(multisession, workers = parallel::detectCores() - 1 )
    rps = future_map( 1:R, \(.) { generate_and_analyze( structure, parameters, ... ) },
                      .options = furrr_options(seed = TRUE) )
    
    plan(sequential)
    gc()
    
  } else {
    rps = map( 1:R, \(.) { generate_and_analyze( structure, parameters, ... ) },
               .progress = TRUE )
  }
  
  summarize_simulation_results( rps, summarize_results = summarize_results )
  
}



