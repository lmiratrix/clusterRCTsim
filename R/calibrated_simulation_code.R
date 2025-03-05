
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
                                        summarize_results = TRUE, ... ) {
  
  
  
  generate_and_analyze <- function( structure = NULL, parameters = NULL, ... ) {
    
    dat = generate_data_from_fingerprint( structure = structure,
                                          parameters = parameters )  
    
    compare_methods( Yobs ~ Z | clusterID | blockID, 
                     data = dat, ... )  
  }
  
  if ( parallel ) {
    library( furrr )
    plan(multisession, workers = parallel::detectCores() - 1 )
    rps = future_map( 1:R, \(.) { generate_and_analyze( structure, parameters, ... ) } )
  } else {
    rps = map( 1:R, \(.) { generate_and_analyze( structure, parameters, ... ) },
               .progress = TRUE )
  }
  
  if ( summarize_results ) {
    rps %>%
      bind_rows() %>%
      group_by( method, weight, biased ) %>%
      summarise( EATE = mean( ATE_hat ),
                 ESEhat = sqrt( mean( SE_hat^2 ) ),
                 sdSEhat = sd( SE_hat ),
                 Q10 = quantile( SE_hat, 0.1 ),
                 Q90 = quantile( SE_hat, 0.9 ),
                 SE = sd( ATE_hat ), .groups = "drop" )
  } else  {
    rps %>%
      bind_rows()
  }
  
}


# ****************************************************************************
# Rerandomization simulation code ----
# ****************************************************************************




#' Run a finite sample simulation by permuting treatment assignment
#' labels within each block.
run_rerandomize_simulation <- function( formula, data, R = 100, parallel = FALSE, ... ) {
  
  
  rerandomize_and_analyze <- function( formula, data, ... ) {
    
    form <- clusterRCT:::deconstruct_var_formula( formula, data = data )
    
    data[form$Z] = rerandomize( data[[form$Z]], data[[form$clusterID]], data[[form$blockID]] )
    
    compare_methods( formula, data, ... )  
  }
  
  if ( parallel ) {
    library( furrr )
    plan(multisession, workers = parallel::detectCores() - 1 )
    rps = future_map( 1:R, \(.) { rerandomize_and_analyze( formula, data, ... ) } )
  } else {
    rps = map( 1:R, \(.) { rerandomize_and_analyze( formula, data, ... ) },
               .progress = TRUE )
  }
  
  
  rps %>%
    bind_rows() %>%
    group_by( method, weight, biased ) %>%
    summarise( EATE = mean( ATE_hat ),
               ESEhat = sqrt( mean( SE_hat^2 ) ),
               sdSEhat = sd( SE_hat ),
               Q10 = quantile( SE_hat, 0.1 ),
               Q90 = quantile( SE_hat, 0.9 ),
               SE = sd( ATE_hat ), .groups = "drop" )
  
}
