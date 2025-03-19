
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



process_results <- function( rps, summarize_results = TRUE ) {
  
  if ( summarize_results == "cross" || summarize_results == "cross-agg" ) {
    # Calculate range across estimators
    rps <- rps %>%
      bind_rows( .id = "runID" )
    
    rr_overall <- rps %>%
      group_by( runID ) %>%
      summarise( rangeATE = max(ATE_hat, na.rm=TRUE) - min(ATE_hat, na.rm=TRUE),
                 na_ATE = sum( is.na( ATE_hat ) ),
                 ratioSE = max( SE_hat, na.rm=TRUE ) / min( SE_hat, na.rm=TRUE ),
                 na_SE = sum( is.na( SE_hat ) ) )
    
    rr_estimand <- rps %>%
      mutate( weight = fct_recode( weight,
                                   C = "Cluster",
                                   P = "Person" ) ) %>%
      group_by( runID, weight ) %>%
      summarise( rangeATE = max(ATE_hat, na.rm=TRUE) - min(ATE_hat, na.rm=TRUE),
                 na_ATE = sum( is.na( ATE_hat ) ),
                 ratioSE = max( SE_hat, na.rm=TRUE ) / min( SE_hat, na.rm=TRUE ),
                 na_SE = sum( is.na( SE_hat ) ),
                 .groups ="drop" )
    
    rr_overall$weight = "overall"
    
    rr = bind_rows( rr_overall, rr_estimand )
    
    if ( summarize_results == "cross-agg" ) {
      ragg <- rr %>%
        group_by( weight ) %>%
        rename( group = weight ) %>%
        summarise( mn_rngATE = mean( rangeATE ),
                   sd_rngATE = sd( rangeATE ),
                   prob_10p = mean( rangeATE > 0.10 ),
                   mn_ratSE = mean( ratioSE ),
                   prob_2x = mean( ratioSE >= 2 ),
                   sd_ratSE = sd( ratioSE ),
                   mn_naATE = mean( na_ATE ),
                   sd_naATE = sd( na_ATE ),
                   mn_naSE = mean( na_SE ),
                   sd_naSE = sd( na_SE ),
                   .groups = "drop" )
      
      ragg
    } else {
      rr
    }
    
  } else if ( summarize_results ) {
    # Calculate performance statistics for each method
    rps %>%
      bind_rows() %>%
      group_by( method, weight, biased ) %>%
      summarise( EATE = mean( ATE_hat, na.rm=TRUE ),
                 SE = sd( ATE_hat, na.rm=TRUE ), 
                 p.noATE = mean( is.na( ATE_hat ) ),
                 p.noSE = mean( is.na( SE_hat ) ),
                 ESEhat = sqrt( mean( SE_hat^2, na.rm=TRUE ) ),
                 sdSEhat = sd( SE_hat, na.rm=TRUE ),
                 Q10 = quantile( SE_hat, 0.1, na.rm=TRUE ),
                 Q90 = quantile( SE_hat, 0.9, na.rm=TRUE ),
                 .groups = "drop" )
  } else {
    rps %>%
      bind_rows( .id = "runID" )
  }
  
}


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
    rps = future_map( 1:R, \(.) { generate_and_analyze( structure, parameters, ... ) },
                      .options = furrr_options(seed = TRUE) )
  } else {
    rps = map( 1:R, \(.) { generate_and_analyze( structure, parameters, ... ) },
               .progress = TRUE )
  }
  
  process_results( rps, summarize_results = summarize_results )
  
}


# ****************************************************************************
# Rerandomization simulation code ----
# ****************************************************************************




#' Run a finite sample simulation by permuting treatment assignment
#' labels within each block.
run_rerandomize_simulation <- function( formula,
                                        data, 
                                        R = 100, 
                                        summarize_results = TRUE,
                                        parallel = FALSE, ... ) {
  
  
  rerandomize_and_analyze <- function( formula, data, 
                                       ... ) {
    
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
  
  process_results( rps, summarize_results = summarize_results )
  
}
