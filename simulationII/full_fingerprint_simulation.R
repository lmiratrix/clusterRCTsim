

# This script runs a simulation that, for each row of study
# characteristics, repeatedly generates data of that form and then
# analyzes it with compare_methods.
#
# Each study is saved as a separate result file.
#
# This is NOT the rerandomization simulation, but instead a simulation
# where the datasets are generated from whole cloth, using the
# parameters estimated from the original empirical data.


library( tidyverse )
library( furrr )
library( cli )
library( glue )
source( here::here( "R/util.R" ) )
source( here::here( "R/calibrated_simulation_code.R" ) )
library( digest )


# create finger_result folder
RESULT_DIR = here::here( "simulationII/results_fingerprint" )
dir.create( RESULT_DIR, showWarnings = FALSE )


# List of studies found in the rerandomization effort -- this is so
# this simulation can be subset so it is as similar as possible to the
# rerandomization simulation as we were finding odd behavior in the
# rerand simulation.

rerandomization_list <- c("Beacon_S8TGPA", "Beacon_S8TOTCR", "Beacon_SRPRST2", "EMERGE_team2__BB", "EMERGE_eclsb2_BB", "EMERGE_team2__BBSSR", "EMERGE_eclsb2_BBSSR", "MPC_eclsb2_scale", "MPC_wj_ap2", "PowerTeaching_TSTSSTAND", "Success_S14WJLWI_SCORE", "TennStar_gktmathss", "TennStar_gktreadss", "TRIAD_rema2", "TRIADScale-up_rema2", "jg_crt_2", "jg_crt_7", "jg_crt_8", "jg_crt_10", "jg_crt_11", "jg_crt_12", "jg_crt_16", "jg_crt_17", "jg_crt_18", "jg_crt_19", "jg_crt_23", "jg_crt_29", "jg_crt_30", "jg_crt_31", "jg_crt_34", "jg_crt_41", "jg_crt_56", "jg_crt_63", "jg_crt_68", "jg_crt_69")

studies = read_csv( here::here( "empirical/results/study_descriptions.csv" ),
                    show_col_types = FALSE )

studies <- filter( studies,
                   blocked == "Blocked",
                   #                   UniqueID %in% rerandomization_list,
                   n < 25000 )

studies


# Drop studies already done.

  # Get all files in results_fingerprint directory that are rds files
  files = list.files( here::here( "simulationII/results_fingerprint" ),
                      pattern = "\\.rds$" )
  files = gsub( "\\.rds$", "", files )
  cli_inform( "Found {length(files)} existing result files" )
  
  studies$have_results = studies$UniqueID %in% files
  


cli_inform( "Working with {nrow(studies)} studies" )




# Code to run single fingerprint for single study ----


string_to_seed <- function(str) {
  ss <- digest(str, algo = "md5")
  strtoi(substr(ss, 2, 7), base = 16)
}



# install.packages("filelock"); install.packages("cli")
log_file <- here::here( "simulationII/all_output.log" )

log_line <- function(text, .envir = parent.frame()) {
  #lk <- filelock::lock(lock_file, timeout = 0)  # returns NULL if already locked
  #on.exit(try(filelock::unlock(lk), silent = TRUE), add = TRUE)
  
  
  #on.exit(try(close(con), silent = TRUE), add = TRUE)
  msg <- glue(text, .envir = .envir)
  
  if ( FALSE ) {
    con <- file(log_file, open = "a", blocking = TRUE)
    cat(sprintf("%s [pid %s] %s\n", 
                format(Sys.time(), "%F %T"),
                Sys.getpid(), 
                msg),
        file = con, append = TRUE)
    flush(con)
  } else {
    cat( msg )
  }
  
}




run_one_finger <- function( k, R = 10, ... ) {
  
  # For parallel processing to function:
  # Load these things to ensure they are global variables in the
  # future worlds
  rerandomize
  
  # Grab relevant row of study
  stud = studies[k,] 
  
  log_line( "Running fingerprint for {k}: {stud$UniqueID} R={R}" )
  log_line( "  n={stud$n} J={stud$J} K={stud$K}" )
  
  # set seed
  set.seed( string_to_seed( stud$UniqueID ) )
  
  # Force covariates to 1 to get "best"
  stud$ncov.1 = 1
  stud$ncov.2 = 1
  
  # Standardize outcome
  stud$sdY0 = 1
  
  generate_and_analyze <- function( parameters = NULL, ... ) {
    dat = generate_data_from_fingerprint( structure = NULL,
                                          parameters = stud )  
    
    ctlX <- compare_methods( Yobs ~ Z | clusterID | blockID, 
                             data = dat,
                             control_formula = ~ X,
                             include_method_characteristics = FALSE,
                             ... )  
    noCtl <- compare_methods( Yobs ~ Z | clusterID | blockID, 
                              data = dat,
                              control_formula = NULL,
                              include_method_characteristics = FALSE,
                              ... )
    bind_rows( X = ctlX, none = noCtl, .id = "control" )
  }
  
  res = map( 1:R, \(.) { generate_and_analyze( stud, ... ) },
             .progress = TRUE )
  
  
  res = bind_rows( res, .id = "runID" )
  
  # save results to file
  fname = here::here( RESULT_DIR, paste0( stud$UniqueID, ".rds" ) )
  saveRDS( as_tibble( res ), file = fname )
  
  log_line( "  Finished process for {k}: {stud$UniqueID}" )
  
  # Success!
  fname
}

if ( FALSE ) {
  debug( run_fingerprint_simulation )
  debug( run_one_finger )
  
  f <- run_one_finger( 5 )
  
  
  res <- read_rds( f )
  res
  table( res$method )
  table( res$weight )
}


# Run the simulation ----

# sort studies from large to small to roughly balance the cluster
# usage

# studies <- arrange( studies, desc(sqrt(n)*J) )

if ( TRUE ) {
  
  studies = filter( studies, !have_results )
  cli_inform( "After dropping existing files, {nrow(studies)} studies remain" )
}



studies %>%
  dplyr::select( UniqueID, have_results, n, nbar, J, Jbar, K ) %>%
  knitr::kable( digits = 1 )

if ( FALSE ) {
  # testing
  run_one_finger( 3, R=200 )
}

if ( FALSE ) {
  cli_h3("Launching simulation to tackle {nrow(studies)} studies in parallel")
  
  library( purrr )
  library( future )
  
  plan(multisession, workers = parallel::detectCores() - 3 )
  
  
  tictoc::tic()
  future_walk( 1:nrow(studies), run_one_finger, 
               R = 200,
               .options = furrr_options( seed = NULL ),
               .progress = TRUE )
  tictoc::toc()
  
  # close parallel
  cli_h3("Closing multiple plans")
  plan(sequential)
  
} else {
  cli_h3("Launching simulation in sequential, not parallel, mode" )
  
  walk( 1:nrow(studies), run_one_finger, R = 200,
        .progress = TRUE )
}

cli_h1( "Script complete" )
