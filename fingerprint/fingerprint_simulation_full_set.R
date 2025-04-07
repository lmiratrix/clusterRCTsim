
# This simulation generates a single dataset for each fingerprint in
# our file, and then calculates the range and whatnot of the
# estimates.

library( tidyverse )
library( furrr )
library( cli )

source( "R/util.R" )
source( "R/calibrated_simulation_code.R" )


# create finger_result folder
dir.create( here::here( "fingerprint_results" ), showWarnings = FALSE )



studies = read_csv( here::here( "empirical/results/study_descriptions.csv" ),
                    show_col_types = FALSE )
studies <- filter( studies,
                   blocked == "Blocked",
                   n < 30000 )
studies

cli_inform( "Working with {nrow(studies)} studies" )

# Code to run single fingerprint for single study

library( digest )
string_to_seed <- function(str) {
  ss <- digest(str, algo = "md5")
  strtoi(substr(ss, 2, 7), base = 16)
}



run_one_finger <- function( k, R = 10 ) {
  
  # Load these things to ensure they are global variables in the
  # future worlds
  rerandomize
  
  # Grab relevant row of study
  stud = studies[k,] 
  
  cat( "Running fingerprint for ", k, ": ", stud$UniqueID, "\n" )
  
  # set seed
  set.seed( string_to_seed( stud$UniqueID ) )
  
  res <- run_fingerprint_simulation( structure = NULL,
                                     parameters = stud, 
                                     R = R,
                                     include_MLM = FALSE,
                                     summarize_results = "cross",
                                     parallel = FALSE )
  
  # save results to file
  saveRDS( res, file = here::here( "fingerprint_results", paste0( stud$UniqueID, ".rds" ) ) )
  
  # Success!
  1
}

if ( FALSE ) {
  run_one_finger( 5 )
}


if ( TRUE ) {
  # run all the fingerprints (in parallel)
  
  library( purrr )
  library( future )
  
  plan(multisession, workers = parallel::detectCores() - 1 )
  
  future_walk( 1:nrow(studies), run_one_finger, R = 20,
                     .options = furrr_options( seed = NULL ),
                     .progress = TRUE )
  
} else {
  
  walk( 1:nrow(studies), run_one_finger, R = 20,
              # .options = furrr_options( seed = NULL ),
              .progress = TRUE )
}


# Load and examine data ----

if ( FALSE ) {
  
  # Read one example result file
  res = read_rds( here::here( "fingerprint_results", paste0( studies$UniqueID[1], ".rds" ) ) )
  res
  
}

dat <- studies$UniqueID %>%
  set_names() %>%
  map( \(file) {
    fname = here::here( "fingerprint_results", paste0( file, ".rds" ) )
    # check if file exists
    if ( file.exists( fname ) ) {
      read_rds( fname )
    } else {
      NULL
    }
  } )  %>% 
  bind_rows( .id = "UniqueID" )

dim( dat )
dat
table( dat$runID)

# Who is missing?
setdiff( studies$UniqueID, dat$UniqueID )
filter( studies, UniqueID %in% setdiff( studies$UniqueID, dat$UniqueID ) )


ranges <- dat %>% 
  rename( group = weight ) %>%
  group_by( group, runID ) %>%
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
             n = n(),
             .groups = "drop" )

sub_dat <- dat %>% 
  filter( runID %in% sample( unique( ranges$runID ),  5 ) )

sub_dat

ggplot( sub_dat, aes( rangeATE ) ) +
  facet_grid( weight ~ runID ) +
  geom_dotplot( method="histodot", binwidth = 0.01 ) +
  labs( title = "5 replicates of fingerprint for ATE range" ) +
  theme_minimal()

summary( sub_dat$ratioSE )
sub_dat$ratioSE[ sub_dat$ratioSE > 5 ] = 5
ggplot( sub_dat, aes( ratioSE ) ) +
  facet_grid( weight ~ runID ) +
  geom_dotplot( method="histodot" ) +
  labs( title = "5 replicates of fingerprint for SE ratio",
        subtitle = "(Ratios windsorized at 5x)" ) +
  theme_minimal()



# Make the above plot, but with the REAL data also!!! ----








