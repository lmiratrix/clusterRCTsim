

# Load all the empirical simulation results and put in unified dataset
# with the empirical results as well

library( tidyverse )
library( clusterRCT )



# Load list of empirical studies ----

studies = read_csv( here::here( "empirical/results/study_descriptions.csv" ),
                    show_col_types = FALSE )
studies <- filter( studies,
                   blocked == "Blocked" )
studies



# Load Josh dataset side ----

# Simulation result directory
# This directory holds the results of some simulation, with one file
# per dataset.
SIM_DIR = here::here( "fingerprint/rerand_sim_results/" )



simdat <- studies$UniqueID %>%
  set_names() %>%
  map( \(file) {
    fname = paste0( SIM_DIR, file, ".rds" )
    # check if file exists
    if ( file.exists( fname ) ) {
      res <- read_rds( fname )
      #summarize_simulation_results( res, "cross" )
      res
    } else {
      NULL
    }
  } )  %>% 
  bind_rows( .id = "UniqueID" ) %>%
  as_tibble()

simdat$runID = paste0( "R", simdat$runID )

# Our simulation results
simdat



# Load MDRC side ----

mdrc_dat <- read_rds( "fingerprint/results_master_finger/master_run_fingerprints.RDS" ) %>%
  as_tibble()
mdrc_dat


mdrc_dat <- mdrc_dat %>%
  mutate( UniqueID     = str_remove_all(paste(Project, Outcome, Trt_arm, sep = '_'), '_NA') )

table( mdrc_dat$Project )
table( mdrc_dat$Outcome )
table( mdrc_dat$Trt_arm )
table( mdrc_dat$UniqueID )  
          


# Combine and clean up ----

names(simdat)
names(mdrc_dat)

simdat = bind_rows( `Curated` = mdrc_dat, `Public` = simdat, 
                    .id="Collection" )
rm( mdrc_dat )

source( here::here( "empirical/R/clean_result_table.R" ) )
simdat <- clean_result_table( simdat )

simdat
table( table(simdat$runID ) )

# Standardize ----

ok = !is.na( simdat$ATE_hat )
simdat <- left_join( simdat,
                     dplyr::select( studies, UniqueID, sdY0 ),
                     by="UniqueID" ) %>%
  mutate( ATE_hat = ATE_hat / sdY0,
          SE_hat = SE_hat / sdY0 ) %>%
  dplyr::select( -sdY0 )

stopifnot( all( !is.na( simdat$ATE_hat[ok] ) ) )



# Who is missing from the fingerprint simulation?
table(simdat$UniqueID)
length( unique( simdat$UniqueID ) )
miss_sets <- setdiff( studies$UniqueID, simdat$UniqueID )

miss_sets


# The missing ones
filter( studies, UniqueID %in% miss_sets ) %>%
  dplyr::select( UniqueID, n, J ) %>%
  knitr::kable()

rm( miss_sets )

# Subset to our core fingerprinted studies
studies <- filter( studies,
                   UniqueID %in% simdat$UniqueID )

# some final cleanup
rm( SIM_DIR )
rm( clean_result_table, make_df_flags )



