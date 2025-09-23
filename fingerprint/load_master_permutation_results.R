

# Load all the empirical permutation simulation results and put in
# unified dataset along with the empirical results (on the og data) as
# well

library( tidyverse )
library( clusterRCT )


core.methods <- c("LRa-FE-db", "LRa-FIcw-db", "LRapw-FE-db", "LRapw-FIpw-db", "MLM-FE", "MLM-FIcw", "MLM-RE")


core.methods.se <- c("LRa-FE-db", "LRa-FE-het", "LRa-FIcw-db", "LRapw-FE-db", 
                     "LRapw-FE-het", "LRapw-FIpw-db", "LRapw-FIpw-het", 
                     "LRicw-FE-crve", "LRicw-FIcw-crve", "MLM-FE", "MLM-FIcw",
                     "MLM-RE")



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
SIM_DIR = here::here( "fingerprint/results_rerand_sim_public/" )



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


# Our simulation results
skimr::skim_without_charts( simdat )

simdat$Covariates = "best"
warning( "Setting public covariate to best by hand." )

# Load MDRC side ----

mdrc_dat <- read_rds( here::here( "fingerprint/results_rerand_master/master_run_fingerprints.RDS" ) ) %>%
  as_tibble()

skimr::skim_without_charts( mdrc_dat )

mdrc_dat <- mdrc_dat %>%
  rename( method = Estimator ) %>%
  mutate( UniqueID     = str_remove_all(paste(Project, Outcome, Trt_arm, sep = '_'), '_NA') )

table( mdrc_dat$Project )
table( mdrc_dat$Outcome )
table( mdrc_dat$Trt_arm )
table( mdrc_dat$UniqueID )  

if ( FALSE ) {
  tt <- filter( mdrc_dat, UniqueID == "Tennessee STAR_gktreadss" ) %>%
    #  dplyr::filter( weight == "overall" ) %>%
    dplyr::select( -Covariates, -best_cov, -max_cor ) 
  View( tt )
  
  t_agg <- tt %>% 
    group_by( runID ) %>%
    summarise( n = n(),
               nna = sum( is.na(ATE_hat ) ),
               EATE = mean( ATE_hat, na.rm=TRUE ),
               max = max( ATE_hat, na.rm=TRUE ),
               min = min( ATE_hat, na.rm=TRUE ),
               del = max - min ) %>%
    arrange( del )
  
  t_agg %>%
    knitr::kable( digits = 2 )
  
  sd( t_agg$EATE, na.rm=TRUE )
}




# Combine and clean up ----

names(simdat)
names(mdrc_dat)

simdat = bind_rows( `Curated` = mdrc_dat, `Public` = simdat, 
                    .id="Collection" )
rm( mdrc_dat )

# Clean up the unique IDs
source( here::here( "empirical/R/clean_result_table.R" ) )
simdat <- clean_result_table( simdat )

simdat$Project = NULL
simdat$Outcome = NULL
simdat$Trt_arm = NULL

table( simdat$Collection, simdat$Covariates, useNA="always" )

skimr::skim( simdat )
table( table(simdat$runID ) )




# Standardize our outcomes ----

ok = !is.na( simdat$ATE_hat )
simdat <- left_join( simdat,
                     dplyr::select( studies, UniqueID, sdY0 ),
                     by="UniqueID" ) %>%
  mutate( ATE_hat = ATE_hat / sdY0,
          SE_hat = SE_hat / sdY0 ) %>%
  dplyr::select( -sdY0 )

# Make sure we have not introduced any NAs in our ATE_hat
stopifnot( all( !is.na( simdat$ATE_hat[ok] ) ) )

names( simdat )

# Who is missing from the fingerprint simulation? ----
table( simdat$UniqueID )
length( unique( simdat$UniqueID ) )
miss_sets <- setdiff( studies$UniqueID, simdat$UniqueID )

miss_sets

# The missing ones
filter( studies, UniqueID %in% miss_sets ) %>%
  dplyr::select( UniqueID, n, J ) %>%
  knitr::kable( caption = "Dropped Datasets" ) %>%
  print()

miss_sets <- filter( studies,
                     UniqueID %in% miss_sets )
# rm( miss_sets )

# Subset list of studies to our core fingerprinted studies
studies <- filter( studies,
                   UniqueID %in% simdat$UniqueID )


# Remove unneeded variables -----
rm( SIM_DIR )
rm( clean_result_table, make_df_flags )


# **********************************************************************
# Now calculate useful aggregations
# **********************************************************************


# Calculate ranges ----


range_dat <- simdat %>%
  group_by( Collection, UniqueID, Covariates ) %>%
  nest() %>%
  mutate( data = map( data, summarize_simulation_results, "cross" ) ) %>%
  unnest( data ) 

range_dat  
table(range_dat$runID )


meths = unique( simdat$method )
fi <- str_detect( meths, "-FI" )
names(fi) <- meths
sort(fi)
range_dat_noFI <- simdat %>%
  dplyr::filter( !str_detect( method, "-FI" ) ) %>%
  group_by( Collection, UniqueID, Covariates ) %>%
  nest() %>%
  mutate( data = map( data, summarize_simulation_results, "cross" ) ) %>%
  unnest( data )

  

# Calculate method performances ----

table( simdat$method, useNA="always" )

meth_dat <- simdat %>%
  filter( runID != "empirical" ) %>%
  group_by( Collection, UniqueID, Covariates ) %>%
  nest() %>%
  mutate( data = map( data, summarize_simulation_results, "method" ) ) %>%
  unnest( data ) %>%
  ungroup()

meth_dat  


# Windsorize the SE ranges ----
range_dat$ratioSE[ range_dat$ratioSE > 5 ] = 5

range_dat <- range_dat %>%
  rename( Estimand = weight )




