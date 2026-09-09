
# First load empirical results

source( here::here( "empirical/load_empirical_results.R" ) )
library( clusterRCT )


# Load the fingerprint simulation results (where we generated data
# that looked like our set of studies).

core.methods <- c("LRa-FE-db", "LRa-FIcw-db", "LRapw-FE-db", "LRapw-FIpw-db", "MLM-FE", "MLM-FIcw", "MLM-RE")


core.methods.se <- c("LRa-FE-db", "LRa-FE-het", "LRa-FIcw-db", "LRapw-FE-db", 
                     "LRapw-FE-het", "LRapw-FIpw-db", "LRapw-FIpw-het", 
                     "LRicw-FE-crve", "LRicw-FIcw-crve", "MLM-FE", "MLM-FIcw",
                     "MLM-RE")


studies <- filter( studies,
                   blocked == "Blocked" )
studies



# Load and examine data ----

cli::cli_h1( "!!! Loading and examining results" )

if ( FALSE ) {
  
  # Read one example result file
  res = read_rds( here::here( "simulationII/results_fingerprint",
                              paste0( studies$UniqueID[1], ".rds" ) ) )
  res
  table( res$runID )  
}


dat <- studies$UniqueID %>%
  set_names() %>%
  map( \(file) {
    fname = here::here( "simulationII/results_fingerprint", paste0( file, ".rds" ) )
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
table( dat$runID )
length( unique( dat$runID ) )

# Hard code that all use a single covariate
dat$Covariates = "best"

# Who is missing?
missing_sets <- setdiff( studies$UniqueID, dat$UniqueID )
filter( studies, UniqueID %in% missing_sets ) %>%
  dplyr::select( UniqueID, n, J )

dat <- left_join( dat, 
                  select( studies, UniqueID, Collection, num_singletons),
                  by = "UniqueID" )



studies <- filter( studies,
                   UniqueID %in% unique( dat$UniqueID ) )





# Zap CRVE estimator SEs for singletons ----


drp <- dat$method %in% c( "LRicw-FIcw-crve",
                          "LRi-FIpw-crve" ) &
  dat$num_singletons > 0
table( drp )

dat$SE_hat[ drp ] <- NA

rm( drp )




# Calculate ranges, etc. ----

range_dat <- dat %>%
  group_by( Collection, UniqueID ) %>%
  nest() %>%
  mutate( data = map( data, summarize_simulation_results, "cross" ) ) %>%
  unnest( data ) 


ranges <- range_dat %>% 
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


set.seed( 403443 )
sub_dat <- range_dat %>% 
  ungroup() %>%
  filter( runID %in% sample( unique( ranges$runID ),  5 ) )

sub_dat

# Distribution of ATEs for all methods across studies
if ( FALSE ) {
  ggplot( dat, aes( UniqueID, ATE_hat ) ) +
    geom_jitter( size=0.25, width=0.1, height=0 ) +
    coord_flip()
}



# Calculate method performances ----

table( dat$method, useNA="always" )
head( dat )

dat <- left_join( dat, 
                  method_characteristics() %>%
                    dplyr::select( method, biased ),
                  by="method" )

meth_dat <- dat %>%
  group_by( Collection, Covariates, UniqueID ) %>%
  nest() %>%
  mutate( data = map( data, summarize_simulation_results, "method" ) ) %>%
  unnest( data ) %>%
  ungroup()

meth_dat  
if ( FALSE ) {
  qplot( meth_dat$p.noSE )
}


# Load the empirical estimates for baseline ---

empirical_estimates <- combined
rm( combined )

rm( outs )
rm( unblocked )
rm( PLOT_PATH )
rm( make_df_flags, plot_path )
