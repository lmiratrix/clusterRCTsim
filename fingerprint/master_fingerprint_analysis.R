
# Analysis of the permutation simulations done on the empirical data

# This script analyzes the permutation simulations done on the
# collection of empirical data we have across both the curated and
# non-curated (with some datasets missing due to access issues).

# This script, among other things, makes plots showing some of the
# fingerprint runs along with the empirical data to see if we are getting
# similar ranges and whatnot in our simulation to empirical.


library( tidyverse )
library( clusterRCT )


source( here::here( "fingerprint/load_master_fingerprint_results.R" ) )



# Calculate ranges for each outcome for each simulation trial ---


range_dat <- simdat %>%
  group_by( UniqueID ) %>%
  nest() %>%
  mutate( data = map( data, summarize_simulation_results, "cross" ) ) %>%
  unnest( data ) 

# Windsorize the SE ranges
range_dat$ratioSE[ range_dat$ratioSE > 5 ] = 5



range_dat <- range_dat %>%
  rename( Estimand = weight )





# Calculate aggregate statistics on ranges across outomes for each run ----


ranges <- range_dat %>% 
  group_by( Estimand, runID ) %>%
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
ranges




# Make the combined plots ----

if ( FALSE ) {
  
  # See the Rmd file for most recent versions of this
  
  set.seed( 4044404 )
  
  selected_sims = sample( unique( range_dat$runID ), 4 )
  selected_sims
  
  # Add empirical data results to simulated data results
  range_dat = bind_rows( range_dat, 
                      filter( summary_stats, Covariates == "none" ) ) %>%
    mutate( empirical = ifelse( runID == "empirical", "empirical", "simulated" ) )
  
  
  sub_dat <- range_dat %>% 
    filter( runID %in% c( selected_sims, "empirical" ) )
  
  
  # Make hidden IDs
  sub_dat$blindRunID = factor( sub_dat$runID,
                               levels = sample( unique(sub_dat$runID) ),
                               labels = 1:5 )
  
  
  # Which is the empirical one??
  sub_dat %>% 
    dplyr::select( runID, blindRunID ) %>%
    distinct() 
  
  
  aterange_plot <- ggplot( sub_dat, aes( rangeATE ) ) +
    facet_grid( Estimand ~ blindRunID ) +
    geom_dotplot( method="histodot" ) +
    labs( title = "5 replicates of fingerprint for ATE ranges",
          subtitle = "(true data also present)" ) +
    theme_minimal()
  aterange_plot
  ggsave( here::here( "fingerprint_ATE_ranges.pdf" ), width = 8, height = 4 ) 
  
  
  seratio_plot <- ggplot( sub_dat, aes( ratioSE ) ) +
    facet_grid( Estimand ~ blindRunID ) +
    geom_dotplot( method="histodot" ) +
    labs( title = "5 replicates of fingerprint for SE ratio",
          subtitle = "(Ratios windsorized at 5x; true data also present)" ) +
    theme_minimal()
  seratio_plot
  
  ggsave( here::here( "fingerprint_SE_ratios.pdf" ), width = 8, height = 4 ) 
  
  
  
  table( range_dat$runID )
  
}





