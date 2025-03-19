

# Make plots showing 4 of the fingerprint runs along with the real
# data to see if we are getting similar ranges and whatnot in our
# simulation to empirical

library( tidyverse )



# Load simulation data ----

studies = read_csv( here::here( "empirical/results/study_descriptions.csv" ),
                    show_col_types = FALSE )
studies <- filter( studies,
                   blocked == "Blocked" )
studies


if ( FALSE ) {
  
  # Read one example result file
  res = read_rds( here::here( "fingerprint_results", paste0( studies$UniqueID[1], ".rds" ) ) )
  res
  
}

simdat <- studies$UniqueID %>%
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


# Who is missing?
miss_sets <- setdiff( studies$UniqueID, simdat$UniqueID )
miss_sets

filter( studies, UniqueID %in% setdiff( studies$UniqueID, simdat$UniqueID ) )

simdat$ratioSE[ simdat$ratioSE > 5 ] = 5

simdat <- simdat %>%
  rename( Estimand = weight )

  




# Load the empirical results ----

# load all data
source( here::here( "empirical/load_empirical_results.R" ) )
source( here::here( "empirical/makeDotplot.R" ) )

nrow( combined )
combined <- filter( combined,
                    !( UniqueID %in% miss_sets ) )
nrow( combined )

source( here::here( "empirical/calculate_all_ranges.R" ) )

rm( combined, combined_exclFI )
rm( outs )
rm( unblocked )

names(summary_stats_estimand)
summary_stats$Estimand = "overall"

summary_stats = bind_rows( summary_stats, summary_stats_estimand )

summary_stats <- summary_stats %>%
  mutate( Estimand = fct_recode( Estimand,
                                 "C" = "Cluster",
                                 "P" = "Person" ) )





summary_stats <- summary_stats %>%
  rename( rangeATE = Range ,
          ratioSE = Ratio )
summary_stats$runID = "real"

simdat = bind_rows( simdat, 
                    filter( summary_stats, Covariates == "none" ) ) %>%
  mutate( real = ifelse( runID == "real", "real", "simulated" ) )





# Make the combined plots ----

set.seed( 4044404 )
sub_dat <- simdat %>% 
  filter( runID %in% c( sample( unique( ranges$runID ), 4 ), "real" ) )


# Make hidden IDs
sub_dat$blindRunID = factor( sub_dat$runID,
                             levels = sample( unique(sub_dat$runID) ),
                             labels = 1:5 )


# Which is the real one??
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



table( simdat$runID )




# Calculate aggregate statistics ----



ranges <- simdat %>% 
  group_by( Estimand, runID, real ) %>%
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




