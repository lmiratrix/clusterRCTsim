

library( tidyverse )



# Ideas
# Scatterplot of number of blocks vs. range in estimates
# Scatterplot of some measure of "baseline DF" vs. range in estimates
#    Baseline DF is the number clusters - 2 * # blocks - # covariates

# Binary outcomes as effect sizes
0.05 / sqrt( 0.25 * 0.75 )


# How do we compare the different groups of 4 estimands?  And how do we track where MLM is in this?


dat = readxl::read_excel( "empirical/Appendix Table 1. Estimates of ITT Effect and its SE.xlsx" )
head( dat ) )
names( dat )
table( dat$Covariates )

dat <- mutate( dat,
               ATE_hat = ATE_hat_sd,
               SE_hat = SE_hat_sd ) %>%
  dplyr::select( -ATE_hat_sd, -SE_hat_sd )

table( dat$Estimator )

agg <- dat %>% 
  group_by( UniqueID, Covariates ) %>%
  summarise( E_ATE = mean( ATE_hat ),
             sd_ATE = sd( ATE_hat ),
             min_ATE = min( ATE_hat ),
             max_ATE = max( ATE_hat ),
             range_ATE = max_ATE - min_ATE,
             E_SE = mean( SE_hat ),
             min_SE = min( SE_hat ),
             max_SE = max( SE_hat ),
             ratio_SE = max_SE / min_SE, 
             .groups = "drop" ) 

agg

summary( agg$range_ATE )

ggplot( agg, aes( sd_ATE ) ) +
  facet_wrap( ~ Covariates ) +
  geom_histogram( bins = 20 ) +
  labs( x = "Range of ATE estimates" ) +
  theme_minimal()


dat$Covariates
d2 <- dat %>% 
  filter( !is.na( ATE_hat ) ) %>%
  group_by( UniqueID, Covariates ) %>%
  mutate( med_ATE = median( ATE_hat, na.rm=TRUE ),
          err = abs( ATE_hat - med_ATE ),
          .groups = "drop" ) %>%
  group_by( Covariates, Estimator ) %>%
  summarise( E_err = mean( err ),
             sd_err = sd( err ),
             min_err = min( err ),
             max_err = max( err ),
             range_err = max_err - min_err,
             .groups = "drop" )
d2

ggplot( d2, aes( Estimator, sd_err ) ) +
  facet_wrap( ~ Covariates ) +
  geom_bar( stat = "identity" ) +
  coord_flip()
