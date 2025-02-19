
# Script purpose:
#
# Identify datasets where Agg and DB are not being the same on when
# they provide standard errors

source( "empirical/load_empirical_results.R" )

combined
table( combined$Estimator )

cc = combined %>%
  filter( Estimator %in% c( "LRa-FIcw-het", "LRa-FIcw-db" ) ) %>%
  dplyr::select( Collection, UniqueID, Covariates, Estimator, df, ATE_hat, SE_hat ) %>%
  pivot_wider( names_from = Estimator, values_from = c( ATE_hat, SE_hat ) ) %>%
  rename_with( ~ str_replace(., "LRa-FIcw-", "" ) ) %>%
  mutate( discATE = is.na( ATE_hat_het ) != is.na( ATE_hat_db ) ) %>%
  mutate( discSE = is.na( SE_hat_het ) != is.na( SE_hat_db ) )
cc

summary( cc$SE_hat_db - cc$SE_hat_het )

sample_n(cc, 20 )
cc %>%
  filter( discSE ) %>%
  print( n = 100 )

cse <- cc %>%
  filter( discSE )
cse

cse %>%
  pull( UniqueID ) %>% unique()








