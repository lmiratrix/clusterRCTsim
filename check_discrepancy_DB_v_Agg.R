

# Identifying datasets where Agg and DB are not being the same on when they provide standard errors

source( "empirical/setup.R" )

combined

cc = combined %>%
  filter( Estimator %in% c( "Agg_FI_Cluster_Cluster", "DB_FI_Cluster_Cluster" ) ) %>%
  dplyr::select( Collection, UniqueID, Covariates, Estimator, df, ATE_hat, SE_hat ) %>%
  pivot_wider( names_from = Estimator, values_from = c( ATE_hat, SE_hat ) ) %>%
  rename_with( ~ str_replace(., "_Cluster_Cluster", "" ) ) %>%
  mutate( discATE = is.na( ATE_hat_Agg_FI ) != is.na( ATE_hat_DB_FI ) ) %>%
  mutate( discSE = is.na( SE_hat_Agg_FI ) != is.na( SE_hat_DB_FI ) )

cc %>%
  filter( discSE ) %>%
  print( n = 100 )

cc %>%
  filter( discSE ) %>%
  pull( UniqueID ) %>% unique()
