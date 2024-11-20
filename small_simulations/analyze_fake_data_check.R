

# Quick script to just compare all methods on a single dataset and
# group the ones with the same point estimates.

library( clusterRCT )
data( fakeCRT )
head( fakeCRT )

cc <- compare_methods( Yobs ~ T.x | S.id | D.id,
                 control_formula = ~ C.ijk + X.jk,
                 data = fakeCRT,
                 include_method_characteristics = TRUE )

cc

cc %>% dplyr::select( -weight, -df )

cc %>%
    arrange( ATE_hat ) %>%
    dplyr::select( -weight, -df, -p_value, -SE_hat ) %>%
    mutate( ATE_hat = round( ATE_hat, digits = 7 ),
            dup = duplicated( ATE_hat ) )

# See if some of the estimators are giving the same point estimate for
# this dataset
source( here::here( "simulations/simulation_helper_functions.R" ) )
cc$runID = 1
group_estimators( cc )



# No covariates ---

cc <- compare_methods( Yobs ~ T.x | S.id | D.id,
                       data = fakeCRT,
                       include_method_characteristics = TRUE )

cc

# See if some of the estimators are giving the same point estimate for
# this dataset
source( here::here( "simulations/simulation_helper_functions.R" ) )
cc$runID = 1
group_estimators( cc )



# No blocking ----

cc <- compare_methods( Yobs ~ T.x | S.id,
                       control_formula = ~ C.ijk + X.jk,
                       data = fakeCRT,
                       include_method_characteristics = FALSE )
cc

