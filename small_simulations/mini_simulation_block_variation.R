


# Simple simulation as sanity check of SEs etc.
#
# These simulations generally have two sizes of clusters, with big
# clusters having larger treatment impacts.
#
#
# This is a "superpopulation" simulation--each iteration is a new
# dataset.


set.seed( 1039 )
library( tidyverse )
library( clusterRCT )

source( here::here( "small_simulations/simulation_DGP.R" ) )



#### Simple simulation ####

# This simulation generates data and compares the estimators to the 4
# different estimands to see which are biased.
#
# The third DGP ("simple") is a special DGP that hand-builds blocks of
# different sizes to ensure that there are correlations with impact
# and size of block and number of blocks in a district, and so forth.
#
# The other DGPs are more standard two-tier DGPs with different
# intraclass correlations and effect size variability.

if ( FALSE ) {
    cat( "Running 'Simple simulation' in mini_simulation_varied.R\n" )

    R = 100

    #    res1 <- run_simple_simulation( ICC.2 = 0.4, omega.2 = 1, R = R )
    #    res2 <- run_simple_simulation( ICC.2 = 0, omega.2 = 0, R = R )


    res3 <- run_simple_simulation( simple = "B", R = R )
    res <- res3 # hack to skip res1 and res2
    res$sim = "simple"

    #res <- bind_rows( varied=res1, not_varied=res2, simple=res3, .id="sim" )

    res %>%
        arrange( EATE ) %>%
        knitr::kable(digits=2)

    res$ESE_hat[ is.na( res$ESE_hat ) ] = 0
    res <- mutate( res, method = fct_reorder2( method, ESE_hat, EATE ) )
    res
    #res <- mutate( res, method = factor( method,
    #                                     levels = rev( sort( levels(method) ) ) ) )

    # Look at third simulation only
    res %>%
        filter( sim == "simple" ) %>%
        ggplot( aes( method, EATE ) ) +
        geom_hline( aes( yintercept = mean( tau_c ),  lty="clean", col="Cluster" ) ) +
        geom_hline( aes( yintercept = mean( tau_p ), lty="clean", col="Person" ) ) +
        geom_hline( aes( yintercept = mean( tau_cb ),  lty = "block", col="Cluster" ) ) +
        geom_hline( aes( yintercept = mean( tau_pb ), lty="block", col="Person" ) ) +
        geom_point(size=3) +
        #        geom_errorbar( aes( ymin = EATE - ESE_hat, ymax = EATE + ESE_hat ), width=0) +
        geom_errorbar( aes( ymin = EATE - SE_E, ymax = EATE + SE_E ), linewidth = 1, width=0) +
        theme_minimal() +
        coord_flip()


}


#### Looking at what estimators align with other estimators ####

# This runs a few iterations
if ( FALSE ) {

    R = 5

    #rps = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = TRUE ),
    #              .id = "runID", .progress=TRUE )

    INC_COV = TRUE

    rps1 = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = FALSE, simple = "B",
                                   include_covariates = INC_COV ),
                   .id = "runID", .progress=TRUE )
    rps2 = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = TRUE, simple = FALSE,
                                   include_covariates = INC_COV ),
                   .id = "runID", .progress=TRUE )
    rps3 = map_df( 1:R, ~ one_run( blocks = TRUE, het_block = TRUE, simple = "A",
                                   include_covariates = INC_COV ),
                   .id = "runID", .progress=TRUE )

    rps = bind_rows( A=rps1, B=rps2, C=rps3, .id = "batch" )
    rps$runID = paste0( rps$batch, rps$runID )
    head( rps )

    saveRDS( rps, here::here( "simulations/pointwise_similarity_results.rds" ) )

    rps = readRDS( here::here( "simulations/pointwise_similarity_results.rds" ) )

    source( here::here( "simulations/simulation_helper_functions.R" ) )
    groups = group_estimators( rps, tolerance = 0.0001 )
    groups


    if ( FALSE ) {
        # Linear regression vs. DB?
        g = group_estimators( rps, tolerance = 0.01, return_matrix = TRUE )
        g[ "LR_FI_CRVE_Person", ]

        rW = filter( rps, method %in% c( "DB_FI_Person_Person", "LR_FI_CRVE_Person" ) ) %>%
            dplyr::select( runID, method, ATE_hat ) %>%
            pivot_wider( names_from = method, values_from = ATE_hat )
        rW$diff = rW$DB_FI_Person_Person - rW$LR_FI_CRVE_Person
        summary( rW$diff )
    }

    #groups = group_estimators( rps, tolerance = 0.05, return_matrix = FALSE )
    #summary( as.numeric( groups ) )

    # List the groups
    map_chr( groups, paste, collapse = ", " ) %>%
        as.data.frame() %>%
        mutate( grp = 1:n() ) %>%
        relocate( grp ) %>%
        knitr::kable()


    if ( FALSE ) {
        # Make list of estimators with group IDs (optional)
        group_df <- map_df(seq_along(groups), ~tibble(estimator = groups[[.x]], group = .x))

        print(group_df)
        table( group_df$group )

    }
}





