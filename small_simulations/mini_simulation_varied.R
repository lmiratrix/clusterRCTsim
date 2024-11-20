


# Simple simulation as sanity check of SEs etc.
#
# This sim has two sizes of clusters, with big clusters having larger
# treatment impacts.
#
# Inspired by Mike Weiss questions on MLM
#
# This is a "superpopulation" simulation--each iteration is a new
# dataset.


set.seed( 1039 )
library( tidyverse )
library( clusterRCT )

source( here::here( "simulations/simulation_DGP.R" ) )


one_run <- function( blocks = TRUE, het_block = FALSE, simple = FALSE, ICC.2 = NULL,
                     include_covariates = FALSE,
                     model.params = model.params.list ) {


    # Make some data
    p = model.params
    if ( !is.null( ICC.2 ) ) {
        p$ICC.2 = ICC.2
    }

    if ( simple != FALSE ) {
        sim.data = make_simple_block_example( version = simple, var_tx = TRUE )
    } else {
        sim.data = make_two_tier_data( n1 = 10, ATE1 = 0,
                                       n2 = 40, ATE2 = 0.4,
                                       params = p,
                                       blocks = blocks, het_block = het_block )
    }



    # Analyze the data
    form = Yobs ~ Tx | S.id | D.id
    if ( !blocks ) {
        form = Yobs ~ Tx | S.id
    }
    control_form = NULL
    if ( include_covariates ) {
        if ( blocks ) {
            control_form = ~ X.jk + C.ijk
        } else {
            control_form = ~ C.ijk
        }
    }

    c1 <- clusterRCT::compare_methods( form, data=sim.data,
                                       control_formula = control_form,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = TRUE,
                                       include_dumb = FALSE,
                                       include_method_characteristics = FALSE )


    # Add the "canonical weight" estimators as a point of comparison
    c2 <- clusterRCT::canonical_weight_estimators( form, sim.data )

    c1 = bind_rows( c1, c2 )

    if ( !blocks ) {
        sim.data$D.id = "sing"
    }


    # Add in the true ATEs
    ds = sim.data %>%
        group_by( D.id, S.id ) %>%
        summarise( tauk = mean( Y1 - Y0 ),
                   n = n(), .groups = "drop" ) %>%
        group_by( D.id ) %>%
        summarise( tau_c = mean( tauk ),
                   tau_p = weighted.mean( tauk, w=n ),
                   J = n(),
                   n = sum( n ) ) %>%
        summarise( tau_cb = mean( tau_c ),
                   tau_c = weighted.mean( tau_c, w=J ),
                   tau_pb = mean( tau_p ),
                   tau_p = weighted.mean( tau_p, w=n ),
                   J = sum( J ),
                   n = sum( n ) )

    c1$tau_c = ds$tau_c
    c1$tau_p = ds$tau_p
    c1$tau_cb = ds$tau_cb
    c1$tau_pb = ds$tau_pb
    c1$J = ds$J
    c1$n = ds$n
    c1$sdY0 = sd( sim.data$Y0 )

    # Return the results
    c1 %>%
        dplyr::select( -any_of( c( "weight", "df" ) ) )
}



if ( FALSE ) {
    # testing
    one_run( simple=TRUE )

    one_run( blocks = TRUE, het_block = TRUE )


    sim.data = make_two_tier_data( n1 = 10, ATE1 = 0,
                                   n2 = 40, ATE2 = 0.4,
                                   params = model.params.list,
                                   blocks = TRUE, het_block = TRUE )

    c1 <- clusterRCT::compare_methods( form = Yobs ~ Tx | S.id | D.id,
                                       data=sim.data,
                                       include_DB = TRUE,
                                       include_MLM = TRUE,
                                       include_LM = TRUE,
                                       include_agg = TRUE,
                                       include_method_characteristics = FALSE,
                                       include_dumb = TRUE )

    c1
}


run_simple_simulation <- function( ICC.2 = 0.4, omega.2 = 1,
                                   simple = FALSE, R = 1000 ) {

    model.params.list$ICC.2 = ICC.2
    model.params.list$omega.2 = omega.2

    if ( TRUE ) {
        library( furrr )
        plan(multisession, workers = parallel::detectCores() - 1 )
        tictoc::tic()
        rps = future_map_dfr( 1:R, ~ one_run( blocks = TRUE,
                                              het_block = TRUE,
                                              simple = simple,
                                              model.params = model.params.list ),
                              .id = "runID", .progress=TRUE,
                              .options = furrr_options( seed=TRUE ) )
        tictoc::tic()
    } else {
        rps = map_df( 1:R, ~ one_run( blocks = TRUE,
                                      het_block = TRUE,
                                      simple = simple,
                                      model.params = model.params.list ),
                      .id = "runID", .progress=TRUE )
    }


    head( rps )
    sd( rps$tau_c)

    res <- rps %>%
        group_by( method ) %>%
        summarise( EATE = mean( ATE_hat ),
                   tau_p = mean( tau_p ),
                   tau_c = mean( tau_c ),
                   tau_cb = mean( tau_cb ),
                   tau_pb = mean( tau_pb ),
                   SE = sd( ATE_hat ),
                   SE_E = sd( ATE_hat ) / sqrt(n()),
                   ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
        mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )

    res
}


if ( FALSE) {

    debug( one_run )
    one_run( simple = TRUE )
}


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




#### Mike's ICC Simulation ####

# Mike wanted to see how the MLM estimators targeted different
# estimands based on the ICC

if ( FALSE ) {

    run_sim <- function( ICC, R = 100 ) {
        cat( "ICC ", ICC, "\n" )

        rps = map_df( 1:R, ~ one_run( ICC.2 = ICC, blocks = FALSE ),
                      .id = "runID", .progress=TRUE )

        head( rps )


        res <- rps %>%
            group_by( method ) %>%
            summarise( EATE = mean( ATE_hat ),
                       tau_p = mean( tau_p ),
                       tau_c = mean( tau_c ),
                       SE = sd( ATE_hat ),
                       ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
            mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )


        res
    }

    #one_run( blocks=FALSE, ICC.2 = 0.8 )

    ICCs = c( 0, 0.2, 0.4, 0.6, 0.8 )
    res <- ICCs %>%
        set_names() %>%
        map_dfr( run_sim, R = 100, .id = "ICC" )

    res

    # Make families of estimators for grouping on the plot
    res$meth_group = str_extract( res$method, "^[^_]+" )

    # Note how all methods are level w.r.t. the ICC, except MLM which
    # goes from person-weighted to cluster-weighted as ICC increases.
    ggplot( res, aes( ICC, EATE ) ) +
        facet_wrap( ~ meth_group, nrow=1 ) +
        scale_y_continuous( limits = c(0, 0.5 ) ) +
        geom_line( aes( col = method, group=method )) +
        geom_hline( aes( yintercept = mean( tau_c ), lty="Cluster" ), col="black" ) +
        geom_hline( aes( yintercept = mean( tau_p ), lty="Person" ), col="black" )


}


# Checking ICC is working
if ( FALSE ) {

    rr = rerun( 100, {

        p = model.params.list
        p$ICC.2 = 0.8
        dat = make_two_tier_data( 10, 0.2, params=p, blocks = FALSE )
        m <- lmer( Yobs ~ 1 + Tx + (1|S.id), data=dat )
        # arm::display(m)
        ICC1 = as.numeric( VarCorr(m)$S.id / (VarCorr(m)$S.id + sigma(m)^2 ) )

        p = model.params.list
        p$ICC.2 = 0
        dat = make_two_tier_data( 10, 0.2, params=p, blocks = FALSE )
        m2 <- lmer( Yobs ~ 1 + Tx + (1|S.id), data=dat )
        # arm::display(m2)
        ICC2 = as.numeric( VarCorr(m2)$S.id / (VarCorr(m2)$S.id + sigma(m2)^2 ) )

        res <- bind_rows( `0.8`=fixef( m ), `0`=fixef( m2 ), .id = "mod" )
        res$ICC = c( ICC1, ICC2 )
        res$tau = mean(dat$Y1 - dat$Y0 )
        res
    } )
    rr = bind_rows( rr )
    rr
    rr %>% group_by( mod ) %>%
        summarise( mean_tx = mean( Tx ),
                   SE = sd( Tx ),
                   SE_Etx = sd( Tx ) / sqrt(n()),
                   mean_ICC = mean(ICC),
                   tau_p = mean(tau),
                   sd_t = sd(tau) )
}



