

model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 30                          # number of schools
    , K = 10                          # number of districts
    , nbar = 10                       # number of individuals per school
    , S.id = NULL                     # N-length vector of school assignments
    , D.id = NULL                     # N-length vector of district assignments
    , Xi0 = 0                         # scalar grand mean outcome under no treatment
    , MDES = 0.125            # minimum detectable effect size
    , R2.3 = 0.1              # percent of district variation
    , ICC.3 = 0.2 # district intraclass correlation
    , omega.3 = 0.1           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.1           # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)



sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
rm( model.params.list )




test_that("DB estimators work", {
    head( sim.data )
    table( sim.data$D.id)
    aa = design_based_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    aa
    expect_true( nrow( aa ) == 4 )

    bb4 = middleton_aronow_estimator( Yobs ~ T.x | S.id | D.id, data=sim.data )
    expect_true( nrow( bb4 ) == 2 )

    # No block-level.
    sim.data$D.id = NULL
    bb = design_based_estimators( Yobs ~ T.x | S.id, data=sim.data, weight = "Cluster" )
    bb
    expect_true( nrow(bb) == 1 )

    bb2 = design_based_estimators( Yobs ~ T.x | S.id, data=sim.data, weight = "Person" )
    bb2
    expect_true( nrow(bb) == 1 )


    bb3 = middleton_aronow_estimator( Yobs ~ T.x | S.id, data=sim.data )

})






test_that("DB weighting works", {

    ss = tibble( blockID = c( 1, 1, 1, 1, 1, 1 ),
                 clusterID = c( 1, 2, 3, 4, 5, 6 ),
                 Z = c( 0, 1, 0, 0, 1, 0 ),
                 Y0 = c( 9, 10, 11, 14, 15, 16 ),
                 tau = c( 0, 0, 0, 10, 10, 10 ) )
    ss = ss[ rep( 1:6, c(2,2,2,4,4,4) ), ]
    ss$Y0 = ss$Y0 + rep( c(-1, 1), nrow(ss)/2 )

    ss2 = ss
    ss2$clusterID = 1:nrow(ss2)
    ss2$tau = ss2$tau * 5 + 50
    ss2$blockID = 2
    ss = bind_rows( ss, ss2, ss2 ) %>%
        mutate( Yobs = Y0 + Z*tau ) %>%
        arrange( blockID, clusterID )
    ss
    table( ss$blockID )
    describe_clusterRCT( Yobs ~ Z | clusterID | blockID, data=ss )

    mean( ss$tau )

    stats <- ss %>%
        group_by( blockID, clusterID ) %>%
        summarise( n = n(),
                   tau = mean( tau ), .groups = "drop" )
    table( stats$n )

    sts <- stats %>% group_by( blockID ) %>%
        summarise( wt = sum( n ),
                   J = n(),
                   tau_cl = mean( tau ),
                   tau_per = weighted.mean( tau, w=n ) )
    sts
    stats %>%
        summarise( wt = sum( n ),
                   J = n(),
                   tau_cl = mean( tau ),
                   tau_per = weighted.mean( tau, w=n ) )

    bl_avgs <- sts %>% summarise(
                       tau_cl = mean(tau_cl),
                       tau_per = mean(tau_per) )
    bl_avgs

    # Weighting within block by person...
    aa = design_based_estimators( Yobs ~ Z | clusterID | blockID, data=ss, weight = "Person",
                                  include_block_estimates = TRUE )
    aa

    blocks = attr(aa, "blocks" )
    blocks
    expect_equal( as.numeric(blocks$ATE_hat), sts$tau_per )

    expect_equal( aa$ATE_hat[[1]], bl_avgs$tau_per )
    expect_equal( aa$ATE_hat[[2]], weighted.mean( sts$tau_per, w=sts$J ) )
    expect_equal( aa$ATE_hat[[3]], weighted.mean( sts$tau_per, w=sts$wt ) )
    expect_equal( aa$ATE_hat[[4]], weighted.mean( sts$tau_per, w=sts$wt ) )

    aa = design_based_estimators( Yobs ~ Z | clusterID | blockID, data=ss, weight = "Cluster",
                                  include_block_estimates = TRUE )
    blocks = attr(aa, "blocks" )
    expect_equal( as.numeric(blocks$ATE_hat), sts$tau_cl )

    aa
    expect_equal( aa$ATE_hat[[1]], bl_avgs$tau_cl )
    expect_equal( aa$ATE_hat[[2]], weighted.mean( sts$tau_cl, w=sts$J ) )
    expect_equal( aa$ATE_hat[[3]], weighted.mean( sts$tau_cl, w=sts$wt ) )
    expect_equal( aa$ATE_hat[[4]], weighted.mean( sts$tau_cl, w=sts$J ) )

    # No block-level.
    ss
    ss$cid = paste( ss$blockID, ss$clusterID, sep="-" )
    bb = design_based_estimators( Yobs ~ Z | cid, data=ss )
    bb
    expect_true( nrow(bb) == 1 )

} )



test_that( "middleton aronow works", {
    ss = tibble( blockID = c( 1, 1, 1, 1, 1, 1 ),
                 clusterID = c( 1, 2, 3, 4, 5, 6 ),
                 Z = c( 0, 1, 0, 0, 1, 0 ),
                 Y0 = c( 9, 10, 11, 14, 15, 16 ),
                 tau = c( 0, 0, 0, 10, 10, 10 ) )
    ss = ss[ rep( 1:6, c(2,2,2,4,4,10) ), ]
    ss$Y0 = ss$Y0 + rep( c(-1, 1), nrow(ss)/2 )

    ss2 = ss
    ss2$blockID = 2
    ss = bind_rows( ss, ss2, ss2 ) %>%
        mutate( Yobs = Y0 + Z*tau ) %>%
        arrange( blockID, clusterID )
    ss

    # Raj estimator, etc
    aa = middleton_aronow_estimator( Yobs ~ Z | clusterID | blockID, data=ss,
                                  include_block_estimates = TRUE )

    blk <- attr( aa, "blocks" )
    blk
    expect_true( nrow( blk ) == 2 )

    expect_true( aa$SE_hat[[2]] < aa$SE_hat[[1]] )

})
