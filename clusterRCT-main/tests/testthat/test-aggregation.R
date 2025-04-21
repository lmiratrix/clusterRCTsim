

test_that( "aggregation works", {
    set.seed( 433434 )

    dd = data.frame( Y = 1:20,
                     clusterID = sample(LETTERS[1:5], 20, replace=TRUE),
                     blockID = sample( LETTERS[1:2], 20, replace=TRUE ),
                     x = 1:20,
                     x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd$Z = randomizr::cluster_ra( clusters = dd$clusterID )
    dd = clusterRCT:::make_canonical_data( Y ~ Z | clusterID | blockID, data=dd,
                                           control_formula = ~ x + x2 )
    dd$Z = randomizr::cluster_ra( clusters = dd$clusterID )
    aa <- clusterRCT:::aggregate_data( data = dd )
    aa
    expect_equal( names(aa), c("blockID", "clusterID","Z", "Ybar", "n"))

    bb <-     clusterRCT:::aggregate_data( dd, control_formula = ~ x + x2 )
    bb
    expect_true( all( c( "x", "x2B", "x2C" ) %in% names(bb) ) )
    expect_true( nrow(bb) == nrow( aa ) )
    expect_true( all( bb$Ybar == aa$Ybar ) )


    dd$blockID = NULL
    a2 = clusterRCT:::aggregate_data( data = dd, control_formula = ~ x + x2 )
    a2

    a4 = clusterRCT:::aggregate_data( data = dd, control_formula = ~ x )
    a4
    expect_true( all( a2$x == a4$x ) )
})


test_that( "aggregation estimators work", {

    data( fakeCRT )
    fakeCRT

    J = length( unique( fakeCRT$S.id ) )
    J

    e2 <- aggregation_estimators( Yobs ~ T.x | S.id | D.id,
                                  #control_formula = ~ X.jk + C.ijk,
                                  data = fakeCRT )
    expect_true( is.data.frame(e2) )

    e2cov <- aggregation_estimators( Yobs ~ T.x | S.id | D.id,
                                     control_formula = ~ X.jk + C.ijk,
                                     data = fakeCRT )
    expect_true( is.data.frame(e2cov) )
    expect_true( all( e2$ATE_hat != e2cov$ATE_hat ) )


    # with covariate interactions (no district)
    e3 <- aggregation_estimators( Yobs ~ T.x | S.id,
                                  control_formula = ~ X.jk + C.ijk,
                                  control_interacted = TRUE,
                                  data = fakeCRT )
    expect_true( is.data.frame(e3) )
})




