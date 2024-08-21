
# Test designs with singleton clusters treated, or doubleton, etc.
#
# generally we should always get an ATE, but we might not get a SE,
# unless we are doing covariate adjustment, which can then bring down
# some of the estimators.
#
# Describe data should note these as unusual designs in the printout.


test_that( "matched pairs designs get handled right", {

    blk = blkvar::generate_blocked_data( n_k = c( 5, 10, 20, 10, 5, 10 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:6 ),
                   D = rep( 1:3, each=2 ),
                   Z = rep( c(0,1), 3 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    blk$X = rnorm( nrow(blk) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk, control_formula = ~ X )
    dsc
    as.data.frame(dsc)
    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_equal( dsc$num_singletons, 6 )
    expect_true( dsc$matched_pairs )
    expect_equal( dsc$num_doubletons, 0 )


    expect_warning( expect_warning( comp <- compare_methods( Yobs ~ Z | B | D, data=blk,
                                                             include_method_characteristics = FALSE ) ) )

    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )
    expect_true( all( comp$df[ comp$method %in% c( "DB_FE_Person", "DB_FE_Cluster" ) ] == 2 ) )

    head(blk)
    suppressWarnings( comp <- compare_methods( Yobs ~ Z | B | D, data=blk, control_formula = ~ X ) )
    comp %>%
        filter( is.na( ATE_hat ) )
    expect_true( any( is.na( comp$ATE_hat ) ) )

})


test_that( "singleton tx or co designs get handled right", {

    set.seed( 4343 )
    blk = blkvar::generate_blocked_data( n_k = c( 5, 10, 20, 10, 11, 13,
                                                  10, 10, 20, 10,
                                                  30, 10, 5, 10, 30, 15 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:16 ),
                   D = rep( 1:3, c(6,4,6) ),
                   Z = randomizr::block_ra(blocks=D, block_m = c(3,1,3) ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )
    expect_true( !dsc$matched_pairs )
    expect_equal( dsc$num_singletons, 1 )

    comp = compare_methods( Yobs ~ Z | B | D, data=blk, include_method_characteristics = FALSE )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )

    # This must depend on amount of cross block variation?  For this
    # seed, df is small. But it can be bigger!
    expect_true( comp$df[ comp$method == "MLM_RIRC" ] <= 3 )
    expect_true( comp$df[ comp$method == "MLM_FIRC" ] <= 3 )


    # One block with 2 clusters
    n_k = c( 5,  10, 10, 20, 30, 40,
             20, 10,
             10, 20, 30, 10, 40, 10,
             40, 10, 5, 10, 20, 30 )
    blk = blkvar::generate_blocked_data( n_k = n_k )
    head( blk )
    blks = tibble( B = paste0( "B", 1:length(n_k) ),
                   D = rep( c( 1, 2, 3, 4 ), c(6,2,6,6) ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )
    expect_equal( dsc$num_singletons, 2 )
    expect_equal( dsc$num_doubletons, 0 )

    comp2 = compare_methods( Yobs ~ Z | B | D, data=blk, include_method_characteristics = FALSE )
    comp2
    expect_true( all( !is.na( comp2$ATE_hat ) ) )
    expect_true( comp$df[ comp$method == "MLM_RIRC" ] <= 3 )
    expect_true( comp$df[ comp$method == "MLM_FIRC" ] <= 3 )

    left_join( comp[ c(1:3) ], comp2[ c(1:3) ], by="method" )
    expect_equal( is.na( comp$SE_hat ), is.na( comp2$SE_hat ) )
})




test_that( "doubleton tx or co designs get handled right", {

    blk = blkvar::generate_blocked_data( n_k = c( 5,  10, 10, 20, 30, 40,
                                                  20, 10,  5, 10, 10,
                                                  10, 20, 30, 10, 40, 10 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:17 ),
                   D = rep( c( 1, 2, 3 ), c(6,5,6) ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )
    expect_equal( dsc$num_doubletons, 1 )

    comp = compare_methods( Yobs ~ Z | B | D, data=blk )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )
    expect_true( all( !is.na( comp$SE_hat ) ) )


    # All blocks have  2 tx and 2 co units in it.
    n_k = c( 5,  10, 30, 40,
             20, 10,  5, 10,
             10, 20, 40, 10 )
    blk = blkvar::generate_blocked_data( n_k = n_k )
    head( blk )
    blks = tibble( B = paste0( "B", 1:length(n_k) ),
                   D = rep( c( 1, 2, 3 ), c(4,4,4) ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    dsc
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )
    expect_equal( dsc$num_singletons, 0 )
    expect_equal( dsc$num_doubletons, 3 * 2 )

    comp = compare_methods( Yobs ~ Z | B | D, data=blk, include_method_characteristics = FALSE )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )
    expect_true( all( !is.na( comp$SE_hat ) ) )

})






test_that( "equal sized designs", {


    blk = blkvar::generate_blocked_data( n_k = rep( 10, 18 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:18 ),
                   D = rep( 1:9, 2 ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    expect_equal( dsc$ncv, 0 )
    expect_equal( dsc$Jcv, 0 )
    expect_equal( dsc$n_block_cv, 0 )
    expect_equal( dsc$tx.cv, 0 )

    expect_true( any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )
    expect_equal( dsc$num_doubletons, 0 )

    expect_warning( expect_warning( comp <- compare_methods( Yobs ~ Z | B | D, data=blk,
                                                             include_method_characteristics = FALSE ) ) )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )





    blk = blkvar::generate_blocked_data( n_k = rep( 10, 16 ) )
    head( blk )
    blks = tibble( B = paste0( "B", 1:16 ),
                   D = rep( 1:4, 4 ),
                   Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
    blk = left_join( blk, blks, by="B" ) %>%
        mutate( Yobs = ifelse( Z, Y1, Y0 ) )

    dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
    expect_equal( dsc$ncv, 0 )
    expect_equal( dsc$Jcv, 0 )
    expect_equal( dsc$n_block_cv, 0 )
    expect_equal( dsc$tx.cv, 0 )

    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "matched pairs" ) ) )
    expect_true( !any( stringr::str_detect( attr( dsc, "notes" ), "singleton" ) ) )
    expect_equal( dsc$num_doubletons, 8 )

    comp <- compare_methods( Yobs ~ Z | B | D, data=blk,
                                                             include_method_characteristics = FALSE )
    comp
    expect_true( all( !is.na( comp$ATE_hat ) ) )

    df = comp$df
    names(df) <- comp$method
    nC = 16
    nB = 4
    df = df[ !stringr::str_detect( names(df), "MLM" ) ]

    A = nC - nB - 1
    B = nC - 2*nB
    expect_equal( as.numeric(df), c( A, B, B, A, A, B, B, B, B, B,B, A, B, B, A, B, B ) )
})



