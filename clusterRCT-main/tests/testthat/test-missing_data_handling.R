
if ( FALSE ) {

    f <- function(x) {
        if ( x < -1 ) {
            warning( "REALLY negative" )
        }
        if (x < 0) {
            warning("*x* is already negative")
            return(x)
        }
        -x
    }

    f(-2)
    f2 <- function(x) {
        f( x * 5 )
    }
    expect_warning(f2(-2))
    expect_warning(f2(-1))
    f2(-1)
    expect_warning(f(-1), "already negative")
    expect_warning(f(1), NA)
    f(-1)
    f(1)
}




test_that("check missing data doesn't crash", {

    data( fakeCRT )
    fakeCRT$Yobs[1:10] = NA

    fakeCRT$Yobs2 = fakeCRT$Yobs + rnorm( nrow(fakeCRT ) )
    fakeCRT$Yobs2[5:20] = NA

    expect_warning( cm3 <- describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=fakeCRT ) )

    cm <- compare_methods( Yobs2 + Yobs ~ T.x | S.id | D.id, data=fakeCRT, include_MLM = FALSE, warn_missing = FALSE )


    fakeCRT$T.x[5:20] = NA
    expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                     include_MLM = FALSE, warn_missing = TRUE ) )

    expect_warning( describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=fakeCRT ) )

    data(fakeCRT)
    fakeCRT$X.jk[30:40] = NA
    expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                     include_MLM = FALSE,
                                     control_formula = ~ X.jk ) )

    data(fakeCRT)
    fakeCRT$S.id[ 100:150 ] = NA
    expect_warning( describe_clusterRCT( Yobs ~ T.x | S.id | D.id,
                                         data=fakeCRT, control_formula = ~ X.jk ) )

    fakeCRT$D.id[ 70:150 ] = NA


    expect_warning( clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id,
                                                      data=fakeCRT, control_formula = ~ X.jk ) )

})


test_that( "patching and missing data works", {

    data( "fakeBrokeCRT" )

    table( Z = fakeBrokeCRT$T.x, Did = fakeBrokeCRT$D.id, useNA = "always" )

    p = patch_singleton_blocks( Yobs ~ T.x | S.id | D.id, data = fakeBrokeCRT,
                                warn_missing = FALSE )
    head(p)
    table( p$T.x, p$D.id )
    pA = patch_data_set( Yobs ~ T.x | S.id | D.id, data=p )
    nrow( pA )
    tb <- table( pA$Z, pA$blockID )
    tb
    expect_true( sum( tb == 0 ) > 0 )

    head( pA )
    table( pA$blockID, pA$Z )

    p2 = clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data = fakeBrokeCRT,
                                           warn_missing = FALSE, patch_data = TRUE )
    head( p2 )
    p3 = patch_singleton_blocks( Yobs ~ Z | clusterID | blockID, data=p2,
                                 warn_missing = FALSE )
    nrow( p3 )


    #expect_equal( nrow(pA), nrow(p3) )
    expect_true( nrow(pA) > nrow(p3) )
})


