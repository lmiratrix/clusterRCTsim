library( testthat )




test_that("formula generator works", {
    dd = data.frame( Y = 1:10, tx = rep(2,10),
                     c = sample(c("A","B"), 10, replace=TRUE),
                       SITE = sample( LETTERS[1:4], 10, replace=TRUE ),
                       x = 1:10,
                       x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd

    rr = make_regression_formula()

    expect_equal( as.character(rr), "Yobs ~ 1 + Z" )

    rr = make_regression_formula( blockID = "SITE", FE = TRUE )
    expect_equal( as.character(rr), "Yobs ~ 0 + Z + SITE" )

    rr = make_regression_formula( blockID = "SITE", interacted = TRUE )
    expect_equal( as.character(rr), "Yobs ~ 0 + Z:SITE + SITE" )

    rr = make_regression_formula( blockID = "SITE", FE = TRUE )
    expect_equal( as.character(rr), "Yobs ~ 0 + Z + SITE" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  blockID = "SITE", interacted = TRUE,
                                  control_formula = ~ x + x2 )
    expect_equal( as.character(rr), "Y ~ 0 + Z:SITE + SITE + x + x2" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  blockID = "SITE", FE = TRUE,
                                  control_formula = ~ x + x2 )
    expect_equal( as.character(rr), "Y ~ 0 + Z + SITE + x + x2" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  blockID = "SITE", FE = FALSE,
                                  control_formula = ~ x + x2 )
    expect_equal( as.character(rr), "Y ~ 1 + Z + x + x2" )


})






test_that("formula generator random effects", {

    dd = data.frame( Y = 1:10, tx = rep(2,10),
                     c = sample(c("A","B"), 10, replace=TRUE),
                     SITE = sample( LETTERS[1:4], 10, replace=TRUE ),
                     x = 1:10,
                     x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd

    rr = make_regression_formula( cluster_RE = TRUE )
    a <- as.character(rr)
    a
    expect_true( a == "Yobs ~ 1 + Z + (1 | clusterID)" )


    rr = make_regression_formula( blockID = "SITE", FE = TRUE, cluster_RE = TRUE )
    a <- as.character(rr)
    a
    expect_equal( a, "Yobs ~ 0 + Z + SITE + (1 | clusterID)" )


})







test_that("making canonical data works", {

    odat = data.frame( Y = 1:10, tx = rep(2,10), c = sample(c("A","B"), 10, replace=TRUE),
                      s = sample( c("X","Y"), 10, replace=TRUE ),
                      x = 1:10,
                      x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    odat
    odat$tx[1:5] = 1

    # dat = make_canonical_data(formula=Y ~ tx | c | s, data=odat)

    dat = clusterRCT:::make_canonical_data( formula = Y ~ tx | c | s, data=odat)

    dat
    expect_true( ncol(dat) == 4 )
    expect_equal( names(dat), c( "Yobs", "Z", "clusterID", "blockID" ) )


    dat2 = clusterRCT:::make_canonical_data(formula=Y ~ tx | c | s, data=odat,
                                            control_formula = ~ x + x2 )
    head( dat2 )
    expect_true( ncol(dat2) == 6 )

    expect_equal( names(dat2), c( "Yobs", "Z", "clusterID", "blockID", "x", "x2" ) )

})




test_that("number_controls works", {

    res = clusterRCT:::number_controls( ~ X1  + X2 + X3 )
    expect_true( res == 3 )

    res = clusterRCT:::number_controls( NULL )
    expect_true( res == 0 )

    # TODO: What if there are interaction terms in the control
    # formula??
} )




test_that("covariate interated with treatment", {

    dd = data.frame( Y = 1:10, tx = rep(2,10),
                     c = sample(c("A","B"), 10, replace=TRUE),
                     SITE = sample( LETTERS[1:4], 10, replace=TRUE ),
                     x = 1:10,
                     x2 = sample( c("A","B","C"), 10, replace=TRUE ) )
    dd


    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  blockID = "SITE", interacted = TRUE,
                                  control_formula = ~ x + x2,
                                  control_interacted = TRUE )
    expect_equal( as.character(rr), "Y ~ 0 + Z:SITE + SITE + Z * (x + x2)" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  blockID = "SITE", FE = TRUE,
                                  control_formula = ~ x + x2,
                                  control_interacted = TRUE )
    rr
    expect_equal( as.character(rr), "Y ~ 0 + Z + SITE + Z * (x + x2)" )

    rr = make_regression_formula( Yobs = "Y", clusterID = "c",
                                  blockID = "SITE", FE = FALSE,
                                  control_formula = ~ x + x2,
                                  control_interacted = TRUE )
    expect_equal( as.character(rr), "Y ~ 1 + Z + Z * (x + x2)" )


} )


