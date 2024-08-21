
library( testthat)

test_that("patch data works", {

    data <- tribble( ~ X, ~ R, ~ Q, ~ W,
                     1,  2, 3, 4,
                     NA, 2, 4, 5,
                     NA, NA,4, 6,
                     2,  3, NA, 7)

    data$T.x = c( 1, 1, 0, 0 )
    data$S.id = c( 10, 5, 10, 5 )
    data = bind_rows( data, data, data, data, data)
    data$Y = rnorm( nrow(data), sd=10 )
    data$X[1:5] = 1:5
    data$W[5:10] = 4
    data$RR = sample( c("A","B","C"), nrow(data), replace=TRUE )
    data$RR[17] = NA
    data$D.id = rep( c("A","B"), c(8, 12) )

    data
    n = nrow(data)
    n
    ptch <- patch_data_set( Y ~ T.x | S.id | D.id, data=data,
                            control_formula = ~ X + R + Q + W + RR )
    ptch
    expect_true( nrow(ptch) == n )

    mn = mean( data$RR == "C", na.rm=TRUE )
    expect_equal( ptch$RRC[17], mn, tolerance = 0.000001 )

    mod = lm( Yobs ~ X_mi + R_mi + Q_mi + RRB_mi,
              data=ptch)
    coef(mod)
    expect_true( all( !is.na( coef(mod) )) )

    #mod = lm( Yobs ~ Z + X + R + Q + W + RRB + RRC + X_mi + R_mi + Q_mi + RRB_mi,
    #          data=ptch)
    #coef(mod)



    data$T.x[9] = NA
    data$Y[12] = NA
    data$D.id[11] = NA
    data$S.id[2] = NA

    ptch2 <- patch_data_set( Y ~ T.x | S.id | D.id, data=data,
                            control_formula = ~ X + R + Q + W + RR )
    ptch2
    expect_true( all( ptch2$Z %in% c(0, 1) ) )


    # Keep rows with missing D.id if we don't use it in the formula.
    ptch3 <- patch_data_set( Y ~ T.x | S.id, data=data,
                             control_formula = ~ X + R + Q + W + RR )
    expect_true( nrow( ptch3 ) > nrow(ptch2) )

    ptch4 <- patch_data_set( Y ~ T.x | S.id, data=data, warn_missing = FALSE )
    ptch4
    expect_true( nrow( ptch3 ) == nrow(ptch4) )
    expect_true( ncol( ptch4 ) == 3 )


    data
    dataC = clusterRCT:::make_canonical_data(Y ~ T.x | S.id | D.id, data=data,
                                            control_formula = ~ X + R + Q + W + RR,
                                            drop_missing = FALSE, warn_missing = FALSE )

    dataC
    ptch5 = patch_data_set( NULL, data=dataC,
                            control_formula = ~ X + R + Q + W + RR )
    expect_equal( ptch5, ptch2 )



    # Testing control formula is updated
    a5 = attr( ptch5, "control_formula" )
    expect_equal( rhs.vars(a5),
                  c( "X", "R", "Q", "W", "RRB", "RRC",
                     "X_mi", "R_mi", "Q_mi", "RRB_mi" ) )



    # Testing canonical data version works
    head( ptch5 )
    nrow(ptch5)
    pp = patch_data_set( data=ptch5,
                         control_formula = ~ X + R + Q )
    expect_equal( nrow(ptch5), nrow(pp) )

    ptch5$Q[1:5] = NA
    ptch5$Z[1:2] = NA
    ptch5$Q_mi = NULL
    pp = patch_data_set( data=ptch5,
                         control_formula = ~ X + R + Q )
    expect_equal( nrow(pp), nrow(ptch5) - 2 )

})



test_that("patch data works with no missing data", {

    data( fakeCRT )

    cf = ~ V.k + X.jk
    ptch = patch_data_set( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                           control_formula = cf )
    can = clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                            control_formula = cf )
    cf2 = attr( ptch, "control_formula" )
    expect_equal( rhs.vars(cf), rhs.vars(cf2) )

    expect_true( all( ptch == can ) )

} )



test_that( "tiny examples", {
    # demo of the patch method itself
    data( "fakeBrokeCRT" )
    dset = fakeBrokeCRT[11:30,] %>%
        dplyr::select(  -V.k, -X.jk )
    dset$C.ijk[1:3] = NA
    dset$X[19:20] = NA
    dset

    expect_warning( patch_data_set( Yobs ~ T.x | S.id | D.id, data=dset,
                    warn_missing = TRUE ) )
    expect_warning( pp <- patch_data_set( Yobs ~ T.x | S.id | D.id, data=dset,
                    control_formula = ~ C.ijk + X,
                    warn_missing = TRUE ) )
    expect_equal( nrow( pp ), 15 )

    # from colin
    data( fakeCRT )
    testdata <- fakeCRT %>%

        mutate(Yobs_miss = ifelse(S.id == 10, NA, Yobs))



    cc <- compare_methods( Yobs_miss ~ T.x | S.id | D.id, data=testdata,

                     include_method_characteristics = FALSE, warn_missing = FALSE )
    expect_true( is.data.frame(cc) )

})
