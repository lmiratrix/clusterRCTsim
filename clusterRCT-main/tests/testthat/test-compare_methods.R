

if ( FALSE ) {
    library( testthat )

}

data( "fakeCRT" )

test_that("compare methods aggregates as expected", {

    mtab <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                             include_method_characteristics = FALSE )
    mtab
    expect_true( is.data.frame(mtab) )

    mtab_cov <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                  control_formula = ~ X.jk + C.ijk,
                                  include_method_characteristics = FALSE )

    # Ensure covariate correction changes the estimates!  (I.e.,
    # control formula is getting propagated)
    mtab$dels = mtab$ATE_hat - mtab_cov$ATE_hat
    mtab
    expect_true( sum( mtab$dels == 0 ) == 2 )

    mtab <- compare_methods( Yobs ~ T.x | S.id, data=fakeCRT,
                             include_method_characteristics = FALSE)

    expect_true( is.data.frame(mtab) )
    mtab

})


test_that("categorical covariates handled", {

    data( fakeCRT )
    cc <- compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                           patch_data = FALSE, include_MLM = FALSE,
                           control_formula = ~ X )
    expect_true( is.data.frame(cc) )

    cc <- compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                           patch_data = TRUE, include_MLM = FALSE,
                           control_formula = ~ X )
    expect_true( is.data.frame(cc) )

})



test_that("missing data handled as desired", {
    data( fakeCRT )
    fake2 = fakeCRT
    head( fake2 )
    fake2$Yobs[1:10] = NA
    fake2$S.id[5:20] = NA
    fake2$X.jk[(1:30)*2] = NA
    fake2$T.x[c(1, 5, 17, 61)] = NA
    fake2$D.id[c( 5, 66 )] = NA

    nrow(fake2)

    expect_warning( mtab <- compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                             patch_data = FALSE,
                                             include_method_characteristics = FALSE) )

    expect_true( is.data.frame(mtab) )

    expect_warning( mtab_cov <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                                  control_formula = ~ X.jk + C.ijk,
                                                  patch_data = FALSE,
                                                  include_method_characteristics = FALSE ) )
    expect_true( is.data.frame(mtab_cov) )

    mtab_cov

    mtab_cov2 <-  compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                   control_formula = ~ X.jk + C.ijk,
                                   include_method_characteristics = FALSE,
                                   warn_missing = FALSE )

    expect_true( is.data.frame(mtab_cov2) )

    # Kept more data, lower standard errors!
    expect_true( mean( mtab_cov2$SE_hat ) < mean( mtab_cov$SE_hat ) )

    # Check that if missingness drops all tx in a district, we stop.
    fake2$T.x[ is.na( fake2$D.id ) | is.na( fake2$T.x ) | (fake2$D.id == 1 & fake2$T.x == 0) ] = NA

    expect_warning( tt <- make_block_table( Yobs ~ T.x | S.id | D.id, data=fake2,
                                            control_formula = ~ X.jk + C.ijk ) )
    expect_true( tt$p.tx[[1]] == 1 )

    # Will throw warning about dropping all-tx blocks when trying to analyze.
    expect_warning( expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                                     control_formula = ~ X.jk + C.ijk,
                                                     patch_data = FALSE,
                                                     include_method_characteristics = FALSE ) ) )

    # Patching data will not help with all tx or all co blocks.
    expect_warning( expect_warning( compare_methods( Yobs ~ T.x | S.id | D.id, data=fake2,
                                                     control_formula = ~ X.jk + C.ijk,
                                                     patch_data = TRUE, warn_missing = TRUE,
                                                     include_method_characteristics = FALSE ) ) )
    #expect_true( is.data.frame(mtab_cov) )
})


