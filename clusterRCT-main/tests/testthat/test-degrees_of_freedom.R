
# Degrees of freedom verification



model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 4                          # number of schools
    , K = 3                          # number of districts
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


set.seed( 4044024 )
ss <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )




test_that( "degrees of freedom calculated right (compare_method)", {

    table( ss$D.id, ss$S.id)
    md = compare_methods( Yobs ~ T.x | S.id | D.id, data = ss ) %>%
        dplyr::select( method, df )
    md

    head( ss )
    md2 = compare_methods( Yobs ~ T.x | S.id | D.id, data = ss,
                           control_formula = ~ X.jk + C.ijk )
    md$twoCov = md2$df

    ss$X = sample( LETTERS[1:5], nrow(ss), replace=TRUE )
    md3 = linear_model_estimators( Yobs ~ T.x | S.id | D.id, data = ss,
                           control_formula = ~ X.jk + C.ijk + X )

    # Throws a warning
    #md3 = aggregation_estimators( Yobs ~ T.x | S.id | D.id, data = ss,
    #                               control_formula = ~ X.jk + C.ijk + X )

    expect_warning( expect_warning( md3 <- compare_methods( Yobs ~ T.x | S.id | D.id, data = ss,
                           control_formula = ~ X.jk + C.ijk + X ) ) )
    md$catCov = md3$df

    md
    expect_true( all( !is.na( md$twoCov ) ) )
    #expect_true( all( md$twoCov <= md$df ) )

    # TODO: What should this table look like?  Is it correct?
})



test_that( "get degrees of freedom hack right (MLM)", {

    table( ss$D.id, ss$S.id)
    md = MLM_estimators( Yobs ~ T.x | S.id | D.id, data = ss )
    expect_equal( md$df[3:5], c( 6,6,6) )

    head( ss )
    md = MLM_estimators( Yobs ~ T.x | S.id | D.id, data = ss,
                         control_formula = ~ X.jk + C.ijk )
    expect_equal( md$df[3:5], c( 4,4,4) )

    ss$X = sample( LETTERS[1:5], nrow(ss), replace=TRUE )
    md = MLM_estimators( Yobs ~ T.x | S.id | D.id, data = ss,
                         control_formula = ~ X.jk + C.ijk + X )
    md
    expect_equal( md$df[3:5], c( 0,0,0) )
})

