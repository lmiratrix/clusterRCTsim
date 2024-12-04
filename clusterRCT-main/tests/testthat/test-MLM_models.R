
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



test_that("MLM estimation works", {

    md = MLM_estimators( Yobs ~ T.x | S.id | D.id, data = sim.data )
    md

    expect_true( is.data.frame( md ) )
    expect_true( all( md$p_value <= 1 ) )
    expect_true( all( md$SE_hat > 0 ) )

    md = MLM_estimators( Yobs ~ T.x | S.id, data = sim.data )
    expect_true( md$method == "MLM" )

    sim.data$Yobs[ sim.data$T.x == 1 ] = sim.data$Yobs[ sim.data$T.x == 1 ] + 0.3
    md = MLM_estimators( Yobs ~ T.x | S.id | D.id, data = sim.data )
    md
    expect_true( all( md$p_value <= 0.05 ) )

})


test_that( "get degrees of freedom hack right", {
    m2 = model.params.list
    m2$K = 3
    m2$J = 4

    ss <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", m2, Tbar = 0.5 )
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


test_that( "warning suppression works", {

    # TODO: Neither of these throws a warning???

    data( fakeCRT )
    mtab_cov = NA

    mtab_cov <-  MLM_estimators( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                 control_formula = ~ X.jk + C.ijk )

    mtab_cov2 = NA
    mtab_cov2 <- MLM_estimators( Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                 control_formula = ~ X.jk + C.ijk,
                                 suppress_warnings = FALSE )

    expect_equal( mtab_cov, mtab_cov2 )
})
