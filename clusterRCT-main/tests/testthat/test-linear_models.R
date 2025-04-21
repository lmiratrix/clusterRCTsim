

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


set.seed( 404404 )

sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )




test_that("linear models works", {

    head( sim.data )
    nrow( sim.data )
    J = length( unique( sim.data$S.id ) )
    lme <- linear_model_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    lme
    expect_true( nrow( lme) == 1 )
    expect_equal( lme$df, J - 10 - 1 )

    head( sim.data )
    lme_cov <- linear_model_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                        control_formula = ~ X.jk + C.ijk )
    lme_cov
    expect_true( lme$ATE_hat != lme_cov$ATE_hat )
    expect_true( lme_cov$df <= lme$df - 3 )
    #expect_equal( lme_cov$df, J - 10 - 1 - 3 )


    agg = aggregation_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    agg
    expect_true( !is.null( agg ) )

})


test_that("interacted linear models works", {

    data( fakeCRT )
    fakeCRT

    J = length( unique( fakeCRT$S.id ) )
    J
    e2 <- interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                              control_formula = ~ X.jk + C.ijk,
                                              data = fakeCRT )
    expect_true( is.data.frame(e2) )
    expect_true( !is.null( e2$method ) )

    # Ensure some sort of df is returned.
    expect_true( all( e2$df > 30 ) )

    e1 <- interacted_linear_model_estimators( Yobs ~ T.x | S.id | D.id,
                                              data = fakeCRT )
    expect_true( nrow(e1) == 3 )

    # No blocks mean no estimators from this. Get a 0 row table.
    e0 <- interacted_linear_model_estimators( Yobs ~ T.x | S.id,
                                              data = fakeCRT )
    expect_true( nrow( e0 ) == 0 )
})
