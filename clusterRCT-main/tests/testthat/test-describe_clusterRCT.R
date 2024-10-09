

library( tidyverse )
library( clusterRCT )

set.seed( 1039 )






test_that( "describer works", {


    data( fakeCRT )

    d <- describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                              control_formula = ~  X.jk + C.ijk )

    d
    expect_true( is.clusterRCTstats(d) )

    dd = as.data.frame(d)
    expect_true( nrow( dd ) == 1 )

    # We can check nesting of School in District
    expect_true( is_nested( fakeCRT$S.id, fakeCRT$D.id ) )
    expect_true( !is_nested( fakeCRT$D.id, fakeCRT$S.id ) )

    # We can calculate statistics for each district
    st <- make_block_table( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                            control_formula = ~  X.jk + C.ijk )
    expect_true( nrow( st ) == 10 )


    # It will trap errors!
    fakeCRT$T.x[1:4] = 2
    expect_error( describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                       control_formula = ~  X.jk + C.ijk ) )


    fakeCRT$T.x[1:4] = NA
    expect_warning( d2 <- describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT ) )

    bb = bind_rows( d, d2 )
    bb

    expect_true( nrow(bb) == 2 )
} )



test_that( "describer works with non-blocked data", {
    data( fakeCRT )

    fakeCRT$bb = fakeCRT$D.id
    d <- describe_clusterRCT( formula = Yobs ~ T.x | S.id, data=fakeCRT,
                              control_formula = ~  X.jk + C.ijk + bb )

    d
    as.data.frame(d)
    expect_true( is.clusterRCTstats(d) )
    expect_equal( d$ncov.1, 1 )
    expect_equal( d$ncov.2, 1 + length( unique( fakeCRT$D.id ) ) )

    d2 <- describe_clusterRCT( formula = Yobs ~ T.x | S.id, data=fakeCRT,
                              control_formula = ~  X.jk + C.ijk )

    expect_equal( d2$ncov.1, d$ncov.1 )
    expect_equal( d2$R2.1, d$R2.1 )
    as.data.frame(d2)
})



test_that( "describer with multiple outcomes works", {


    data( fakeCRT )
    nrow( fakeCRT )
    fakeCRT$Y2 = fakeCRT$Yobs + rnorm( nrow(fakeCRT) )
    fakeCRT$Y2[1:200] = NA
    fakeCRT$X.jk2 = fakeCRT$X.jk + fakeCRT$Yobs + rnorm( nrow(fakeCRT), sd=0.2 )

    expect_warning( d <- describe_clusterRCT( formula = Yobs + Y2 ~ T.x | S.id | D.id, data=fakeCRT,
                                              control_formula = ~  X.jk + C.ijk + X.jk2 ) )

    d
    expect_true( is.data.frame(d) )
    expect_true( nrow(d) == 2 )

} )



test_that( "more general describer tests", {

    set.seed( 1039 )

    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of clusters
        , K = 10                          # number of blocks
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of block assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of block variation
        , ICC.3 = 0.2             # block intraclass correlation
        , omega.3 = 0.2           # ratio of block effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )

    data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

    data = slice_sample( data, n = nrow(data) / 2 )
    dd = sample( unique(data$S.id), 10 )
    data = filter( data, !( S.id %in% dd ) )
    data <- clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=data )
    head( data )

    dd <- describe_clusterRCT( data=data )
    dd
    expect_true( is.clusterRCTstats( dd ) )
    expect_true( is.data.frame( dd ) )


    # Check R2 calculations
    head( fakeCRT )
    dsc = describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                               control_formula = ~  X.jk + C.ijk )

    expect_true( dsc$R2.1 > 0 )
    expect_true( dsc$R2.2 > 0 )





})


test_that( "ICC calcs work", {

    set.seed( 4054053)
    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 40                          # number of schools
        , K = 60                          # number of districts
        , nbar = 20                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of district variation
        , ICC.3 = 0.5             # district intraclass correlation
        , omega.3 = 0.2           # ratio of district effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )

    data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

    head( data )

    dd <- describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=data,
                               control_formula = ~ X.jk )
    dd
    expect_equal( dd$block_ICC, 0.5, tolerance = 0.1 )
    expect_equal( dd$cluster_ICC, 0.2, tolerance = 0.1 )
    expect_equal( dd$R2.2, 0.1, tolerance = 0.1 )

    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 90                          # number of schools
        , K = 40                          # number of districts
        , nbar = 20                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.0              # percent of district variation
        , ICC.3 = 0             # district intraclass correlation
        , omega.3 = 0.2           # ratio of district effect size variability
        , R2.2 = 0.0              # percent of school variation
        , ICC.2 = 0.05             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.8    # percent of indiv variation explained
    )

    data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

    head( data )


    dd <- describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=data,
                               control_formula = ~ C.ijk )
    dd
    expect_equal( dd$block_ICC, 0.05, tolerance = 0.1 )
    expect_equal( dd$cluster_ICC, 0, tolerance = 0.1 )
    expect_equal( dd$R2.1, 0.8, tolerance = 0.1 )
})



test_that( "R2 calcs work", {

    dd = data.frame( cid = rep( 1:10, 1:10 ) )
    dd = mutate( dd,
                 Ybar = cid,
                 Z = as.numeric( cid %% 3 == 1 ) )
    dd2 = bind_rows( dd, dd )
    dd$cid = dd$cid + 100
    dd3 = bind_rows( dd, dd2 )
    dd$cid = dd$cid - 100
    dd2$Ybar = dd2$Ybar + 1
    dd3$Ybar = dd3$Ybar + 2
    dd = bind_rows( A = dd, B = dd2, C = dd3, .id = "sid" )
    dd$Y = dd$Ybar + rnorm( nrow(dd) )
    dd <- mutate( dd,
                  Z = ifelse( cid <= as.numeric(as.factor(sid))*2, 1, 0 ) )
    head(dd)

    bt <- make_block_table(  Y ~ Z | cid | sid, data=dd )
    expect_true( nrow( bt ) == 3 )

    desc <- describe_clusterRCT( Y ~ Z | cid | sid, data=dd)
    dd <- dd %>%
        group_by( sid ) %>%
        mutate(   X1 = Y + rnorm(n(), sd=3),
                  X2 = cid,
                  X3 = mean(Y),
                  X4 = sample( c( "A", "B", "C"), n(), replace=TRUE ),
                  Y = Y + 1 * Z ) %>%
        ungroup()

    data = clusterRCT:::make_canonical_data( Y ~ Z | cid | sid, data=dd,
                                             control_formula = ~ X1 + X2 + X3 + X4 )
    head( data )
    expect_true( ncol(data) == 8 )

    a <- clusterRCT:::calc_covariate_R2s( data )
    a
    expect_true( a$ncov.1 == 3 )
    expect_true( a$ncov.2 == 5 )


    set.seed( 445040 )
    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 50                          # number of schools
        , K = 5                          # number of districts
        , nbar = 30                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.5              # percent of district variation
        , ICC.3 = 0.2             # district intraclass correlation
        , omega.3 = 0.2           # ratio of district effect size variability
        , R2.2 = 0.75              # percent of school variation
        , ICC.2 = 0.4             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.25    # percent of indiv variation explained
    )
    sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
    kp = rbinom( nrow(sim.data), 1, prob = (50+as.numeric(sim.data$S.id)) / 300 )
    mean(kp)
    head( sim.data )
    sim.data = sim.data[ kp == 1, ]
    cdat = clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=sim.data, control_formula = ~ C.ijk + X.jk )
    desc = clusterRCT:::calc_covariate_R2s( cdat )
    expect_equal( desc$R2.1, 0.25, tolerance = 0.2 )
    expect_equal( desc$R2.2, 0.75, tolerance = 0.2 )

})






test_that( "Missing data handled", {


    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of clusters
        , K = 10                          # number of blocks
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of block assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of block variation
        , ICC.3 = 0.2             # block intraclass correlation
        , omega.3 = 0.2           # ratio of block effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )

    set.seed( 1039 )
    data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
    head( data )
    n = nrow(data)
    n

    ms = sample( n, 1500 )
    data$C.ijk[ms] = NA
    data$T.x[1:200] = NA
    data$Yobs[ms] = data$Yobs[ms] + 20

    dd = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=data,
                              warn_missing = FALSE )
    dd
    expect_equal( dd$n, n - 200 )

    dd = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=data,
                              control_formula = ~ X.jk,
                              warn_missing = FALSE )
    dd
    expect_equal( dd$n, n - 200 )

    as.data.frame(dd)

    dd2 = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=data,
                               control_formula = ~ C.ijk + X.jk,
                               warn_missing = FALSE )
    dd2
    as.data.frame( dd2 )
    expect_true( dd2$R2.1 > 0.5 )

    # Extra covariate for missing indicator
    expect_equal( dd2$ncov.1, 2 )

})




#test_that( "singleton blocks and matched pairs works", {

    # See test_unusual_designs.R

#})

