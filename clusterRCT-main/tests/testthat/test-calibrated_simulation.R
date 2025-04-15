

test_that("calibrated simulation works and summary methods work", {

    data( fakeCRT )

    rs <- run_rerandomize_simulation( Yobs ~ T.x | S.id | D.id,
                                      data = fakeCRT,
                                      control_formula = ~ X.jk + C.ijk,
                                      R = 3,
                                      include_MLM = FALSE,
                                      summarize_results = FALSE )

    expect_true( is.data.frame(rs) )
    tt = table( rs$runID )
    table(rs$weight)
    expect_equal( length(tt), 4 )
    tt

    rss = summarize_simulation_results( rs, summarize_results = "method" )
    rss
    expect_true( is.data.frame(rss) )
    expect_true( "method" %in% names(rss) )

    rs2 = summarize_simulation_results( rs, summarize_results = "cross" )
    rs2
    expect_true( "runID" %in% names(rs2) )

    rs3 = summarize_simulation_results( rs, summarize_results = "cross-agg" )
    rs3
    expect_equal( nrow(rs3), 3 )
    expect_equal( rs3$group, c( "Cluster", "Person", "Overall" ) )



    # Missing data stuff
    data( "fakeBrokeCRT" )
    fakeBrokeCRT <- fakeBrokeCRT %>%
        rename( YYY = Yobs,
                TTT = T.x )
    expect_warning( expect_warning(
        rs <- run_rerandomize_simulation( YYY ~ TTT | S.id | D.id,
                                      data = fakeBrokeCRT,
                                      control_formula = ~ X.jk + C.ijk + X,
                                      R = 1,
                                      summarize_results = FALSE )
    ))

    # No warnings check, also simple blocking
    rs <- run_rerandomize_simulation( YYY ~ TTT | S.id,
                                      data = fakeBrokeCRT,
                                      control_formula = ~ X.jk + C.ijk + X,
                                      R = 1,
                                      warn_missing = FALSE,
                                      summarize_results = "cross" )
    expect_equal( nrow( rs ), 6 )
})




