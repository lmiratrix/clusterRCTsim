

test_that("check data integrity works", {

    data( fakeCRT )
    head( fakeCRT )

    expect_true( is_nested( fakeCRT$S.id, fakeCRT$D.id ) )
    expect_true( !is_nested( fakeCRT$D.id, fakeCRT$S.id ) )

    dat1 = clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    rs <- clusterRCT:::check_data_integrity( dat1 )
    expect_true( rs )

    dat1$Z[ 1:3 ] = 1
    expect_error( clusterRCT:::check_data_integrity( dat1 ) )

    dat1$Z[1:3] = 2
    expect_error( clusterRCT:::check_data_integrity( dat1 ) )

    d2 = dat1
    d2$Z[ d2$SiteID == d2$SiteID[[1]] ] = 0
    expect_error( clusterRCT:::check_data_integrity( d2 ) )

    dat1 = dat1[ -c(1:3), ]

    dat1$clusterID[1] = dat1$clusterID[2]
    expect_error( clusterRCT:::check_data_integrity( dat1 ) )
    dat1 = dat1[ -c(1:3), ]

    dat2 = clusterRCT:::make_canonical_data( Yobs ~ T.x | D.id | S.id, data=fakeCRT )
    expect_error( clusterRCT:::check_data_integrity( dat2 ) )
})






test_that( "error messages look good", {

    data( fakeCRT )
    #head( fakeCRT )

    # The following error messages should all list the missing variables.
    expect_error( compare_methods( Yobs ~ Z | S.id, data = fakeCRT ) )
    expect_error( compare_methods( Yobs2 ~ Z | S.id, data = fakeCRT ) )
    expect_error( compare_methods( Yobs ~ Z | S.id | boo, data = fakeCRT,
                                   control_formula = ~ doggy ) )
    expect_error( compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                                   control_formula = ~ doggy ) )
    cc <- compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                           patch_data = FALSE,
                           control_formula = ~ X )
    expect_true( is.data.frame(cc) )
})






test_that( "combining blocks works", {

    data(fakeCRT)
    a = identify_singleton_blocks( Yobs ~ T.x |  S.id | D.id, data=fakeCRT  )
    expect_true( is.null( a ) )

    fakeCRT$T.x[ fakeCRT$D.id == 2 ] = 1
    fakeCRT$T.x[ fakeCRT$D.id == 9 ] = 0
    fakeCRT <- rename( fakeCRT,
                       alt.id = D.id )
    expect_error( compare_methods( Yobs ~ T.x | S.id | alt.id, data=fakeCRT,
                                   handle_singleton_blocks = "fail" ) )


    # Now try to identify bad blocks
    a = identify_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=fakeCRT  )
    expect_true( is.data.frame(a) )
    expect_equal( a$pTx, c(1,0) )

    # And now patch!
    ptched = patch_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=fakeCRT,
                                drop_data = FALSE, warn_missing = FALSE )
    tb = table( ptched$T.x, ptched$alt.id )
    expect_true( all( tb > 0 ) )

    aa = ptched %>% group_by( alt.id ) %>%
        summarise( n = n(),
                   nu = length( unique( S.id ) ) )
    expect_true( aa$nu[[1]] == 2 )
    cnt = sum( fakeCRT$alt.id %in% c( "2", "9" ) )
    expect_true( aa$n[[1]] == cnt)


    # don't pool clusters
    a = patch_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=fakeCRT,
                                drop_data = FALSE,
                                pool_clusters = FALSE, warn_missing = FALSE  )

    n_c = length( unique( ptched$S.id ) )
    n_c2 = length( unique( a$S.id ) )
    expect_true( n_c < n_c2 )
    expect_equal( ptched$alt.id, a$alt.id )

    tb = table( a$T.x, a$alt.id )
    expect_true( all( tb > 0 ) )

    aa = a %>% group_by( alt.id ) %>%
        summarise( n = n(),
                   nu = length( unique( S.id ) ) )
    expect_true( aa$nu[[1]] > 2 )
    expect_true( aa$n[[1]] == cnt)

    # test dropping
    a_drp = patch_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=fakeCRT,
                                drop_data = TRUE, warn_missing = FALSE )
    expect_true( nrow( fakeCRT ) > nrow(a_drp) )

    ds <- identify_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=a_drp )
    expect_true( is.null( ds ) )

    # If only a single singeton block?
    ff = dplyr::filter( fakeCRT, alt.id != "2" )
    a = identify_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=ff  )
    a
    expect_equal( nrow(a), 1 )
    expect_true( a$pTx == 0 )

    expect_warning( b <- patch_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=ff,
                                drop_data = FALSE, warn_missing = FALSE ) )

    dd <- identify_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=b )
    expect_equal( dd$blockID, ".pooled" )

    c <- patch_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=ff,
                                                 drop_data = TRUE, warn_missing = FALSE )
    c_id <- identify_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=c )
    expect_true( is.null( c_id ) )


    # Missing block ID?

    # Make little dataset with missing block ID
    f2 = fakeCRT
    f2$alt.id[ f2$alt.id == "9" ] = NA
    a = identify_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=f2  )
    a
    expect_equal( nrow(a), 2 )

    # now pool
    b_null = patch_singleton_blocks( Yobs ~ T.x | S.id | alt.id, data=f2,
                                drop_data = FALSE, warn_missing = FALSE )
    expect_equal( nrow(b_null), nrow(ptched) )
    tb = table( b_null$T.x, b_null$alt.id )
    tb2 = table( ptched$T.x, ptched$alt.id )
    expect_equal( tb, tb2 )
})
