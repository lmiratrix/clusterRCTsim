

test_that("make canonical data works", {

    data( fakeCRT )
    fakeCRT$X[1:10] = NA

    a = clusterRCT:::make_canonical_data(Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    head( a )
    expect_true( ncol( a ) == 4 )

    a = clusterRCT:::make_canonical_data(Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                         control_formula = ~ X, warn_missing = FALSE)
    head( a )
    expect_true( nrow( a ) == nrow(fakeCRT) - 10 )

    a = clusterRCT:::make_canonical_data(Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                         control_formula = ~ X, drop_missing = FALSE, warn_missing = FALSE )
    head( a )
    expect_true( nrow( a ) == nrow(fakeCRT)  )


    fakeCRT$T.x[101:120] = NA
    a = clusterRCT:::make_canonical_data(Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                         control_formula = ~ X, drop_missing = FALSE, warn_missing = FALSE )
    expect_true( nrow( a ) == nrow(fakeCRT)  )

    a = clusterRCT:::make_canonical_data(Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                                         control_formula = ~ X, drop_missing = TRUE, warn_missing = FALSE )
    sum( a$Z )
    expect_equal( nrow( a ), nrow(fakeCRT) - 30 )
})
