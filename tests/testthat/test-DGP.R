
# This file tests the data generating code in data_generation_code.R


library( tidyverse )
library( testthat )
library( clusterRCT )

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_crt.R" ) )


#### Test the single block DGP code ####

test_that( "single block DGP code works", {
  
  df = blkvar::generate_multilevel_data_no_cov( n.bar=200, J=50,
                                                tau.11.star = 0.1^2,
                                                ICC = 0.20,
                                                p = 0.70,
                                                variable.p = FALSE,
                                                cluster.rand = TRUE,
                                                size.impact.correlate = TRUE,
                                                finite.model = FALSE )
  head( df )
  expect_equal( length( unique( df$sid ) ), 50 )
  
  op <- capture_output( a <- describe.data(df) )
  expect_true( !is.null( op ) )
  expect_equal( nrow( a ), 2 )
  
})


test_that( "district generation function works", {
  
  df <- generate_districts(
    K = 3, 
    J.bar = 20,
    p.bar = 0.3,
    ATE.district = 0.3,
    sd.district = 0.2,
    ICC.district = 0.3,
    variable.p = TRUE,
    
    # cluster-level params
    n.bar = 30, 
    sd.school = 0.4, 
    ICC.school = 0.2,
    
    # size and correlation params
    cv.district = 1,
    cv.school = 1,
    size.impact.cor.district = TRUE,
    size.intercept.cor.district = TRUE,
    size.impact.cor.school = TRUE,
    proptx.impact.cor.district = TRUE )
  
  
  districts = attr( df, "district_key" )
  districts
  expect_equal( nrow(districts), 3 )
  
})


#### Check induce correlation code #####

test_that( "induce_cor() function", {
  
  x = 10*rnorm( 10000 )
  y = 10*rnorm( 10000 )
  yn = induce_cor( x, y, cor = 1 )
  
  cor( x, yn )
  
  expect_equal( cor( x, yn ), 1, tolerance = 0.05 )
  
  yn = induce_cor( x, y, cor = -0.5 )
  cor( x, yn )
  expect_equal( cor( x, yn ), -0.5, tolerance = 0.05 )
  
  
  
})



test_that( "generating R2 covariates works", {
  
  set.seed( 40444 )
  Y = rnorm( 1000, mean=2, sd = 0.5 )
  
  X = generate_covariate_block( Y, 0.5^2, R2 = 0 )  
  M = summary( lm( Y ~ X$X1 ) )
  expect_equal( M$r.squared, 0, tolerance = 0.01 )
  
  X = generate_covariate_block( Y, 0.5^2, R2 = 0.5 )  
  M = summary( lm( Y ~ X$X1 ) )
  expect_equal( M$r.squared, 0.5, tolerance = 0.03 )
  expect_equal( sd( X$X1 ), 1, tolerance = 0.05 )
  
  
  Y = rnorm( 1000, mean=2, sd = 3 )
  X = generate_covariate_block( Y, 3^2, R2 = 0.75 )  
  M = summary( lm( Y ~ X$X1 ) )
  expect_equal( M$r.squared, 0.75, tolerance = 0.05 )
  expect_equal( sd( X$X1 ), 1, tolerance = 0.05 )
  expect_equal( cor( X$X1, Y ), sqrt(0.75), tolerance = 0.02 )
  
  X = generate_covariate_block( Y, 0.5^2, R2 = 1 )  
  expect_equal( cor( X$X1, Y ), 1 )
  
  Y = rnorm( 3000, mean=2, sd = 6 )
  X = generate_covariate_block( Y, 6^2, R2 = 0.3, k = 4 )  
  expect_equal( sd( X$X1 ), 1, tolerance = 0.05  )
  expect_equal( sd( X$X2 ), 1, tolerance = 0.05  )
  X$Y = Y
  M = summary( lm( Y ~ ., data=X ) )
  expect_equal( M$r.squared, 0.30, tolerance = 0.05 )
  
  Y = rep( 2, 100 )
  X = generate_covariate_block( Y, 6^2, R2 = 0.3, k = 4 )
  sd( X$X1 )
  sd( X$X4 )
  expect_true( sd( X$X1 ) > 0 )
  expect_equal( ncol( X ), 4 )
  
  
  Y = runif( 1000 )
  X = generate_covariate_block( Y, 0.5^2, R2 = 0.5 )
  M = summary( lm( Y ~ X$X1 ) )
  # The above will not have right R2 since not normal data.  :-(
  
  Y = rep( 2, 1000 )
  X = generate_covariate_block( Y, 0, R2 = 0, k = 3 )
  sds <- apply( X, 2, sd )
  expect_equal( as.numeric(sds), rep(1,3), tolerance = 0.05 )
  
})



#### Check building set of districts ####

test_that( "build_disticts() works", {
  
  
  set.seed( 40404 )
  grid <- build_districts( K = 10, J.bar = 200, p.bar = 0.7,
                           ATE = 0.2, sd.ATE = 0.1,
                           n.bar = 20,
                           ICC.district = 0.2,
                           variable.p = FALSE, 
                           R2 = 0.7)
  
  
  grid
  
  expect_true( sd( grid$Xd1 ) > 0 )
  expect_true( sd( grid$J ) > 0 )
  expect_equal( grid$p * grid$J - round( grid$p * grid$J ), rep( 0, 10 ) )
  
  expect_equal( grid$p , rep( 0.7, 10 ), tolerance = 1/100 )
  
  
  
  set.seed( 430404 )
  
  grid <- build_districts( K = 10, J.bar = 100, p.bar = 0.7,
                           ATE = 0.2, sd.ATE = 0.1,
                           ICC.district = 0,
                           cv.district = 0.6,
                           variable.p = TRUE,
                           n.bar = 40,
                           num.cov.district = 5,
                           proptx.impact.cor = 1.0 )
  
  
  grid
  cor( grid$p, grid$district.tau ) 
  expect_true( cor( grid$district.tau, grid$p ) > 0.9 )
  expect_true( !is.null( grid$Xd5 ) )
  expect_equal( sd( grid$Xd5 ), 1, tolerance=0.3 )  
  
  expect_equal( grid$p * grid$J - round( grid$p * grid$J ), rep( 0, 10 ) )
  
})



#### Check limiting variation of districts ####

test_that( "limiting variation of districts works", {
  
  set.seed( 40404 )
  grid <- build_districts( K = 10, J.bar = 200, p.bar = 0.7,
                           ATE = 0.2, sd.ATE = 0.1,
                           ICC.district = 0.2,
                           n.bar = 35,
                           variable.p = FALSE,
                           cv.district = 0 )
  
  expect_true( sd( grid$J ) == 0 )
  
  set.seed( 40404 )
  grid <- build_districts( K = 10, J.bar = 200, p.bar = 0.7,
                           ATE = 0.2, sd.ATE = 0.1,
                           n.bar = 35,
                           ICC.district = 0.2,
                           variable.p = FALSE,
                           cv.district = 0.5 )
  
  expect_true( sd( grid$J ) > 0 )
  
  
  
  df <- generate_districts(
    K = 3, 
    J.bar = 20,
    p.bar = 0.2,
    ATE.district = 0.3,
    sd.district = 0,
    ICC.district = 0,
    variable.p = TRUE,
    
    # cluster-level params
    n.bar = 30, 
    sd.school = 0.4, 
    ICC.school = 0.2,

    # size and correlation params
    cv.district = 1,
    cv.school = 1,
    size.impact.cor.district = FALSE,
    size.intercept.cor.district = FALSE,
    size.impact.cor.school = TRUE,
    proptx.impact.cor.district = TRUE )
  
  
  districts = attr( df, "district_key" )
  districts
  
  expect_equal( nrow(districts), 3 )
  expect_equal( sd( districts$district.tau ), 0 )
  expect_equal( sd( districts$district.alpha ), 0 )
  
})



#### Overall variation is correct ####

test_that( "overall variation is 1", {
  
  df <- generate_districts(
    K = 25,
    J.bar = 50, 
    p.bar = 0.5,
    ATE.district = 0.2, 
    sd.district = 0.1,
    ICC.district = 0.5,
    variable.p = T,
    
    # correlation params at cluster level
    cv.district = 0.6, 
    size.impact.cor.district = 0.2,
    proptx.impact.cor.district = 0.2,
    size.intercept.cor.district = 0.2,
    
    # school-level params
    n.bar = 50, 
    sd.school = 0.1, 
    ICC.school = 0.2,
    
    # correlation params at individual level
    cv.school = 0.6,
    size.impact.cor.school = 0
  )
  
  dist <- attr( df, "district_key")
  expect_equal( nrow(dist), 25 )
  
  expect_equal( sd( df$Y0 ), 1, tolerance = 0.05 )  
  
})



#### School and individual level covariates ####

test_that( "School- and individual-level covariates generate", {
  
  df <- generate_districts(
    K = 25,
    J.bar = 50, 
    p.bar = 0.5,
    ATE.district = 0.2, 
    sd.district = 0.1,
    ICC.district = 0.5,
    variable.p = T,
    
    num.cov.district = 5,
    num.cov.school = 2,
    num.cov.indiv = 3,
    R2.district = 0.5,
    R2.school = 0.3,
    R2.indiv = 0.7,
    
    # correlation params at cluster level
    cv.district = 0.6, 
    size.impact.cor.district = 0.2,
    proptx.impact.cor.district = 0.2,
    size.intercept.cor.district = 0.2,
    
    # school-level params
    n.bar = 50, 
    sd.school = 0.1, 
    ICC.school = 0.2,
    
    # correlation params at individual level
    cv.school = 0.6,
    size.impact.cor.school = 0  )
  
  dist <- attr( df, "district_key")
  expect_equal( nrow(dist), 25 )
  head( dist )
  expect_true( all( paste0( "Xd", 1:5 ) %in% colnames(dist) ) )
  expect_true( all( paste0( "Xd", 1:5 ) %in% colnames(df ) ))
  head( df )  
  
  expect_true( all( paste0( "X", 1:3 ) %in% colnames(df) ) )
  expect_true( all( paste0( "W", 1:2 ) %in% colnames(df) ) )
  
  ICC.d = 0.5
  ICC.s = 0.3
  R2.d = 0
  R2.s = 0
  R2.i = 0.5
  
  
  df <- generate_districts(
    
    ICC.district = ICC.d,
    ICC.school = ICC.s,
    R2.district = R2.d,
    R2.school = R2.s,
    R2.indiv = R2.i,
    
    K = 25, J.bar = 50, p.bar = 0.5,
    ATE.district = 0.2, sd.district = 0.1,
    variable.p = T,
    cv.district = 0.6, 
    size.impact.cor.district = 0.2,
    proptx.impact.cor.district = 0.2,
    size.intercept.cor.district = 0.2,
    n.bar = 50, sd.school = 0.1, 
    cv.school = 0.6,
    size.impact.cor.school = 0  )
  
  dist <- attr( df, "district_key")
  expect_equal( nrow(dist), 25 )
  head( dist )
  expect_true( "Xd1" %in% colnames(dist) )
  expect_true( "Xd1" %in% colnames(df) )
  head( df )  
  
  
  expect_true( all( c( "X", "W" ) %in% colnames(df) ) )
  
} )


test_that( "School- and individual-level covariates capture ICC and R2", {
  
  set.seed( 50550 )
  ICC.d = 0.5
  ICC.s = 0.3
  R2.d = 0
  R2.s = 0
  R2.i = 0.5
  
  # Given these parameters, gen data, try to recover parameters
  one_test <- function( ICC.d, ICC.s, R2.d, R2.s, R2.i ) {
    
    rrr <- replicate( 5, 
                      { 
                        df <- generate_districts(
                          
                          ICC.district = ICC.d,
                          ICC.school = ICC.s,
                          R2.district = R2.d,
                          R2.school = R2.s,
                          R2.indiv = R2.i,
                          
                          K = 40, J.bar = 40, p.bar = 0.5,
                          ATE.district = 0.2, sd.district = 0.1,
                          variable.p = T,
                          cv.district = 0.6, 
                          size.impact.cor.district = 0.2,
                          proptx.impact.cor.district = 0.2,
                          size.intercept.cor.district = 0.2,
                          n.bar = 50, sd.school = 0.1, 
                          cv.school = 0.6,
                          size.impact.cor.school = 0 )
                        
                        expect_equal( sd( df$Y0 ), 1, tolerance = 0.15 )
                        
                        df
                      }, simplify=FALSE )
    
    Ms <- map_dbl( rrr, \( ddd ) { 
      M = lm( Y0 ~ 1 + X + W + Xd1, data=ddd )
      Ms = summary( M )
      Ms$r.squared
    } )
    Ms
    mean( Ms )
    targ_R2 <- ICC.d*R2.d + ICC.s*R2.s + (1-ICC.d-ICC.s)*R2.i
    targ_R2
    
    expect_equal( mean(Ms), targ_R2, tolerance = 0.05 )
    
    
    Ms <- map( rrr, \( df ) {
      Mmlm = lmer( Y0 ~ 1 + X + W + Xd1 + (1|sid) + (1|did), data=df )
      CI = confint(Mmlm, method="Wald")[ -(1:3), ]
      CI %>%
        as.data.frame( ) %>%
        rownames_to_column(var = "term")
    } ) 
    Ms <- bind_rows( Ms )
    gg <- Ms %>%
      group_by( term ) %>%
      summarise( low = mean( `2.5 %` ),
                 high = mean( `97.5 %` ) )
    gg
    params = tribble( ~ term, ~ param,
                      "(Intercept)", 0,
                      "W", sqrt( R2.s * ICC.s ),
                      "X", sqrt( R2.i * (1-ICC.d-ICC.s) ),
                      "Xd1", sqrt( R2.d * ICC.d ) )
    params
    CI = left_join( gg, params, by="term" )
    CI
    
    expect_true( all( CI$low <= CI$param ) )
    expect_true( all( CI$high >= CI$param ) )
  }
  
  set.seed( 40440 )
  
  one_test( 0.5, 0.3, 0, 0, 0 )
  
  # WARNING: These tests are a bit dodgy, and the CIs are not
  # capturing the truth the right probability of times.
  set.seed( 4413 )
  one_test( 0.5, 0.3, 0, 0, 0.5 )
  
  set.seed( 44134 )
  one_test( 0.6, 0.3, 0, 0.7, 0 )
  
  set.seed( 55054 )
  one_test( 0.5, 0.2, 0.4, 0, 0 )
  
  set.seed( 550223 )
  expect_warning( one_test( 0.4, 0.3, 0.9, 0.9, 0.9 ) )
  
  set.seed( 4333 )
  one_test( 0.3, 0.2, 0.99, 0.99, 0.99 )
  
})




#### School ICC = 0 testing ####


test_that( "ICC = 0 works", {
  
  set.seed( 505505 )
  df = blkvar::generate_multilevel_data_no_cov( n.bar=200, J=50,
                                                tau.11.star = 0.1^2,
                                                ICC = 0,
                                                p = 0.70,
                                                variable.p = TRUE,
                                                size.impact.correlate = TRUE,
                                                proptx.impact.correlate = TRUE,
                                                finite.model = FALSE )
  
  head( df )
  
  library( lme4 )
  M = lmer( Y0 ~ 1 + (1|sid), data=df )
  expect_true( VarCorr( M )$sid < 0.03 )
  
  
  df <- generate_districts(
    K = 25,
    J.bar = 50, 
    p.bar = 0.5,
    ATE.district = 0.2, 
    sd.district = 0.1,
    ICC.district = 0.5,
    variable.p = T,
    
    # correlation params at cluster level
    cv.district = 0.6, 
    size.impact.cor.district = 0.2,
    proptx.impact.cor.district = 0.2,
    size.intercept.cor.district = 0.2,
    
    # school-level params
    n.bar = 50, 
    sd.school = 0.1, 
    ICC.school = 0,
    
    # correlation params at individual level
    cv.school = 0.6,
    size.impact.cor.school = 0
  )
  
  dist <- attr( df, "district_key")
  expect_equal( nrow(dist), 25 )
  
  library( lme4 )
  M = lmer( Y0 ~ 1 + (1|sid) + (1|did), data=df )
  if ( FALSE ) {
    arm::display( M )
    VarCorr( M )$sid
    VarCorr( M )$did
  }
  expect_true( VarCorr( M )$sid < 0.03 )
  expect_true( VarCorr( M )$did > 0.4 )
  
  
} )


