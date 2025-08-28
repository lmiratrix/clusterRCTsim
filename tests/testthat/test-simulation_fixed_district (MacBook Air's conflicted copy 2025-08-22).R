


# Test the simulation support functions
library( testthat )
library( tidyverse )

rm( list=ls() )

source( here::here( "R/simulation_support_code_fixed_district.R" ) )


##### Test the DGP code #####

test_that( "single.CRT.fixed.district() works", {

  grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.ATE = 0.1,
                           ICC.district = 0.2,
                           cv.district = 0.2,
                           variable.p = TRUE, cv.nbar = 0.2,
                           num.cov.school = 3, size.intercept.cor = -0.8,
                           size.impact.cor = 0.75, proptx.impact.cor = 0.75, nbar.J.cor = 0.9, 
                           empirical = TRUE )
  
  grid
  describe_district_grid( grid )
  
  grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.ATE = 0.1,
                           ICC.district = 0.2,
                           cv.district = 0.2,
                           variable.p = TRUE, cv.nbar = 0.2,
                           num.cov.district = 3, size.intercept.cor = -0.8,
                           size.impact.cor = 0.75, proptx.impact.cor = 0.75, nbar.J.cor = 0.9, 
                           empirical = TRUE )
  
  grid
  describe_district_grid( grid )
  
  rr = single.CRT.fixed.district( grid = grid, n.runs = 3 )
  rr
  expect_true( is.data.frame(rr) )
  expect_equal( n_distinct(rr$subrun), 3 )
  
  head( rr )
  
  expect_equal( sd( rr$ATE.finite.dist.school ), 0 )
  expect_equal( sd( rr$ATE.finite.school ), 0 )
  expect_equal( sd( rr$ATE.finite.student ), 0 )
  expect_equal( sd( rr$ATE.finite.dist.student ), 0 )
  
  g <- rr %>% group_by( method ) %>%
    summarize( sdATE = sd(ATE.hat),
               n = n() )
  expect_true( all( g$n == 3 ) )
  expect_true( all( g$sdATE > 0.001 ) )
  
  
  rr2 = single.CRT.fixed.district( grid, just.describe.data = TRUE )
  expect_equal( rr2$ATE.distr[[1]], 0.2 )
  
  rr2 = single.CRT.fixed.district( grid, just.return.data = TRUE )
  expect_true( is.data.frame( rr2 ) )
  expect_true( all( c( "Xd1", "Z", "Yobs" ) %in% colnames(rr2) ) )
  
  rr2 = single.CRT.fixed.district( grid, n.runs = 2 )
  
  expect_true( is.data.frame( rr2 ) )
  expect_true( !is.null( rr2$method ) )
  
  df = run.scenario.fixed.district( grid, n.runs = 3, R = 2, parallel = FALSE )  
  
  expect_equal( n_distinct( df$run ), 2 )
  expect_equal( n_distinct( df$subrun ), 2*3 )
  
    
})



