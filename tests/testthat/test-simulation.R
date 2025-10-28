


# Test the simulation support functions
library( testthat )
library( tidyverse )

source( here::here( "R/simulation_support_code_crt.R" ) )
source( here::here( "R/data_generation_code.R" ) )


##### Test the DGP code #####

test_that( "single.CRT() works", {
  
  set.seed( 44440 )
  rr = single.CRT( K = 10, 
                   J.bar = 10,
                   p.bar = 0.5,
                   ATE.block = 0.2, 
                   sd.block = 0.1,
                   ICC.block = 0.2,
                   variable.p = T,
                   
                   # cluster-level params
                   n.bar = 50, 
                   sd.cluster = 0.1, 
                   ICC.cluster = 0.2,
                   cv.cluster = 0.2,
                   cv.block = 0.2,

                    
                   # size and correlation params
                   size.impact.cor.block = 0,
                   size.intercept.cor.block = 0.5,
                   size.impact.cor.cluster = 0,
                   proptx.impact.cor.block = 0,
                   
                   include_method_characteristics=FALSE)
  
  
  head( rr )
  
  expect_equal( length( unique( rr$subrun ) ), 3 )
  expect_equal( sd( rr$ATE.finite.dist.cluster ), 0 )
  expect_equal( sd( rr$ATE.finite.cluster ), 0 )
  expect_equal( sd( rr$ATE.finite.student ), 0 )
  expect_equal( sd( rr$ATE.finite.dist.student ), 0 )
  
  g <- rr %>% group_by( method ) %>%
    summarize( sdATE = sd(ATE.hat),
               n = n() )
  expect_true( all( g$n == 3 ) )
  expect_true( all( g$sdATE > 0.001 ) )
  
  
  rr2 = single.CRT( K = 10, 
                    J.bar = 10,
                    p.bar = 0.5,
                    ATE.block = 0.2, 
                    sd.block = 0.1,
                    ICC.block = 0.2,
                    variable.p = T,
                    
                    # cluster-level params
                    n.bar = 50, 
                    sd.cluster = 0.1, 
                    ICC.cluster = 0.2,
                    cv.cluster = 0.2,
                    
                    # size and correlation params
                    size.impact.cor.block = 0,
                    size.intercept.cor.block = 0.5,
                    size.impact.cor.cluster = 0,
                    proptx.impact.cor.block = 0,
                    
                    include_MLM=FALSE)
  
  expect_true( !( "MLM" %in% rr2$method ) )
  expect_true( nrow( rr2 ) < nrow(rr) )
  
})



