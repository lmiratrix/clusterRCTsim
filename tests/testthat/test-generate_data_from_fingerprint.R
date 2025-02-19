
source( here::here( "R/generate_data_from_fingerprint.R" ) )

test_that("running sim from fingerprint", {
  
  library( clusterRCT )
  
  # Make 5 cluster dataset
  structure = data.frame( blockID = c("A", "A", "B", "B", "B"),
                          clusterID = c("1", "2", "1", "2", "3"),
                          Z = c( 0, 1, 0, 1, 0),
                          n = c( 3, 2, 2, 1, 2 ) )
  structure
  
  parameters = list( block_ICC = 0.1, cluster_ICC = 0.2, R2.1 = 0.1, R2.2 = 0.2,
                     ncov.1 = 1, ncov.2 = 1 )
  dd <- generate_data_from_fingerprint( structure, parameters )
  dd
  
  expect_equal( nrow(dd), sum( structure$n ) )
  
  # Generate data "like" fakeCRT
  data( fakeCRT )
  head( fakeCRT )
  fakeCRT$X2 = rnorm( nrow(fakeCRT) )
  dd = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                            control_formula = ~ X.jk + C.ijk + X2 )
  dd
  
  struc <- get_structure( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
  struc
  
  # make 1 row for each individual observation, based on structure
  make_skeleton(struc)
  nrow(fakeCRT)
  
  fakeStar <- generate_data_from_fingerprint( struc, dd )
  head( fakeStar )
  star_stat <- describe_clusterRCT( Yobs ~ Z | clusterID | blockID, 
                                    control_formula = ~ X1 + X2 + W,
                                    data=fakeStar )
  
  bb <- bind_rows( fakeCRT = dd, fakeStar = star_stat, .id="dataset" )
  bb
  expect_equal( bb[1,2:8], bb[2,2:8] )

})
