
# Compare all estimators for different designs

library( clusterRCT )
library( tidyverse )

set.seed( 40440 )


# Make a series of interesting datasets ----
datasets = list()


# Matched pairs
blk = blkvar::generate_blocked_data( n_k = c( 5, 10, 20, 10, 5, 10 ) )
blks = tibble( B = paste0( "B", 1:6 ),
               D = rep( 1:3, each=2 ),
               Z = rep( c(0,1), 3 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk, control_formula = ~ X )
dsc

datasets$pair1 = blk


# singleton tx or co design

set.seed( 4343 )
blk = blkvar::generate_blocked_data( n_k = c( 5, 10, 20, 10, 11, 13,
                                              10, 10, 20, 10,
                                              30, 10, 5, 10, 30, 15 ) )
head( blk )
blks = tibble( B = paste0( "B", 1:16 ),
               D = rep( 1:3, c(6,4,6) ),
               Z = randomizr::block_ra(blocks=D, block_m = c(3,1,3) ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$one_sing = blk


# One block with 2 clusters
n_k = c( 5,  10, 10, 20, 30, 40,
         20, 10,
         10, 20, 30, 10, 40, 10,
         40, 10, 5, 10, 20, 30 )
blk = blkvar::generate_blocked_data( n_k = n_k )
head( blk )
blks = tibble( B = paste0( "B", 1:length(n_k) ),
               D = rep( c( 1, 2, 3, 4 ), c(6,2,6,6) ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$one_pr = blk


# One block with 3 clusters
n_k = c( 5,  10, 10, 20, 30, 40,
         20, 10,
         10, 20, 30, 10, 40, 10,
         40, 10, 5, 10, 20, 30, 10 )
blk = blkvar::generate_blocked_data( n_k = n_k )
head( blk )
blks = tibble( B = paste0( "B", 1:length(n_k) ),
               D = rep( c( 1, 2, 3, 4 ), c(6,3,6,6) ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$one_triad = blk



# One block with 2 tx and 2 co units in it.

set.seed( 50550 )

blk = blkvar::generate_blocked_data( n_k = c( 5,  10, 10, 20, 30, 40,
                                              20, 10,  5, 10, 10,
                                              10, 20, 30, 10, 40, 10 ) )
head( blk )
blks = tibble( B = paste0( "B", 1:17 ),
               D = rep( c( 1, 2, 3 ), c(6,5,6) ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$one_dbl = blk


# Now cut down to smaller experiment and add some covariates!
# table( blk$B, blk$D, blk$Z )
b2 <- blk %>%
  filter( B %in% paste0( "B", c( 1,4, 11,7, 12,13,
                                 2,3,5,10,8,14,15,17) ) )
describe_clusterRCT( Yobs ~ Z | B | D, data=b2 )

datasets$dbl2 = b2



# All blocks have  2 tx and 2 co units in it.
n_k = c( 5,  10, 30, 40,
         20, 10,  5, 10,
         10, 20, 40, 10 )
blk = blkvar::generate_blocked_data( n_k = n_k )
head( blk )
blks = tibble( B = paste0( "B", 1:length(n_k) ),
               D = rep( c( 1, 2, 3 ), c(4,4,4) ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$dbl2 = blk




# degrees of freedom at margin

n_k = c( 5,  10, 30, 40, 10, 10,
         20, 10,  5, 10, 10, 20,
         10, 20, 40, 10, 30, 20 )
blk = blkvar::generate_blocked_data( n_k = n_k )
head( blk )
blks = tibble( B = paste0( "B", 1:length(n_k) ),
               D = rep( c( 1, 2, 3 ), c(6,6,6) ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc
datasets$big = blk




## all triads -----

n_k = c( 5,  10, 30, 40, 10, 10,
         20, 10,  5, 10, 10, 20 )
blk = blkvar::generate_blocked_data( n_k = n_k )
head( blk )
blks = tibble( B = paste0( "B", 1:length(n_k) ),
               D = rep( c( 1, 2, 3, 4 ), 3 ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$triad = blk




# Matched pairs and matched doubles designs

blk = blkvar::generate_blocked_data( n_k = rep( 10, 18 ) )
head( blk )
blks = tibble( B = paste0( "B", 1:18 ),
               D = rep( 1:9, 2 ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$pair2 = blk


blk = blkvar::generate_blocked_data( n_k = rep( 10, 16 ) )
head( blk )
blks = tibble( B = paste0( "B", 1:16 ),
               D = rep( 1:4, 4 ),
               Z = randomizr::block_ra(blocks=D, prob=0.5 ) )
blk = left_join( blk, blks, by="B" ) %>%
  mutate( Yobs = ifelse( Z, Y1, Y0 ) )

blk$X = rnorm( nrow(blk) )

dsc = describe_clusterRCT( Yobs ~ Z | B | D, data=blk )
dsc

datasets$dbl2 = blk



# Assemble all datasets ----

length( datasets )
datasets <- datasets[ order( names(datasets) ) ]
names(datasets)

if ( FALSE ) {
  # Add in a group level covariate
  datasets <- map( datasets,
                   \(.x) {
                     .x %>%
                       mutate( X = Y0 * 0.5 + rnorm(n()) ) %>%
                       group_by( B, D ) %>%
                       mutate( Yobs = Yobs + rnorm(1, sd=1 ),
                               W = rnorm(1) ) %>%
                       ungroup() %>%
                       group_by( D ) %>%
                       mutate( Yobs = Yobs + rnorm(100) ) %>%
                       ungroup()

                   } )
} else {
  datasets <- map( datasets,
                   \(.x) {
                     .x %>%
                       group_by( B, D ) %>%
                       mutate( W = rnorm(1) ) %>%
                       ungroup()
                   } )
}

# Describe all datasets ----
descs <- map_df( datasets, describe_clusterRCT,
                 formula = Yobs ~ Z | B | D,
                 control_formula = ~ W, .id = "dataset" )
descs

# Print table
descs %>%
  dplyr::select( dataset, n, J, K, Jbar,
                 matched_pairs, num_singletons, num_doubletons,
                 cluster_ICC, block_ICC, R2.1, R2.2 ) %>%
  rename( pairs = matched_pairs,
          num_sing = num_singletons,
          num_dbl = num_doubletons ) %>%
  mutate( df = J - 2*K  ) %>%
  relocate( df, .after = Jbar ) %>%
  knitr::kable( digits = 2 )

# Functions to analyze all datasets ----

analyze_dat <- function( datasets, ... ) {
  res_X <- map_df( datasets,
                   \(.x) {
                     compare_methods( Yobs ~ Z | B | D,
                                      data=.x,
                                      ... )
                   },
                   .id = "dataset" ) %>%
    mutate( ATE = ifelse( is.na( ATE_hat ), "*", "." ),
            SE = ifelse( is.na( SE_hat ), "*", "." ) ) %>%
    group_by( method ) %>%
    mutate( n_ATE = sum( ATE == "*" ),
            n_SE = sum( SE == "*" ) ) %>%
    ungroup()

}

ATE_table <- function( res ) {
  res %>%
    dplyr::select( dataset, method, n_ATE, ATE ) %>%
    pivot_wider( names_from = dataset, values_from = ATE ) %>%
    arrange( n_ATE, method ) %>%
    knitr::kable( )

}

SE_table <- function( res ) {
  res %>%
    dplyr::select( dataset, method, SE, n_SE ) %>%
    pivot_wider( names_from = dataset, values_from = SE ) %>%
    arrange( n_SE, method ) %>%
    knitr::kable( )

}

# Analyze all datasets ----

res <- analyze_dat( datasets )

ATE_table( res )
SE_table( res )


res_X <- analyze_dat( datasets,
                      control_formula = ~ X )
ATE_table( res_X )
SE_table( res_X )

if ( FALSE ) {
  compare_methods( Yobs ~ Z | B | D,
                   datasets$pair1, control_formula = ~ X )
}


add_table <- function( res, res2 ) {
  res <- dplyr::select( res,
                        dataset, method, ATE, SE,
                        n_ATE, n_SE )
  res <- left_join( res, res2,
                    by=c( "dataset", "method" ) ) %>%
    mutate( ATE = paste0( ATE.x, ATE.y ),
            SE = paste0( SE.x, SE.y ),
            n_ATE = n_ATE.x + n_ATE.y,
            n_SE = n_SE.x + n_SE.y )

  res %>%
    dplyr::select( dataset, method, ATE, SE,
                   n_ATE, n_SE )
}


resJnt <- add_table( res, res_X )


ATE_table( resJnt )
SE_table( resJnt )



## Group level covariate ----

res_W <- analyze_dat( datasets,
                      control_formula = ~ W )
ATE_table( res_W )
SE_table( res_W )

res_XW <- analyze_dat( datasets,
                      control_formula = ~ X + W )
ATE_table( res_XW )
SE_table( res_XW )

resJnt2 <- add_table( resJnt, res_W ) %>%
  add_table( res_XW )


ATE_table( resJnt2 )
SE_table( resJnt2 )


if ( FALSE ) {
  compare_methods( Yobs ~ Z | B | D,
                   datasets$pair1, control_formula = ~ W )
  describe_clusterRCT( Yobs ~ Z | B | D,
                   datasets$pair1, control_formula = ~ X )

  describe_clusterRCT( Yobs ~ Z | B | D,
                       datasets$pair1, control_formula = NULL )

  describe_clusterRCT( Yobs ~ Z | B | D,
                   datasets$pair2, control_formula = NULL )

  compare_methods( Yobs ~ Z | B | D,
                   datasets$pair2, control_formula = NULL )

  compare_methods( Yobs ~ Z | B | D,
                   datasets$pair2, control_formula = ~ W + X  )

}

# Looking at when point estimates are same or different ----


source( here::here( "group_estimators.R" ) )

res$runID = res$dataset
group_estimators(res)
group_estimators(res, var="SE_hat")

#debug( group_estimators )
res_X$runID = res_X$dataset
group_estimators(res_X)

res_W$runID = res_W$dataset
group_estimators(res_W)


mrg <- bind_rows( nocov = res,
                  X = res_X,
                  W = res_W,
                  XW = res_XW, .id = "covs" ) %>%
  mutate( runID = paste0( dataset, "_", covs ) )
group_estimators( mrg )
group_estimators( mrg, var="SE_hat" )
