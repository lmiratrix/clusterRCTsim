

#' Generate data from a "fingerprint" of a structure (list of cluster
#' IDs with associated sample sizes) and parameters (values of the
#' ICC, etc., as you would get from a `describe_data()` call).
#'
#' @param structure A data frame with columns `blockID`, `clusterID`,
#'   and `Z` (treatment assignment).
#' @param parameters A list with the following elements: `block_ICC`,
#'   `cluster_ICC`, `R2.1`, `R2.2`, `ncov.1`, `ncov.2`.
#' @param tau The treatment effect (considered constant across dataset
#'   at this time).
#'
#' @return A data frame with the generated data.
generate_data_from_fingerprint <- function( structure, parameters, 
                                            tau = 0.20 ) {
  
  K = length( unique( structure$blockID ) )
  J = length( unique( structure$clusterID ) )  
  n = nrow( structure )
  structure <- mutate( structure,
                       blockID = as.character(blockID),
                       clusterID = as.character(clusterID ) )
  
  block_int = rnorm( K, mean = 0, sd = sqrt( parameters$block_ICC ) )
  names( block_int ) = unique( structure$blockID )
  
  cluster_int = rnorm( J, mean = 0, sd = sqrt( parameters$cluster_ICC ) )
  names( cluster_int ) = unique( structure$clusterID )
  
  W = rnorm( J )
  names(W) = names( cluster_int )
  X = rnorm( n )
  
  beta_W = sqrt( parameters$R2.2 * parameters$cluster_ICC )
  beta_X = sqrt( (1-parameters$block_ICC-parameters$cluster_ICC) * parameters$R2.1 )
  
  sigma_e = sqrt( ( 1 - parameters$block_ICC - parameters$cluster_ICC ) * (1-parameters$R2.1 ) )
  
  
  # TODO: Possibly group mean center individual level covariate.  Good
  # idea?
  dat <- structure %>% 
    mutate( e = rnorm( n(), mean = 0, sd = sigma_e ),
            X = X ) %>%
    group_by( clusterID ) %>%
    mutate( X = X - mean(X) ) %>%
    ungroup() %>%
    mutate(
      Y = block_int[ blockID ] + cluster_int[ clusterID ] + Z * tau + e + W[clusterID] * beta_W + X * beta_X )
  
  # Make it so we end up with right number of covariates
  if ( parameters$ncov.2 > 0 ) {
    dat$X = dat$X + W[ dat$clusterID ]
  }
  
  X = blkvar:::make_covariate_set( dat$X, k = parameters$ncov.1, 
                                   var_prefix = "X" )
  dat$X = NULL
  dat <- bind_cols( dat, X )
  
  if ( parameters$ncov.2 > parameters$ncov.1 ) {
    W = blkvar:::make_covariate_set( rnorm(J), 
                                     k = parameters$ncov.2 - parameters$ncov.1, 
                                     var_prefix = "W" )
    W$clusterID = names( cluster_int )
    
    dat <- left_join( dat, W, by="clusterID" )
  }
  
  dat
  
}


# Some testing code for the above
if (FALSE) {
  structure = data.frame( blockID = c("A", "A", "B", "B", "B"),
                          clusterID = c("1", "2", "1", "2", "3"),
                          Z = c( 0, 1, 0, 1, 0) )
  parameters = list( block_ICC = 0.1, cluster_ICC = 0.2, R2.1 = 0.1, R2.2 = 0.2,
                     ncov.1 = 1, ncov.2 = 1 )
  generate_data_from_fingerprint( structure, parameters )
  
  
  # Generate data "like" fakeCRT
  data( fakeCRT )
  dd = describe_clusterRCT( Yobs ~ Z | cid | did, data = fakeCRT,
                            control_formula = ~ d )
  generate_data_from_fingerprint( dd$structure, dd$parameters )
  
}

