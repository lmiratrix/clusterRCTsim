

source( here::here( "R/util.R" ) )


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
#' @param rerandomize If `TRUE`, then the data will be re-randomized,
#'   keeping the same proportions.  FALSE means it will keep the
#'   original randomization.
#' @return A data frame with the generated data.
generate_data_from_fingerprint <- function( structure, parameters, 
                                            rerandomize = TRUE,
                                            tau = 0.20,
                                            names = c( "Yobs", "Z", "clusterID", "blockID") ) {
  
  if ( !all( c("blockID", "clusterID", "Z") %in% colnames( structure ) ) ) {
    stop( "structure must have columns blockID, clusterID, and Z" )
  }
  if ( "n" %in% colnames( structure ) ) {
    structure = make_skeleton(structure)
  }
  
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
  
  if ( rerandomize ) {
    structure$Z = rerandomize( structure$Z, structure$clusterID, structure$blockID )
  }
  
  # TODO: Possibly group mean center individual level covariate.  Good
  # idea?
  dat <- structure %>% 
    mutate( e = rnorm( n(), mean = 0, sd = sigma_e ),
            X = X ) %>%
    group_by( clusterID ) %>%
    mutate( X = X - mean(X) ) %>%
    ungroup() %>%
    mutate(
      Y0 = block_int[ blockID ] + cluster_int[ clusterID ] + e + W[clusterID] * beta_W + X * beta_X,
      tau = tau,
      Y1 = Y0 + tau,
      Yobs = Z * Y1 + (1-Z) * Y0 )
  
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
  
  dat <- dat %>%
    rename( !!names[1] := Yobs,
            !!names[2] := Z,
            !!names[3] := clusterID,
            !!names[4] := blockID ) %>%
    dplyr::select( -e )
  
  dat
}


#' Given data, return a data frame with the structure of the data
get_structure <- function( formula, data ) {
  data = clusterRCT:::make_canonical_data( formula=formula, data=data )
  data %>% 
    group_by( blockID, clusterID ) %>%
    summarize( Z = mean( Z == 1 ), 
               n = n(),
               .groups="drop" ) %>%
    ungroup() 
}


make_skeleton <- function( structure ) {
  K = length( unique( structure$blockID ) )
  J = length( unique( structure$clusterID ) )
  
  tibble( blockID = rep( structure$blockID, structure$n ),
          clusterID = rep( structure$clusterID, structure$n ),
          Z = rep( structure$Z, structure$n ) )
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
  head( fakeCRT )
  fakeCRT$X2 = rnorm( nrow(fakeCRT) )
  dd = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                            control_formula = ~ X.jk + C.ijk + X2 )
  dd
  
  struc <- get_structure( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
  struc
  
  make_skeleton(struc)
  nrow(fakeCRT)
  
  fakeStar <- generate_data_from_fingerprint( struc, dd )
  head( fakeStar )
  star_stat <- describe_clusterRCT( Y ~ Z | clusterID | blockID, 
                                    control_formula = ~ X1 + X2 + W,
                                    data=fakeStar )
  
  bind_rows( dd, star_stat )
}

