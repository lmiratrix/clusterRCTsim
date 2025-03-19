

source( here::here( "R/data_generation_code.R" ) )

#' Generate data from a "fingerprint" of a structure (list of cluster
#' IDs with associated sample sizes) and parameters (values of the
#' ICC, etc., as you would get from a `describe_data()` call).
#'
#' There is no correlation of varying cluster or block sizes with
#' outcome or treatment impact or proportion treated.  Everything is
#' independent and, thus, well behaved in principle.
#'
#' @param structure A data frame with columns `blockID`, `clusterID`,
#'   and `Z` (treatment assignment).  If `n` is also present, then the
#'   data will be expanded to have the correct number of rows.  If
#'   structure is NULL, then it will make data with roughly the same
#'   number of clusters, etc., as specified by values in the
#'   parameters list.
#' @param parameters A list with the following elements: `block_ICC`,
#'   `cluster_ICC`, `R2.1`, `R2.2`, `ncov.1`, `ncov.2`.  Additional
#'   elements will be ignored.
#' @param tau The treatment effect (considered constant across dataset
#'   at this time).
#' @param rerandomize If `TRUE`, then the data will be re-randomized,
#'   keeping the same proportions.  FALSE means it will keep the
#'   original randomization.
#' @param names List of what to call outcome, treatment, cluster, and
#'   block ID, in that order.
#'
#' @return A data frame with the generated data.
generate_data_from_fingerprint <- function( structure = NULL, parameters = NULL, 
                                            rerandomize = TRUE,
                                            tau = 0.0,
                                            names = c( "Yobs", "Z", "clusterID", "blockID") ) {
  
  stopifnot( !is.null( parameters ) )
  
  # Make the dataset of block and cluster IDs
  if ( is.null( structure ) ) {
    stopifnot( !is.null( parameters ) )
    stopifnot( all( c("K", "n", "n.25", "J", "J.25", "nbar", "Jcv") %in% names(parameters) ) )
    K = parameters$K
    min.size = 2
    
    num_clusters <- cluster_distn(K, 
                                  J.bar = parameters$Jbar, 
                                  cv = parameters$Jcv,
                                  min.size = min.size )
    

    # This should allow variation in treatment
    probs <- generate_treatment_probabilities( K, 
                                               parameters$tx.avg,
                                               parameters$tx.cv * parameters$tx.avg )
    num_Z = round( num_clusters * probs )
    num_Z = ifelse( num_Z == 0, 1, ifelse( num_Z == num_clusters, num_clusters-1, num_Z  ) )
    
    Z = rep( rep( c(0,1), K ),
             as.numeric( t( cbind( num_Z, num_clusters - num_Z ) ) ) )
    stopifnot( length(Z) == sum(num_clusters) )
    structure = tibble( blockID = rep( paste0( "B", 1:K ), num_clusters ),
                        clusterID = rep( paste0( "C", 1:sum(num_clusters) ) ),
                        Z = Z )
    
    structure$n <- cluster_distn( nrow(structure), 
                                  J.bar = parameters$nbar, 
                                  cv = parameters$ncv, 
                                  min.size = 1 )
    
    # Check we have valid random assignment
    tt = table( structure$blockID, structure$Z )
    stopifnot( all( tt > 0 ) )
  } 
  
  if ( !all( c("blockID", "clusterID", "Z") %in% colnames( structure ) ) ) {
    stop( "structure must have columns blockID, clusterID, and Z" )
  }
  if ( "n" %in% colnames( structure ) ) {
    structure = make_skeleton(structure)
  }
  structure$clusterID = paste0( structure$blockID, ".", structure$clusterID)
  
  K = length( unique( structure$blockID ) )
  J = length( unique( structure$clusterID ) )  
  n = nrow( structure )
  structure <- mutate( structure,
                       blockID = as.character(blockID),
                       clusterID = as.character(clusterID ) )
  
  
  # Generate block and cluster intercepts
  block_int = rnorm( K, mean = 0, sd = sqrt( parameters$block_ICC ) )
  names( block_int ) = unique( structure$blockID )
  
  cluster_int = rnorm( J, mean = 0, sd = sqrt( parameters$cluster_ICC ) )
  names( cluster_int ) = unique( structure$clusterID )
  
  # Generate our covariates
  W = rnorm( J )
  names(W) = names( cluster_int )
  X = rnorm( n )
  
  beta_W = sqrt( parameters$R2.2 * parameters$cluster_ICC )
  
  beta_X = sqrt( (1-parameters$block_ICC-parameters$cluster_ICC) * parameters$R2.1 )
  sigma_e = sqrt( ( 1 - parameters$block_ICC - parameters$cluster_ICC ) * (1-parameters$R2.1 ) )
  
  if ( rerandomize ) {
    structure$Z = rerandomize( structure$Z, structure$clusterID, structure$blockID )
  }
  
  # TODO: I am group mean centering the individual level covariate.
  # Good idea?
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
  
  # Make it so we end up with right number of covariates by combining
  # our covariates.  The compare code, etc., will then generate group
  # mean centered covariates to split them back into two.
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
  
  
  # Make the final dataset from the above results
  dat <- dat %>%
    rename( !!names[1] := Yobs,
            !!names[2] := Z,
            !!names[3] := clusterID,
            !!names[4] := blockID ) %>%
    dplyr::select( -e )
  
  dat
}


#' Given a structure, make a full skeleton
#'
#' A "skeleton" has one row per level 1 observation.  A structure has
#' cluster-level statistics, including number of units and number of
#' treated units.
#'
#' @return dataframe, one row per level 1 unit.
make_skeleton <- function( structure ) {
  
  stopifnot( is.data.frame( structure ) )
  stopifnot( all( c("blockID", "clusterID", "Z", "n") %in% colnames( structure ) ) )
  
  #K = length( unique( structure$blockID ) )
  #J = length( unique( structure$clusterID ) )
  
  tibble( blockID = rep( structure$blockID, structure$n ),
          clusterID = rep( structure$clusterID, structure$n ),
          Z = rep( structure$Z, structure$n ) )
}


