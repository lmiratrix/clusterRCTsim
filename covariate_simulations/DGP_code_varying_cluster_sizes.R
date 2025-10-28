

#' Randomly drop units so the cluster sizes vary.
#'
#' @param corr If nonzero, will drop units to induce correlation with
#'   the cluster-level impact or control outcome mean and the cluster
#'   size.
tilt_balanced_PUMP_data <- function( data, tilt=0.3, corr=0, byY0=TRUE ) {
  
  stats <- data %>% 
    group_by( S.id, D.id ) %>%
    summarise( Ybar0 = mean( Y0 ),
               tau = mean(Y1 - Y0 ),
               n = n(), .groups = "drop" )
  
  if ( byY0 ) {
    stats$tau = stats$Ybar0
  }
  
  stats$rem = round( runif( nrow(stats), tilt, 1 ) * stats$n )
  stats$rem = induce_cor( stats$tau, stats$rem, corr )
  
  # Now drop a random subset of 'rem' rows of each cluster in the data
  data2 <- data %>%
    left_join( stats %>%
                 select( S.id, D.id, rem ), 
               by = c("S.id", "D.id") )
  #summary( data2$rem )
  
  data2 <- data %>%
    left_join(select(stats, S.id, D.id, rem), by = c("S.id", "D.id")) %>%
    group_by(S.id, D.id) %>%
    slice(sample.int(
      n(),
      size = pmax( 1, coalesce(first(rem), n()) ),
      replace = FALSE
    )) %>%
    ungroup()
  
  
  data2
}


#' Make blocked cluster RCT data
#'
#' Given a list of parameters, make a dataset for a blocked cluster
#' RCT with a covariate that explains some of the variation in
#' outcomes. The covariate is generated to be independent of treatment
#' assignment.
#'
#' @param param_list A list with the following elements:
#'   - K: Number of blocks
#'   - J: Number of clusters
#'   - nbar: Average cluster size
#'   - MDES: Minimum detectable effect size (treatment effect)
#'   - ICC.3: Block-level ICC
#'   - ICC.2: Cluster-level ICC
#'   - R2.2: Proportion of cluster-level variance explained by covariate
#'   - R2.1: Proportion of individual-level variance explained by covariate
#' @param scale_W Scale parameter for the covariate so the X
#'   regression coefficient will differ from the W one (default 0.5)
#' @return A data frame with the generated data.
#' 
make_covariate_data <- function( param_list, scale_W = 0.5 ) {
  
  gamma.k = rnorm( param_list$K, 
                   0, sqrt( param_list$ICC.3 ) )
  alpha.jk = rnorm( param_list$J*param_list$K, 
                    0, sqrt( (1-param_list$R2.2)*param_list$ICC.2 ) )
  e.ijk = rnorm( param_list$J * param_list$K * param_list$nbar,
                 0, 
                 sqrt( (1-param_list$R2.1)*(1 - param_list$ICC.2 - param_list$ICC.3 ) ) )
  
  W = rnorm( param_list$J * param_list$K, 0, 1 ) * scale_W
  X = rnorm( param_list$J * param_list$K * param_list$nbar, 0, 1 )
  
  beta_W = sqrt( param_list$R2.2*param_list$ICC.2 ) / scale_W
  
  beta_X = sqrt( param_list$R2.1*(1 - param_list$ICC.2 - param_list$ICC.3 ) )
  
  clust = tibble( S.id = 1:(param_list$J*param_list$K),
                  D.id = rep( 1:param_list$K, each=param_list$J ),
                  gamma.k = gamma.k[D.id],
                  alpha.jk = alpha.jk,
                  W = W ) %>%
    expand_grid( i = 1:param_list$nbar )
  
  clust <- clust %>%
    mutate( X = X ) %>%
    group_by( D.id, S.id  ) %>%
    mutate( Xbar = mean(X),
            Xc = X - mean(X) ) %>%
    ungroup() %>%
    mutate(
      Y0 = gamma.k + alpha.jk + beta_W*W + beta_X*Xbar + beta_X*Xc + e.ijk,
      Y1 = Y0 + param_list$MDES,
      T.x = randomizr::block_and_cluster_ra( blocks=D.id,
                                             clusters=S.id,
                                             prob=0.5),
      Yobs = ifelse( T.x, Y1, Y0 ) )
  
  
  clust
}


