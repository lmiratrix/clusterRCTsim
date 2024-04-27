

# This script has all the support code to run large simulations to
# look at stability of the different SE estimators and true
# performances of the ATE estimators.
#
# This script has functions to run a single trial, generate a list of
# all scenarios, and a few other things.

library( tidyverse )
library( purrr )
library( multiwayvcov )

library( clusterRCT )


#### Utility functions #####


# Print things to console using sprintf.
scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}



#### Data-generating process ####

# generate variable cluster sizes
#  - analogue to blkvar::block_distn()
cluster_distn <- function(K, J.bar, size.ratio) {
  stopifnot(J.bar > 5)
  nk <- blkvar::block_distn(J=K, n.bar=J.bar-4, size.ratio)
  nk[nk<4] <- 4
  
  nk
}

# preserving order of x, rearrange y to correlate x
#  via a magic hacky fuzzy sort
#  - `cor` sets strength of correlation
#  - returns: rearranged y
induce_cor <- function(x, y, cor) {
  
  # New version: Copula!
  rnk = MASS::mvrnorm( length(y), mu=c(0,0), Sigma=matrix(c(1,cor,cor,1),2,2) )
  
  rnk <- rnk[ order(rnk[,1]), ]
  rnk <- rnk[ rank( x ), ]
  
  yn <- sort(y)[ rank( rnk[,2] ) ]
  
  #cor( x, yn )
  yn
  # Old version
  #temp <- y[order(cor^2 * y + (1-cor^2) * rnorm(length(y), sd=sd(y)))]
  #y <- temp[ rank( sign(cor)*x ) ]
}

if ( FALSE ) {
  
  cor = 1
  x = 10*rnorm( 10000 )
  y = 10*rnorm( 10000 )
  yn = induce_cor( x, y, cor = 1 )
  
  cor( x, yn )
  
  expect_equal( cor( x, yn ), 1, tolerance = 0.05 )
  
  yn = induce_cor( x, y, cor = -0.5 )
  cor( x, yn )
  expect_equal( cor( x, yn ), -0.5, tolerance = 0.05 )
  
}



#' Generate set of covariates that have given R2 to outcome y.
#' 
#' @param k Number of covariates
generate_covariate_block <- function( y, var_y, R2, k = 1,
                                      var_prefix = "X" ) {
  
  stopifnot( k >= 1 )
  stopifnot( R2 >= 0 && R2 <= 1 )
  #stopifnot( var_y > 0 )
  
  n = length( y )
  
  var_e = 1
  var_x = 1
  if ( var_y > 0 ) {
    var_e = (1-R2)*var_y / (R2 + (1-R2)*var_y)
    var_x = var_y + var_e - var_y * var_e
  }
  
  X1 = ( sqrt(1-var_e) * y + rnorm( n, mean=0, sd=sqrt( var_e ) ) ) / sqrt(var_x)
  if ( k > 2 ) {
    X = as.data.frame( cbind( X1, matrix( rnorm( n * (k-1), mean=0, sd=1 ), n, k-1 ) ) ) %>%
      as_tibble()
  } else {
    X = tibble( X1 = X1 )
  }
  colnames( X ) = paste0( var_prefix, 1:k )
  return ( X )
}


#### Code to build districts ####

#' Build set of districts with cluster sizes and proportion treated,
#' etc.
#'
#' @param variable.p Variable proportion treated per district, uniform
#'   around the p.bar value, reaching 2/3 the distance from 0 or 1
#'   (whichever is closer).
#' @param size.ratio.district IF variable.J is TRUE, use this as a
#'   measure of how much variation there is in number of clusters per
#'   district.
#'
#' @return dataframe, one row per district.
build_districts <- function( K, J.bar, p.bar=0.5,
                             ATE, sd.ATE,
                             ICC,
                             variable.J = FALSE,
                             variable.p = FALSE,
                             size.ratio = 0.6, 
                             size.impact.cor = 0,
                             proptx.impact.cor = 0,
                             size.intercept.cor = 0,
                             R2 = NULL,
                             num.covariates = 0 + !is.null(R2) ) {
  
  stopifnot("size.impact.cor must be in (-1,1)" = 
              size.impact.cor >= -1 & size.impact.cor <= 1)
  stopifnot("proptx.impact.cor must be in (-1,1)" = 
              proptx.impact.cor >= -1 & proptx.impact.cor <= 1)
  stopifnot("size.intercept.cor must be in (-1,1)" = 
              size.intercept.cor >= -1 & size.intercept.cor <= 1)
  
  # district J values
  if (variable.J) {
    num_clusters <- cluster_distn(K, J.bar, size.ratio=size.ratio)
  } else {
    num_clusters <- rep(J.bar, K)
  }
  
  # district proportions treated
  if (variable.p) {
    pp <- (2/3) * min(p.bar, 1-p.bar)
    prop_tx_clusters <- runif(K, min=p.bar-pp, max=p.bar+pp)
  } else {
    prop_tx_clusters <- rep(p.bar, K)
  }
  
  # district-level parameters
  #  - note: ensure var of district-level effects is exactly sd.ATE^2, ICC
  grid <- tibble(
    did = 1:K,
    # district.tau = rnorm(K, mean=ATE, sd=sd.ATE))
    district.tau = MASS::mvrnorm(K, mu=ATE, Sigma=sd.ATE^2, 
                                 empirical=TRUE) %>% 
      as.vector(),
    district.alpha = MASS::mvrnorm(K, mu=0, Sigma=ICC,
                                   empirical=TRUE) %>% 
      as.vector())
  
  # generate correlations between district.tau and num_clusters & prop_tx_clusters
  #  - uses some magic hacks; idea is to keep order of district.tau,
  #    and to reorder num_clusters & prop_tx_clusters, to avoid automatically
  #    correlating the two of them
  if (size.impact.cor != 0) {
    stopifnot("Cannot correlate district sizes and average impacts if sizes do not vary" = 
                variable.J)
    num_clusters <- induce_cor( x = grid$district.tau,
                                y = num_clusters,
                                cor = size.impact.cor)
  }
  if (proptx.impact.cor != 0) {
    stopifnot("Cannot correlate district proportions treated and average impacts if proportions treated do not vary" = 
                variable.p)
    prop_tx_clusters <- induce_cor(x = grid$district.tau,
                                   y = prop_tx_clusters,
                                   cor = proptx.impact.cor)
    # p_temp <- prop_tx_clusters[order(
    #   proptx.impact.cor^2 * prop_tx_clusters + 
    #     (1-proptx.impact.cor^2) * rnorm(K, sd=sd(prop_tx_clusters)) )]
    # prop_tx_clusters <- p_temp[ order(order(sign(proptx.impact.cor)*grid$district.tau)) ]
  }
  if (size.intercept.cor != 0) {
    stopifnot("Cannot correlate district sizes and intercepts if ICC = 0" = 
                (ICC != 0))
    grid$district.alpha <- induce_cor(x = num_clusters,
                                      y = grid$district.alpha,
                                      cor = size.intercept.cor)
  }
  
  grid$J <- num_clusters
  grid$p <- prop_tx_clusters
  
  # force rerandomize.data() to assign at least two treated/control clusters
  #  to each district, by applying floor and ceiling to p, depending on J
  #   - i.e., force 2 <= J*p <= J-2
  grid <- grid %>% 
    mutate(p = ifelse(J*p < 2, 2/J, p),
           p = ifelse(J*p > J-2, (J-2)/J, p))
  
  if ( num.covariates > 0 ) {
    if ( is.null(R2 ) ) {
      R2 = 0
    }
    cov_block <- generate_covariate_block( grid$district.alpha,
                                           ICC, R2, k = num.covariates, 
                                           var_prefix = "Xd" )
    grid = bind_cols( grid, cov_block )
  }
  
  return( grid )
}



if ( FALSE ) {
  grid <- build_districts( K = 10, J.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.ATE = 0.1,
                           ICC = 0.2,
                           variable.J = TRUE, variable.p = TRUE )
  
  grid
}





#### Code to simulate full data set ####


# Generate a cluster-randomized trial, with multiple districts
generate_districts <- function(
    # district-level params
  K, J.bar, p.bar=0.5,
  ATE.district, sd.district,
  ICC.district,
  variable.J = FALSE,
  variable.p = FALSE,
  R2.district = NULL,
  num.cov.district = 0 + !is.null(R2.district),
  
  # correlation params at district level
  size.ratio.district = 0.6, 
  size.impact.cor.district = 0,
  proptx.impact.cor.district = 0,
  size.intercept.cor.district = 0,
  
  # school-level params
  n.bar, sd.school, ICC.school = 0.20,
  variable.n = FALSE,
  R2.school = NULL,
  num.cov.school = 0 + !is.null(R2.school),
  
  # correlation params at school level
  size.ratio.school = 0.6,
  size.impact.cor.school = 0,
  cor.strength.school = 0.75,
  
  # Individual-level params
  R2.indiv = NULL,
  num.cov.indiv = 0 + !is.null(R2.indiv) ) {
  
 
  stopifnot( ICC.district + ICC.school < 1 )
  
  stopifnot(size.impact.cor.school %in% c(-1,0,1))
  
  
  grid <- build_districts( K = K, J.bar = J.bar, p.bar = p.bar,
                           ATE = ATE.district, sd = sd.district,
                           ICC = ICC.district,
                           variable.J = variable.J, variable.p = variable.p,
                           size.ratio = size.ratio.district,
                           size.impact.cor = size.impact.cor.district,
                           proptx.impact.cor = proptx.impact.cor.district,
                           size.intercept.cor = size.intercept.cor.district,
                           R2 = R2.district,
                           num.covariates = num.cov.district )
  
  
  # make student-level dataset
  sim.data <- grid %>%
    rowwise() %>%
    mutate(sim_data = list(
      blkvar::generate_multilevel_data( 
        n.bar = n.bar, 
        J = J,
        p = p,
        gamma.00 = district.alpha,
        gamma.10 = district.tau,
        tau.11 = sd.school^2,
        ICC = ICC.school / (1-ICC.district),
        rho2.0W = R2.school,
        rho2.1W = 0,
        R2.X = R2.indiv,
        num.W = num.cov.school,
        num.X = num.cov.indiv,
        variable.n = variable.n,
        finite.model = FALSE,
        size.impact.correlate = size.impact.cor.school,
        correlate.strength = cor.strength.school,
        size.ratio = size.ratio.school,
        proptx.impact.correlate = FALSE,
        zero.corr = TRUE,
        varY0 = 1 - ICC.district
      ))) %>%
    tidyr::unnest(cols=sim_data) %>% 
    dplyr::select(-J, -p, -district.tau, -district.alpha)
  
  # sim.data %>% filter(Z==0) %>% summarize(sd = sd(Yobs))
  
  # ensure that sids are unique
  sid_key <- sim.data %>%
    dplyr::select(did, sid) %>%
    unique() %>%
    mutate(sid_new = 1:n())
  
  sim.data <- sim.data %>%
    left_join(sid_key, by=c("did", "sid")) %>%
    mutate(sid = sid_new) %>%
    dplyr::select(-sid_new)
  
  # record J, p, district.tau for each district
  attr(sim.data, "district_key") <- grid
  
  # Not needed since randomization happens within generate_multilevel_data_no_cov:
  # assign cluster-randomized treatment
  # sim.data <- rerandomize.data(sim.data)
  
  sim.data
}



if (FALSE) {
  df <- generate_districts(
    K = 25,
    J.bar = 50, 
    p.bar = 0.5,
    ATE.district = 0.2, 
    sd.district = 0.1,
    ICC.district = 0.5,
    variable.J = T,
    variable.p = T,
    
    # correlation params at cluster level
    size.ratio.district = 0.6, 
    size.impact.cor.district = 0.2,
    proptx.impact.cor.district = 0.2,
    size.intercept.cor.district = 0.2,
    
    # school-level params
    n.bar = 50, 
    sd.school = 0.1, 
    ICC.school = 0.2,
    variable.n = T,
    
    # correlation params at individual level
    size.ratio.school = 0.6,
    size.impact.cor.school = 0,
    cor.strength.school = 0.75
  )
  
  describe_clusterRCT(Yobs ~ Z | sid | did, df)
  df %>% 
    group_by(did) %>% 
    summarize(mn = mean(Y0), n=n()) %>% 
    summarize(cor = cor(mn, n))
}

