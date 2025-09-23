

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

source( here::here( "R/util.R" ) )
source( here::here( "R/generate_nonnegative_integers.R" ) )

#### Utility functions #####


# Print things to console using sprintf.
scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}



#### Data-generating process ####

# generate variable cluster sizes
cluster_distn <- function(K, J.bar, cv, min.size = 4) {
  stopifnot(J.bar >= min.size)
  sd = cv * J.bar
  nk <- min.size + generate_nonnegative_integers( n = K, 
                                                  mean = J.bar-min.size, 
                                                  sd = sd )
  #nk[nk<min.size] <- min.size
  stopifnot( all( nk >= min.size ) )
  nk
}



# preserving order of x, rearrange y to correlate x
#  via a magic hacky fuzzy sort
#  - `cor` sets strength of correlation
#  - returns: rearranged y
induce_cor <- function(x, y, cor, empirical = FALSE ) {
  
  # New version: Copula!
  rnk = MASS::mvrnorm( length(y), 
                       mu=c(0,0), 
                       Sigma=matrix(c(1,cor,cor,1),2,2), 
                       empirical = empirical)
  
  rnk <- rnk[ order(rnk[,1]), ]
  rnk <- rnk[ rank( x, ties.method="first" ), ]
  
  yn <- sort(y)[ rank( rnk[,2] ) ]
  
  #cor( x, yn )
  yn
  # Old version
  #temp <- y[order(cor^2 * y + (1-cor^2) * rnorm(length(y), sd=sd(y)))]
  #y <- temp[ rank( sign(cor)*x ) ]
}


####      testing code ----
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
  
  
  x = sample( 10 )
  y = 10 + sample( 10 )
  x = qnorm( seq( 0.01, 0.99, length.out = length(x) ) )
  y = qnorm( seq( 0.01, 0.99, length.out = length(y) ) )
  cor( x, y )
  
  yn = induce_cor( x, y, cor = 0.80 )
  yn
  cor( x, yn )
  
  yn = induce_cor( x, y, cor = 0.80, empirical = TRUE )
  cor( x, yn )
  
  p1 = replicate( 1000, {
    yn = induce_cor( x, y, cor = 0.80, empirical = FALSE )
    cor( x, yn )
  } )
  
  p2 = replicate( 1000, {
    yn = induce_cor( x, y, cor = 0.80, empirical = TRUE )
    cor( x, yn )
  } )
  
  mean( p1 )
  sd( p1 )
  mean( p2 )
  sd( p2 )
  
  
  cor = 0.72
  rnk = MASS::mvrnorm( 100, mu=c(0,0), Sigma=matrix(c(1,cor,cor,1),2,2), empirical = TRUE)
  cor( rnk[,1], rnk[,2] )  
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


#' Describe a set of district parameters by printing to screen
describe_district_grid <- function( grid ) {
  
  K = nrow( grid )
  J.bar = mean( grid$J )
  sd.J = sd( grid$J )
  p.bar = mean( grid$p )
  sd.p = sd( grid$p )
  n.bar = weighted.mean( grid$n.bar, grid$J )
  sd.nbar = sd( grid$n.bar )
  
  if ( K == 1 ) {
    sd.J = 0
    sd.p = 0
    sd.nbar = 0
  }
  ATE.district = mean( grid$district.tau )
  sd.district = sd( grid$district.tau )
  ICC.district = sd( grid$district.alpha )
  
  size.impact.cor.district = ifelse( sd.district != 0 && sd.J != 0, 
                                     cor( grid$J, grid$district.tau ), NA )
  size.intercept.cor = ifelse( sd.J != 0 && ICC.district != 0,
                               cor( grid$J, grid$district.alpha ), NA )
  proptx.impact.cor.district = ifelse( sd.district != 0 && sd( grid$p ) != 0, 
                                       cor( grid$p, grid$district.tau ), NA )
  if ( sd.J > 0 && sd.nbar > 0 ) {
    nbar.cor.J = cor( grid$n.bar, grid$J )
  } else {
    nbar.cor.J = NA
  }
  
  taubarJ = weighted.mean( grid$district.tau, w=grid$J )
  taubarn = weighted.mean( grid$district.tau, w=grid$J * grid$n.bar )
  
  
  #         district ATE(J)=%.2f\tdistrict ATE(n)=%.2f
  
  scat( "K=%d\tJ.bar=%.2f (%.2f, cv=%.2f)\tn.bar=%.2f (%.2f cv=%.2f)\tp.bar=%.2f (%.2f)\n",
        K, 
        J.bar, sd.J, sd.J / J.bar, 
        n.bar, sd.nbar, sd.nbar / n.bar, 
        p.bar, sd.p )
  scat( "District stats:\n\tdistrict ATE=%.2f (%.2f)\tdistrict ICC=%.2f\t
\tcor(J,tau)=%.2f\t\tcor(ptx,tau)=%.2f
\tcor(J,alpha)=%.2f\t\tcor(nbar, J)=%.2f\n",
        ATE.district, sd.district, 
        ICC.district,
        #  taubarJ, taubarn,
        size.impact.cor.district, proptx.impact.cor.district,
        size.intercept.cor,
        nbar.cor.J )
  
  
  # get some empirical estimates of some quantities
  dst = nrow(grid)
  dmx = sum(grid$J * grid$n.bar)
  rps = pmax( 1, round( min( 400/dst, 100000 / dmx ) ) )
  
  
  dd = map_dfr( 1:rps, ~generate_data_for_districts( grid ), .id="run" )
  dd <- mutate( dd,
                did = paste0( run, "_", did ),
                sid = paste0( run, "_", sid ) )
  #dsc = describe_clusterRCT(Yobs ~ Z | sid | did, dd)
  #dsc
  
  
  # Print out school level stuff
  dagg <- dd %>% group_by( did, sid ) %>%
    summarise( tau = mean( Y1 - Y0 ),
               n = n(), .groups = "drop" )
  if ( sd( dagg$tau ) > 0 && sd( dagg$n ) > 0 ) {
    corntau = cor( dagg$tau, dagg$n )
  } else {
    corntau = NA
  }
  
  J = length( unique( dd$sid) )
  
  
  cat( "Cluster stats:\n")
  scat( "\tschool ICC=%.2f (%.2f)\n\tavg cluster size = %.1f (range %.0f--%.0f, sd = %.1f, cv=%.2f)\n\tcor(n,tau) = %.2f (across blocks)\n",
        mean(grid$ICC.school), sd(grid$ICC.school),
        mean( dagg$n ), 
        min(dagg$n), max(dagg$n), sd(dagg$n), sd(dagg$n) / mean(dagg$n),
        #quantile(dagg$n, 0.25), quantile(dagg$n, 0.75), sd(dagg$n),
        corntau )
  
  #M = lmer( tau ~ 1 + (1|did), data=dagg )
  #arm::display(M)
  
  daggagg <- dagg %>%
    group_by( did ) %>%
    summarise( sd_tau = sd( tau ),
               tau = mean( tau ),
               n = mean( n ) )
  
  
  cat( "ATE information:\n" )
  scat( "\tWithin-district cluster ATE and variation: %.2f (%.2f)\n",
        mean(daggagg$tau),
        sqrt(mean(daggagg$sd_tau^2)) )
  scat( "\tTotal cluster ATE: %.2f (%.2f)\n",
        mean( dagg$tau ), 
        sd( dagg$tau ) )
  scat( "\tTotal individual ATE: %.2f (%.2f)\n",
        mean( dd$Y1 - dd$Y0 ), 
        sd( dd$Y1 - dd$Y0 ) )
  
  #cat( glue::glue( "\t\n" ) )
  #attr(dd, "district_key") <- grid
  #d2 = describe.data( dd, print = FALSE )
  #print( knitr::kable( d2, digits = 2 ) )
  cat( "\n\t" )
  cat( glue::glue( " * Some stats estimated on {nrow(dd)} generated observations across {J} clusters.") )
  cat( "\n" )
  
  
  stats = tibble( K = K, J = sum(grid$J), 
                  mean.J = J.bar, sd.J = sd.J, cv.J = sd.J / J.bar,
                  mean.p = p.bar, sd.p = sd.p,
                  mean.nbar = n.bar, sd.nbar = sd.nbar,
                  ATE.district = ATE.district, 
                  sdATE.district = sd.district,
                  ICC.district = ICC.district,
                  mean.n = mean( dagg$n ),
                  sd.n = sd( dagg$n ),
                  cv.n = sd( dagg$n ) / mean( dagg$n ),
                  mean.total.cluster.ATE = mean( dagg$tau ),
                  sd.total.cluster.ATE = sd( dagg$tau ),
                  mean.cluster.ATE = mean( daggagg$tau ),
                  avg.sd.cluster.ATE = sqrt( mean( daggagg$sd_tau^2 ) ),
                  mean.total.indiv.ATE = mean( dd$Y1 - dd$Y0 ),
                  sd.total.indiv.ATE = sd( dd$Y1 - dd$Y0 ),
                  cor.size.impact.district = size.impact.cor.district,
                  cor.size.intercept.district = size.intercept.cor,
                  cor.proptx.impact.district = proptx.impact.cor.district,
                  cor.nbar.J = nbar.cor.J,
                  cor.n.tau = corntau
  )
  invisible(stats)
}


#' Build set of districts with cluster sizes and proportion treated,
#' etc.
#'
#' @param sd.district Standard deviation of the district-average ATEs
#'   (i.e., the target cluster-average ATE for each district).  The
#'   district average ATE is what is passed to the code to generate
#'   the clusters for the district as the "ATE."
#' @param cv.school This is the coefficient of variation of cluster
#'   sizes.
#' @param cv.district This is the coefficient of variation of the
#'   number of clusters in each district, across districts.
#' @param cv.nbar This is the coefficient of variation of how much
#'   variation there is in the cluster sizes within district, across
#'   district (i.e., how much average cluster size `nbar` varies
#'   across district)
#' @param empirical If TRUE generate district variation to get as
#'   close as possible to the set targets.  In particular, if the
#'   sd.district = 0.1, then the standard deviation of the actual generated
#'   ATEs will equal 0.1 exactly if empirical = TRUE.  I.e. TRUE means
#'   try to get the target correlations more closely by connecting
#'   things with a calibrated connection with the given correlation
#'   (correlation will still not be exact).
#' @param variable.p Variable proportion treated per district, uniform
#'   around the p.bar value, reaching 2/3 the distance from 0 or 1
#'   (whichever is closer).
#'
#' @return dataframe, one row per district.
build_districts <- function( K, J.bar, n.bar, p.bar=0.5,
                             ATE = 0, sd.district = 0,
                             ICC.district = 0,
                             variable.p = FALSE,
                             cv.district = 0.6, 
                             cv.nbar = 0,
                             min.size = 1,
                             size.impact.cor = 0,
                             nbar.J.cor = 0,
                             proptx.impact.cor = 0,
                             size.intercept.cor = 0,
                             R2 = NULL,
                             num.cov.district = 0 + !is.null(R2),
                             
                             # school-level params
                             sd.school = 0, ICC.school = 0.20,
                             R2.school = NULL,
                             num.cov.school = 0 + !is.null(R2.school),
                             
                             # correlation params at school level
                             cv.school = 0.6,
                             size.impact.cor.school = 0,
                             
                             # Individual-level params
                             R2.indiv = NULL,
                             num.cov.indiv = 0 + !is.null(R2.indiv),
                             
                             empirical = FALSE ) {
  
  
  if ( min.size > 4 ) {
    warning( "min.size > 4: please note that #tx and #co units >= min.size/2, which might surprise you.")
  }
  stopifnot("size.impact.cor must be in (-1,1)" = 
              size.impact.cor >= -1 & size.impact.cor <= 1)
  stopifnot("proptx.impact.cor must be in (-1,1)" = 
              proptx.impact.cor >= -1 & proptx.impact.cor <= 1)
  stopifnot("size.intercept.cor must be in (-1,1)" = 
              size.intercept.cor >= -1 & size.intercept.cor <= 1)
  stopifnot("nbar.J.cor must be in (-1,1)" = 
              nbar.J.cor >= -1 & nbar.J.cor <= 1)
  stopifnot("ICC.district must be in [0,1]" = 
              ICC.district >= 0 & ICC.district <= 1)
  
  
  # district J values
  if (cv.district > 0) {
    num_clusters <- cluster_distn( K = K, 
                                   J.bar = J.bar, 
                                   cv=cv.district,
                                   min.size = min.size )
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
  #  - note: ensure var of district-level effects is exactly sd.district^2, ICC
  grid <- NA
  if ( K == 1 ) {
    grid = tibble( did = 1, district.tau = ATE, district.alpha = 0 )  
  } else {
    grid <- tibble(
      did = 1:K,
      # district.tau = rnorm(K, mean=ATE, sd=sd.district))
      district.tau = MASS::mvrnorm(K, mu=ATE, Sigma=sd.district^2, 
                                   empirical=empirical) %>% 
        as.vector(),
      district.alpha = MASS::mvrnorm(K, mu=0, Sigma=ICC.district,
                                     empirical=empirical) %>% 
        as.vector())
  }
  
  # generate correlations between district.tau and num_clusters & prop_tx_clusters
  #  - uses some magic hacks; idea is to keep order of district.tau,
  #    and to reorder num_clusters & prop_tx_clusters, to avoid automatically
  #    correlating the two of them
  if (size.impact.cor != 0 || (cv.district>0 && empirical) ) {
    stopifnot("Cannot correlate district sizes and average impacts if sizes do not vary" = 
                cv.district > 0)
    num_clusters <- induce_cor( x = grid$district.tau,
                                y = num_clusters,
                                cor = size.impact.cor,
                                empirical = empirical )
  }
  
  if (proptx.impact.cor != 0 || (variable.p && empirical)) {
    stopifnot("Cannot correlate district proportions treated and average impacts if proportions treated do not vary" = 
                variable.p)
    prop_tx_clusters <- induce_cor(x = grid$district.tau,
                                   y = prop_tx_clusters,
                                   cor = proptx.impact.cor,
                                   empirical = empirical )
    # p_temp <- prop_tx_clusters[order(
    #   proptx.impact.cor^2 * prop_tx_clusters + 
    #     (1-proptx.impact.cor^2) * rnorm(K, sd=sd(prop_tx_clusters)) )]
    # prop_tx_clusters <- p_temp[ order(order(sign(proptx.impact.cor)*grid$district.tau)) ]
  }
  if (size.intercept.cor != 0 || (ICC.district!=0 && empirical)) {
    stopifnot("Cannot correlate district sizes and intercepts if ICC.district = 0" = 
                (ICC.district != 0))
    grid$district.alpha <- induce_cor(x = num_clusters,
                                      y = grid$district.alpha,
                                      cor = size.intercept.cor,
                                      empirical = empirical )
  }
  
  #if num_clusters is less than 2 force it to be two 
  grid$J <- pmax(num_clusters, 2)

  
  
  # Make the proportions correspond to number of clusters
  grid$p = round( prop_tx_clusters * grid$J ) / grid$J
  
  
  # Make average cluster sizes for districts
  if (cv.nbar > 0 ) {
    n.bar <- cluster_distn(K, n.bar, cv=cv.nbar)
    if ( nbar.J.cor != 0 ) {
      n.bar <- induce_cor( x = num_clusters,
                           y = n.bar,
                           cor = nbar.J.cor,
                           empirical = empirical )
    }
  } else {
    stopifnot( "cannot correlate number of clusters and average cluster size if average cluster sizes (cv.nbar) do not vary" = 
                 (nbar.J.cor==0) )
  }
  grid$n.bar = n.bar
  
  
  # make sure eventual randomization has at least min.size/2 treated/control clusters
  #  in each district, by applying floor and ceiling to p, depending on J
  #   - i.e., force 2 <= J*p <= J-2
  #print("this is the value of min.size")
  #print(min.size)
  stopifnot( "All districts should have at least max( min.size, 2 ) clusters" = all( grid$J >= max( 2, min.size) ) )
  ms2 = pmax( 1, round( min.size/2 ) )
  grid <- grid %>% 
    mutate(p = ifelse(J*p < ms2, ms2/J, p),
           p = ifelse(J*p > J-ms2, (J-ms2)/J, p))
  
  # Add district level covariates if needed.
  if ( num.cov.district > 0 ) {
    if ( is.null( R2 ) ) {
      R2 = 0
    }
    cov_block <- generate_covariate_block( grid$district.alpha,
                                           ICC.district, R2, k = num.cov.district, 
                                           var_prefix = "Xd" )
    grid = bind_cols( grid, cov_block )
  }
  
  
  # Add in the school level stuff
  grid$sd.school = sd.school
  grid$ICC.school = ICC.school
  grid$num.cov.school = num.cov.school
  grid$R2.school = R2.school
  if ( num.cov.school > 0 && is.null( R2.school ) ) {
    grid$R2.school = 0
  }
  grid$cv.school = cv.school
  grid$size.impact.cor.school = size.impact.cor.school
  
  # Add in the individual level stuff
  grid$R2.indiv = R2.indiv
  grid$num.cov.indiv = num.cov.indiv
  
  return( grid )
}



## testing code ----
if ( FALSE ) {
  grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.district = 0.1,
                           ICC.district = 0.2,
                           variable.p = TRUE, 
                           num.cov.district = 3 )
  
  grid
}



#### Code to build data given a description of the districts ####

generate_data_for_districts <- function(
    # datafrane if district-level characteristics
  grid,
  ICC.district = var( grid$district.alpha ) ) {
  
  if ( nrow( grid ) == 1 ) {
    ICC.district = 0
  }
  
  stopifnot( all( ICC.district + grid$ICC.school < 1 ) )
  
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
        finite.model = FALSE,
        size.impact.correlate = sign( size.impact.cor.school ),
        correlate.strength = size.impact.cor.school,
        size.ratio = cv.school * 0.6,
        proptx.impact.correlate = FALSE,
        zero.corr = TRUE,
        varY0 = 1 - ICC.district,
        cluster.rand = TRUE,
        min.size = 1
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
  
  
  sim.data
}


## Testing code ----
if (FALSE) {
  
  grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                           ATE = 0.2, sd.district = 0.1,
                           ICC.district = 0.2,
                           variable.p = TRUE, 
                           num.cov.district = 3,
                           
                           # correlation params at cluster level
                           cv = 0.6, 
                           size.impact.cor = 0.2,
                           proptx.impact.cor = 0.2,
                           size.intercept.cor = 0.2, 
                           
                           # school-level params
                           sd.school = 0.5, 
                           ICC.school = 0.7,

                           # correlation params at individual level
                           cv.school = 0.1,
                           size.impact.cor.school = 0.75 )
  
  grid
  
  # Make data from grid
  df <- generate_data_for_districts( grid )
  
  df
  
  describe_clusterRCT(Yobs ~ Z | sid | did, df)
  
  df %>% 
    group_by(did) %>% 
    summarize(mn = mean(Y0), n=n()) %>% 
    summarize(cor = cor(mn, n))
}




#### Code to simulate full data set ####


#' Generate a cluster-randomized trial, with multiple districts.  This
#' will build a grid of the districts (stats on each of the
#' districts), and then fill in the grid with individual schools and
#' students.
#' 
generate_districts <- function(
    # district-level params
  K, J.bar, p.bar=0.5,
  ATE.district, sd.district,
  ICC.district,
  variable.p = FALSE,
  R2.district = NULL,
  num.cov.district = 0 + !is.null(R2.district),
  
  # variation in number of clusters and average cluster size
  cv.district = 0.6,
  cv.nbar = 0,
  min.size = 1,
  
  # correlation params at district level
  size.impact.cor.district = 0,
  proptx.impact.cor.district = 0,
  size.intercept.cor.district = 0,
  
  # school-level params
  n.bar, sd.school,
  ICC.school = 0.20,
  R2.school = NULL,
  num.cov.school = 0 + !is.null(R2.school),
  
  # correlation params at school level
  cv.school = 0.6,
  size.impact.cor.school = 0.75,
  
  # Individual-level params
  R2.indiv = NULL,
  num.cov.indiv = 0 + !is.null(R2.indiv),
  
  empirical = FALSE ) {
  
  stopifnot( ICC.district + ICC.school < 1 )
  
  grid <- build_districts( K = K, J.bar = J.bar, n.bar = n.bar, 
                           p.bar = p.bar,
                           ATE = ATE.district, sd.district = sd.district,
                           ICC.district = ICC.district,
                           variable.p = variable.p, 
                           cv.district = cv.district,
                           cv.nbar = cv.nbar,
                           size.impact.cor = size.impact.cor.district,
                           proptx.impact.cor = proptx.impact.cor.district,
                           size.intercept.cor = size.intercept.cor.district,
                           R2 = R2.district,
                           num.cov.district = num.cov.district,
                           min.size = min.size,
                           
                           # school-level params
                           sd.school = sd.school, ICC.school = ICC.school,
                           R2.school = R2.school,
                           num.cov.school = num.cov.school,
                           cv.school = cv.school,
                           size.impact.cor.school = size.impact.cor.school,
                           
                           # Individual-level params
                           R2.indiv = R2.indiv,
                           num.cov.indiv = num.cov.indiv,
                           
                           empirical = empirical )
  
  
  sim.data <- generate_data_for_districts( grid, 
                                           ICC.district = ICC.district )
  
  attr(sim.data, "district_key") <- grid
  
  sim.data
}



if (FALSE) {
  
  set.seed( 40404 )
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
    size.impact.cor.school = 0.75
  )
  df
  
  describe_clusterRCT(Yobs ~ Z | sid | did, df)
  
  df %>% 
    group_by(did) %>% 
    summarize(mn = mean(Y0), n=n()) %>% 
    summarize(cor = cor(mn, n))
}

