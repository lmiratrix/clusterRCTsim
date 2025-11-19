
library( tidyverse )

source( here::here( "R/data_generation_code.R"))

#' Randomly drop units within cluster so the cluster sizes vary.
#'
#' @param corr If nonzero, will drop units to induce correlation with
#'   the cluster-level impact or control outcome mean and the cluster
#'   size.
#' @param tilt 1 - Minimum proportion of units to retain in each
#'   cluster. Tilt of 0 means don't subsample.  Tilt of 1 means
#'   maximum.
tilt_balanced_PUMP_data <- function( data, tilt=0.3, corr=0, byY0=TRUE ) {
  
  stopifnot( tilt >= 0 && tilt <= 1 )
  tilt = 1 - tilt
  
  stats <- data %>% 
    group_by( S.id, D.id ) %>%
    summarise( Ybar0 = mean( Y0 ),
               tau = mean( Y1 - Y0 ),
               n = n(), 
               .groups = "drop" )
  
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
#' The strenght of the covariates is controlled by R2.2 and R2.1,
#' which are the proportions of variance explained at the cluster and
#' individual level.  The combined covariate C is a weighted sum of
#' the individual-level covariate X and the cluster-level covariate W,
#' where C = X + scale_W * W.
#'
#' X explains only individual level variation.  W explains only level
#' 2 variation.
#'
#' @param param_list A list with the following elements:
#'   - K: Number of blocks
#'   - J: Number of clusters
#'   - nbar: Average cluster size
#'   - MDES: Minimum detectable effect size (treatment effect)
#'   - ICC.3: Block-level ICC
#'   - ICC.2: Cluster-level ICC
#'   - R2.2: Proportion of cluster-level variance explained by true W covariate.
#'   - R2.1: Proportion of individual-level variance explained by true X covariate.
#' @param scale_W Scale parameter for the combined covariate C, where
#'   C = X + scale_W * W.
#'
#' @return A data frame with the generated data.
#' 
make_covariate_data <- function( param_list, scale_W = 0.5 ) {
  
  stopifnot( c( "K", "J", "nbar", "MDES", 
                "ICC.3", "ICC.2", "R2.2", "R2.1" ) %in% names( param_list ) )
  
  gamma.k = rnorm( param_list$K, 
                   0, sqrt( param_list$ICC.3 ) )
  alpha.jk = rnorm( param_list$J*param_list$K, 
                    0, sqrt( (1-param_list$R2.2)*param_list$ICC.2 ) )
  e.ijk = rnorm( param_list$J * param_list$K * param_list$nbar,
                 0, 
                 sqrt( (1-param_list$R2.1)*(1 - param_list$ICC.2 - param_list$ICC.3 ) ) )
  
  W = rnorm( param_list$J * param_list$K, 0, 1 )
  X = rnorm( param_list$J * param_list$K * param_list$nbar, 0, 1 )
  
  beta_W = sqrt( param_list$R2.2*param_list$ICC.2 )
  
  beta_X = sqrt( param_list$R2.1*(1 - param_list$ICC.2 - param_list$ICC.3 ) )
  
  clust = tibble( S.id = 1:(param_list$J*param_list$K),
                  D.id = rep( 1:param_list$K, each=param_list$J ),
                  gamma.k = gamma.k[D.id],
                  alpha.jk = alpha.jk,
                  W = W ) %>%
    expand_grid( i = 1:param_list$nbar )
  
  clust <- clust %>%
    mutate( X = X,
            Y0 = gamma.k + alpha.jk + beta_W*W + beta_X*X + e.ijk,
            Y1 = Y0 + param_list$MDES,
            Z = randomizr::block_and_cluster_ra( blocks=D.id,
                                                 clusters=S.id,
                                                 prob=0.5),
            C = X + scale_W * W,
            Yobs = ifelse( Z, Y1, Y0 ) )
  
  clust <- clust %>%
    dplyr::select( -gamma.k, -alpha.jk, -i )
  
  clust 
}




if ( FALSE ) {
  
  pl = list( K = 20,
             J = 20,
             nbar = 20,
             MDES = 0.2,
             ICC.3 = 0.05,
             ICC.2 = 0.15,
             R2.2 = 0.5,
             R2.1 = 0.3 )
  
  dat = make_covariate_data( pl, scale_W = 0.5 )
  dat
  
  library( lme4 )
  M <- lmer( Yobs ~ 1 + Z + W + X + (1|D.id) + (1|S.id), data=dat )
  arm::display(M)
  
  clusterRCT::describe_clusterRCT( Yobs ~ Z | S.id | D.id, 
                                   data=dat,
                                   control_formula = ~ W + X )  
  
  
  dat <- tilt_balanced_PUMP_data( dat, tilt=0.3, corr=0.5, byY0=TRUE )
  
  clusterRCT::describe_clusterRCT( Yobs ~ Z | S.id | D.id, 
                                   data=dat,
                                   control_formula = ~ W + X )
  
  mean( dat$Y1 - dat$Y0 )
}



if ( FALSE ) {
  
  model.params.list <- list(
    J = 4                          # number of schools
    , K = 20                          # number of districts
    , nbar = 20                       # number of individuals per school
    , MDES = 0.2            # minimum detectable effect size
    , ICC.3 = 0.2             # district intraclass correlation
    , R2.2 = 0.7              
    , ICC.2 = 0.2             # school intraclass correlation
    , R2.1 = 0.8    
  )
  
  model.params.list$J = 2
  
  #sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- make_covariate_data( model.params.list, scale_W = 0.5 )
  
  table( dat$D.id )
  table( dat$S.id )
  table( dat$D.id, dat$Z )
  
  clusterRCT::describe_clusterRCT( Yobs ~ Z | S.id | D.id, 
                                   data=sim.data,
                                   control_formula = ~ C )
  
  sim.data <- tilt_balanced_PUMP_data(sim.data, 0.8, 0.8 )
  clusterRCT::describe_clusterRCT( Yobs ~ Z | S.id | D.id, 
                                   data=sim.data,
                                   control_formula = ~ C )
  
}




#  Final DGP code ----



gen_sim_data <- function( exp_type = c( "matched pairs", "double cluster", "quad" ),
                          X_on = TRUE,
                          W_on = TRUE,
                          K = NULL,
                          R2.2 = 0.8,
                          R2.1 = 0.8,
                          ICC.2 = 0.1,
                          ICC.3 = 0.1 ) {
  
  
  # exp_type = match.arg( exp_type )
  
  model.params.list <- list(
    J = 4                          # number of schools
    , K = 20                          # number of districts
    , nbar = 20                       # number of individuals per school
    , MDES = 0.2            # minimum detectable effect size
    , ICC.3 = ICC.3             # district intraclass correlation
    , R2.2 = R2.2              
    , ICC.2 = ICC.2             # school intraclass correlation
    , R2.1 = R2.1   
  )
  
  
  if (!X_on ) {
    model.params.list$R2.1 = 0.0
  }
  if ( !W_on ) {
    model.params.list$R2.2 = 0.0
  }
  
  if ( exp_type == "matched pairs" ) {
    model.params.list$J = 2
  } else if ( exp_type == "double cluster" ) {
    model.params.list$J = 2
  } else if ( exp_type == "quad" ) {
    model.params.list$J = 4
    model.params.list$K = 3
  }
  
  if ( !is.null( K ) ) {
    model.params.list$K = K
  }
  
  #sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- make_covariate_data( model.params.list, scale_W = 0.5 )
  sim.data <- tilt_balanced_PUMP_data(sim.data, 0.3, 0.8 )
  
  # Make a double cluster
  if ( exp_type == "double cluster" ) {
    sim.data$D.id[ sim.data$D.id == 1 ] <- 2
  }
  
  sim.data
}


