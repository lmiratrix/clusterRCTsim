
# Simulation to investigate the weighting properties of the MLM
# approaches

source( here::here( "small_simulations/simulation_DGP.R" ) )


#### Mike's ICC Simulation ####

# Inspired by Mike Weiss questions on MLM

# Mike wanted to see how the MLM estimators targeted different
# estimands based on the ICC


run_sim_MLM <- function( ICC, blocks = FALSE, ..., R = 100 ) {
  cat( "ICC ", ICC, "\n" )
  
  rps = map_df( 1:R, \(r) one_run( ICC.2 = ICC, blocks = blocks, ...,
                                add_canonical_weights = FALSE ),
                .id = "runID", .progress=TRUE )
  
  head( rps )
  
  
  res <- rps %>%
    group_by( method ) %>%
    summarise( EATE = mean( ATE_hat ),
               tau_p = mean( tau_p ),
               tau_c = mean( tau_c ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )
  
  
  res
}

if ( FALSE ) {
  one_run( ICC.2 = 0.2, blocks=FALSE )
  
  debug( one_run )
  one_run( ICC.2 = 0.2, blocks=TRUE )
  
  run_sim_MLM( 0.2, R=4 )
}


if ( FALSE ) {
  

  
  #one_run( blocks=FALSE, ICC.2 = 0.8 )
  
  ICCs = c( 0, 0.2, 0.4, 0.6, 0.8 )
  setup_parallel_if_needed()
  
  res <- ICCs %>%
    set_names() %>%
    future_map_dfr( run_sim_MLM, R = 100, .id = "ICC", .progress=TRUE )
  
  res
  
  # Make families of estimators for grouping on the plot
  res$meth_group = str_extract( res$method, "^[^_-]+" )
  table( res$meth_group )
  
  # Note how all methods are level w.r.t. the ICC, except MLM which
  # goes from person-weighted to cluster-weighted as ICC increases.
  ggplot( res, aes( ICC, EATE ) ) +
    facet_wrap( ~ meth_group, nrow=1 ) +
    scale_y_continuous( limits = c(0, 0.5 ) ) +
    geom_line( aes( group=method )) +
    geom_hline( aes( yintercept = mean( tau_c ), lty="Cluster" ), col="black" ) +
    geom_hline( aes( yintercept = mean( tau_p ), lty="Person" ), col="black" )
  
  
}





#### Mike's ICC Simulation extended to blocks ####

if ( FALSE ) {
  library( clusterRCT)
  # Our first simulation has all the blocks the same size and shape.
  dat <- one_run( ICC.2 = 0.4, blocks = TRUE,
           add_canonical_weights = FALSE, return_data = TRUE )
  describe_clusterRCT( Yobs ~ Tx | S.id | D.id, data=dat )
  
  
  model.params.spec = model.params.list
  model.params.spec$omega.3 = 0
  model.params.spec$ICC.3 = 0.2
  model.params.spec$nbar = 20
  model.params.spec$K = 6
  
  mm = model.params.spec[ !map_lgl( model.params.spec, is.null ) ]
  as.data.frame(mm)
  
  ICCs = c( 0, 0.025, 0.05, 0.1, 0.2, 0.4 )
  setup_parallel_if_needed()
  
  res <- ICCs %>%
    set_names() %>%
    future_map_dfr( run_sim_MLM, R = 100, blocks = TRUE,
                    model.params = model.params.spec,
                    .id = "ICC", .progress=TRUE )
  
  res
  
  # Make families of estimators for grouping on the plot
  res$meth_group = str_extract( res$method, "^[^_-]+" )
  table( res$meth_group )
  
  # Note how all methods are level w.r.t. the ICC, except MLM which
  # goes from person-weighted to cluster-weighted as ICC increases.
  ggplot( res, aes( ICC, EATE ) ) +
    facet_wrap( ~ meth_group, nrow=1 ) +
    #facet_wrap( ~ method ) +
    scale_y_continuous( limits = c(0, 0.5 ) ) +
    geom_line( aes( group=method )) +
    geom_hline( aes( yintercept = mean( tau_c ), lty="Cluster" ), col="black" ) +
    geom_hline( aes( yintercept = mean( tau_p ), lty="Person" ), col="black" )
  
  
  
  # And now let's do a messier version ---- 
  one_run( ICC.2 = 0.75, blocks = TRUE,
           add_canonical_weights = FALSE )
  
  ICCs = c( 0, 0.2, 0.4, 0.6, 0.75 )
  setup_parallel_if_needed()
  
  res <- ICCs %>%
    set_names() %>%
    future_map_dfr( run_sim_MLM, R = 100, blocks = TRUE, 
                    .id = "ICC", .progress=TRUE )
  
  res
  
  if ( FALSE ) {
    # For debugging
    res <- ICCs %>%
      set_names() %>%
      map_dfr( run_sim_MLM, R = 100, blocks = TRUE, 
                      .id = "ICC", .progress=TRUE )
  }
  
  # Make families of estimators for grouping on the plot
  res$meth_group = str_extract( res$method, "^[^_-]+" )
  table( res$meth_group )
  
  # Note how all methods are level w.r.t. the ICC, except MLM which
  # goes from person-weighted to cluster-weighted as ICC increases.
  ggplot( res, aes( ICC, EATE ) ) +
    facet_wrap( ~ meth_group, nrow=1 ) +
    scale_y_continuous( limits = c(0, 0.5 ) ) +
    geom_line( aes( group=method )) +
    geom_hline( aes( yintercept = mean( tau_c ), lty="Cluster" ), col="black" ) +
    geom_hline( aes( yintercept = mean( tau_p ), lty="Person" ), col="black" )
  
  
}






# Checking ICC is working ----

# Run a sim to make sure we can recover a DGP parameter of the ICC.
if ( FALSE ) {
  
  check_it <- function( ICC ) {
    p = model.params.list
    p$ICC.2 = ICC
    p$K = 10
    p$omega.2 = 0
    p$omega.3 = 0
    dat = make_two_tier_PUMP_data( 10, 0.2, params=p, blocks = TRUE )
    m <- lmer( Yobs ~ 1 + Tx + (1|S.id) + (1|D.id), data=dat )
    
    c( ICC = ICC, 
       fixef( m ), 
       S.id = VarCorr(m)$S.id,
       D.id = VarCorr(m)$D.id,
       sigma2 = sigma(m)^2,
       ICC_hat = ICC1,
       tau = mean(dat$Y1 - dat$Y0 )
    )
  }
  
  rr = map( 1:40, ~ {
    ICC1 = check_it( 0 )
    ICC2 = check_it( 0.7 )
    res <- bind_rows( ICC1, ICC2 )
    res
  } )
  
  rr = bind_rows( rr )
  rr
  rr %>% 
    mutate( sd2Y = (S.id+D.id+sigma2),
            ICC_hat = S.id / sd2Y ) %>%
    group_by( ICC ) %>%
    summarise( mean_tx = mean( Tx ),
               SE = sd( Tx ),
               SE_Etx = sd( Tx ) / sqrt(n()),
               meanS = mean(S.id),
               meanD = mean(D.id),
               meanSigma = mean(sigma2),
               meansd2Y = mean(sd2Y),
               mean_ICC = mean(ICC_hat),
               sd_ICC = sd(ICC_hat),
               tau_p = mean(tau),
               sd_t = sd(tau) )

  
  
}




