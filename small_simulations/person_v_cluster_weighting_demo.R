
# person_v_cluster_weighting_demo.R
#
# Demonstrates a thought experiment: in a study with mixed large/small
# clusters, targeting the person-weighted ATE effectively halves the
# number of informative clusters vs. the cluster-weighted ATE.
# Examines how SE scales with ICC across estimators.
#
# Output: line plot of SE vs. ICC by method.
# Status: Runs with N_REPS = 5 as of 2026-06-04

# Simulation to implement following thought experiment


#@Luke MMiratrix I can’t think of how to cleanly explain what is going on here.  I might understand.

# Imagine a study with 40 clusters:
#   10 T of size 1000, 
# 10 C of size 1000, 
# 10 T of size 20, 
# 10 C of size 20.
# 

# When we target the effect for the average cluster, we have 40
# cluster RCT, with some small and big clusters.  Weighting the small
# clusters the same as the big clusters hurts us a little bit, but we
# have 40 observations, and unless the clusters are REALLY small, this
# doesn’t hurt us much at all.
#
# When we target the effect for the average person, we basically are
# dropping 20 clusters.  So any benefit we might have gotten by only
# caring about the big clusters is basically offset by what is now
# basically a 20 cluster trial.
#
# Thoughts?
# 

library( simhelpers )

source( here::here( "small_simulations/simulation_DGP.R" ) )

N_REPS <- 5  # Set low for testing; increase for a full run



my_dgp <- function( ICC, mpl = model.params.list ) {
  source( here::here( "small_simulations/simulation_DGP.R" ) )
  
  p = mpl
  p$ICC.2 = ICC
  
  p$K = 1
  p$J = 40
  sim.data = make_two_tier_PUMP_data( n1 = 20, ATE1 = 0,
                                      n2 = 100, ATE2 = 0,
                                      params = p,
                                      blocks = FALSE )
  sim.data
}



my_estimator <- function( sim.data ) {
  
  blocks = n_distinct(sim.data$D.id) > 1
  
  # Analyze the data
  form = Yobs ~ Tx | S.id | D.id
  if ( !blocks ) {
    form = Yobs ~ Tx | S.id
  }
  control_form = NULL
  include_covariates = TRUE
  if ( include_covariates ) {
    if ( blocks ) {
      control_form = ~ X.jk + C.ijk
    } else {
      control_form = ~ C.ijk
    }
  }
  
  c1 <- clusterRCT::compare_methods( form, data=sim.data,
                                     control_formula = control_form,
                                     include_DB = TRUE,
                                     include_MLM = FALSE,
                                     include_LM = TRUE,
                                     include_agg = TRUE,
                                     include_dumb = FALSE,
                                     include_method_characteristics = FALSE )
  
  
  
  if ( !blocks ) {
    sim.data$D.id = "sing"
  }
  
  
  # Add in the true ATEs
  ds = sim.data %>%
    group_by( D.id, S.id ) %>%
    summarise( tauk = mean( Y1 - Y0 ),
               n = n(), .groups = "drop" ) %>%
    group_by( D.id ) %>%
    summarise( tau_c = mean( tauk ),
               tau_p = weighted.mean( tauk, w=n ),
               J = n(),
               n = sum( n ) ) %>%
    summarise( tau_cb = mean( tau_c ),
               tau_c = weighted.mean( tau_c, w=J ),
               tau_pb = mean( tau_p ),
               tau_p = weighted.mean( tau_p, w=n ),
               J = sum( J ),
               n = sum( n ) )
  
  c1$tau_c = ds$tau_c
  c1$tau_p = ds$tau_p
  #c1$tau_cb = ds$tau_cb
  #c1$tau_pb = ds$tau_pb
  c1$J = ds$J
  c1$n = ds$n
  c1$sdY0 = sd( sim.data$Y0 )
  
  # Return the results
  c1 %>%
    dplyr::select( -any_of( c( "weight", "df" ) ) )  
}



if ( FALSE ) {
  dat = my_dgp( 0.1, mpl = model.params.list )
  res = my_estimator(dat)
  res
  #   2 * 10 * (20 + 1000 )
}



sim_func = bundle_sim(my_dgp, my_estimator, calc_performance )

args( sim_func )
sim_func( 3, ICC = 0.1 )



ICCs = c( 0, 0.05, 0.1, 0.2, 0.4 )
setup_parallel_if_needed()

res <- ICCs %>%
  set_names() %>%
  future_map_dfr( sim_func, reps = N_REPS,
                  .id = "ICC", .progress=TRUE,
                  .options = furrr_options( seed=TRUE ) )

res

res$ICC = as.numeric( res$ICC )
res
mc = clusterRCT::method_characteristics()
res <- left_join( res, mc, by="method" )

ggplot( res, aes( ICC, SE, col=method, lty=weight )) +
  geom_line() +
  scale_y_continuous( limits = c(0, 0.5 ) ) +
  labs( x = "ICC", y = "Estimated ATE", col = "Method" ) +
  theme_minimal() +
  theme( legend.position = "bottom" )




