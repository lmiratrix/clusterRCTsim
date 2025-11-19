

# Simulation of matched pairs experiments

# How much do individual level covariates help?


library( tidyverse )
library( clusterRCT )


source( here::here( "covariate_simulations/DGP_code_varying_cluster_sizes.R" ) )


# Simulation functions ----


gen_sim_data <- function( exp_type = c( "matched pairs", "double cluster", "quad" ),
                          X_noise = 0.05 ) {
  
  # exp_type = match.arg( exp_type )
  
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
  
  if ( exp_type == "matched pairs" ) {
    model.params.list$J = 2
  } else if ( exp_type == "double cluster" ) {
    model.params.list$J = 2
  } else if ( exp_type == "quad" ) {
    model.params.list$J = 4
    model.params.list$K = 3
  }
  
  #sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- make_covariate_data( model.params.list, scale_W = 0.5 )
  sim.data <- tilt_balanced_PUMP_data(sim.data, 0.3, 0.8 )
  
  # Make a double cluster
  if ( exp_type == "double cluster" ) {
    sim.data$D.id[ sim.data$D.id == 1 ] <- 2
  }
  
  if (TRUE ) {
    sim.data <- sim.data %>% 
      group_by( D.id, S.id ) %>%
      mutate( X = W + rnorm( n(), 0, X_noise ) ) %>%
      ungroup()
  } else {
    stop( "Unfinished" )
    # for seeing how adjustment of aggregate vs individual gives different SEs
    sim.data <- sim.data %>% 
      group_by( D.id, S.id ) %>%
      mutate( X = X + W,
              W = mean( X ) ) %>%
      ungroup()
  }
  
  sim.data
}


if ( FALSE ) {
  a = gen_sim_data( "matched pairs", X_noise = 0.5 )
  cor( a$X, a$W )
  a = gen_sim_data( "matched pairs", X_noise = 0.05 )
  cor( a$X, a$W )  
}



analyze_sim_data <- function( sim.data, inc_MLM = FALSE ) {
  
  c1 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = FALSE,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  c2 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                     control_formula = ~ X,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = FALSE,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  c3 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                     control_formula = ~ W,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = FALSE,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  c1 = bind_rows( raw=c1, X = c2, W = c3, .id="adjustment" )
  c1$tau = mean( sim.data$Y1 - sim.data$Y0 )
  c1$sdY0 = sd( sim.data$Y0 )
  
  c1
}


library( simhelpers )

sim_function <- bundle_sim( gen_sim_data, analyze_sim_data )
args( sim_function )



# Run the simulations ----

params = expand.grid(
  exp_type = c( "matched pairs", "double cluster", "quad" ),
  X_noise = c( 0.05, 0.5 )
)
params = bind_rows( params, params, params, params )
params$seed = 24244 + 1:nrow(params) * 13
params$reps = 100
params
params = slice_sample( params, n = nrow(params) )
params

library( furrr )

if ( FALSE ) {
  cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "FORK")
  plan(cluster, workers = cl )
  parallel::clusterEvalQ(cl, {
    source( "R/data_generation_code.R" )
    source( "DGP_code.R" )
    NULL
  })
}
plan( multisession, workers = parallel::detectCores() - 1 )

tictoc::tic()
objs <- mget(ls(), envir = .GlobalEnv)
names(objs)
params$res = future_pmap( params,
                          .f = sim_function,
                          .options = furrr_options(seed = NULL,
                                                   packages = c( "tidyverse", "clusterRCT" ),
                                                   globals = objs), .progress = TRUE )
tictoc::tic()
plan(sequential)

params <- as_tibble(params)
params

rps <- params %>%
  unnest( cols = c(res) )
head( rps )

write_rds( rps, "covariate_exploration_simulation_results.rds" )

# Analyze the results ----

rps = read_rds( "covariate_exploration_simulation_results.rds" )

res <- rps %>%
  group_by( exp_type, X_noise, method, adjustment ) %>%
  summarise( SE = sd( ATE_hat ), .groups = "drop" ) %>%
  pivot_wider( names_from=adjustment, values_from=SE )

res
table( res$method )

res <- rps %>%
  group_by( exp_type, X_noise, method, adjustment ) %>%
  summarise(## Etau = mean( tau ),
    # sdtau = sd( tau ),
    EATE = mean( ATE_hat ),
    SE = sd( ATE_hat ),
    ESE_hat = sqrt( mean( SE_hat^2 ) ),
    sd_SEhat = sd( SE_hat ),
    #bias = mean( EATE - tau ),
    #SEb = sd(ATE_hat-tau) / sqrt(n()),
    #bias_t = bias / SEb, 
    .groups = "drop" ) %>%
  mutate( calib = sqrt( ESE_hat^2 / SE^2 ) ) 


res %>%
  #arrange( calib ) %>%
  arrange( method, adjustment ) %>%
  knitr::kable(digits=2) %>%
  print()




res
ggplot( res, aes( method, SE, pch = as.factor(X_noise), col=adjustment ) ) +
  #facet_wrap( ~ X_noise ) +
  geom_point( position=position_dodge(width=0.5) ) +
  geom_hline( yintercept = 0 ) +
  coord_flip()


res$SE[ is.na(res$SE) ] = 0
ggplot( res, aes( adjustment, SE ) ) +
  facet_grid( X_noise ~ exp_type ) +
  geom_jitter( width = 0.25 )


res$calibNA = is.na(res$calib)
res$calib %na<-% 0.6
ggplot( res, aes( adjustment, calib, col = calibNA ) ) +
  facet_grid( X_noise ~ exp_type ) +
  geom_jitter( width = 0.25 ) +
  geom_hline( yintercept = 1, col="red" ) +
  theme_minimal()

library( roperators )
res$sd_SEhat %na<-% 0
res$sdSEna = res$sd_SEhat == 0
ggplot( res, aes( adjustment, sd_SEhat, col=sdSEna ) ) +
  facet_grid( X_noise ~ exp_type ) +
  geom_jitter( width = 0.25 )




