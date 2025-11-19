

# Simulation to understand cluster and individual covariates better.


library( tidyverse )
library( clusterRCT )


source( here::here( "covariate_simulations/DGP_code_varying_cluster_sizes.R" ) )
source( here::here( "covariate_simulations/sim_analysis_function.R" ) )

# Test/illustrate the simulation functions ----


if ( FALSE ) {
  a = gen_sim_data( "matched pairs", X_on = FALSE, W_on = TRUE )
  cor( a$X, a$W )
  
  clusterRCT::describe_clusterRCT( Yobs ~ Z | S.id | D.id, data= a,
                                   control_formula = ~ C )
  
  # TODO: Getitng all tx and all co blocks.  Fix!
  
}


# The analysis function for the simulations ----




if ( FALSE ) {
  dta = gen_sim_data( "quad" )
  analyze_sim_data( dta, inc_MLM = FALSE )
  
  describe_clusterRCT( Yobs ~ Z | S.id | D.id, data= dta,
                         control_formula = ~ C )
  
  sim.data <- dta %>%
    group_by( D.id, S.id ) %>%
    mutate( C_mn = mean(C),
            C_gc = C - C_mn ) %>%
    ungroup()
  
  clusterRCT::compare_methods( Yobs ~ Z | S.id | D.id, data=sim.data,
                               control_formula = ~ C_mn + C_gc,
                               include_MLM = FALSE,
                               include_LM = TRUE,
                               include_agg = TRUE,
                               include_gee = FALSE,
                               include_method_characteristics = FALSE )
  
}



library( simhelpers )

sim_function <- bundle_sim( gen_sim_data, analyze_sim_data )
args( sim_function )



# Run the simulations ----

params = expand.grid(
  exp_type = c( "matched pairs", "double cluster", "quad" ),
  X_on = c( TRUE, FALSE ),
  W_on = c( TRUE, FALSE )
)

params = bind_rows( params, params, params, params )
params = bind_rows( params, params )
cli::cli_inform( "We have {nrow(params)} jobs" )
params$seed = 24244 + 1:nrow(params) * 13
params$reps = 50
params
params = slice_sample( params, n = nrow(params) )
params <- as_tibble(params)
params
cli::cli_inform( "We have {sum(params$reps)} iterations, total" )

library( furrr )

if ( FALSE ) {
  # Not working code---very sadness stuff
  # > Trying to get the parallel code to load things right
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
                                                   globals = objs),
                          .progress = TRUE )
tictoc::tic()
plan(sequential)

params

rps <- params %>%
  unnest( cols = c(res) )
rps

write_rds( rps, 
           here::here( "covariate_simulations/results/covariate_exploration_simulation_results.rds" ) )






# Analyze the results ----

rps = read_rds( here::here( "covariate_simulations/results/covariate_exploration_simulation_results.rds" ) )

res <- rps %>%
  group_by( exp_type, X_on, W_on, method, adjustment ) %>%
  summarise( SE = sd( ATE_hat, na.rm=TRUE ), .groups = "drop" )

resW <- res %>%
  pivot_wider( names_from=adjustment, values_from=SE )

resW %>%
  knitr::kable( digits = 2 )
table( res$method )

library( roperators )

res$SE %na<-% -0.01

res$SE[ res$SE > 0.5 ] = 0.5


res <- filter( res,
               !str_detect(method, "-db" ) ) %>%
  mutate( method = str_replace( method, "-(crve|het)$", "" ) )

raws = filter( res, adjustment == "raw" )

# Plot by method across the three scenarios and different adjustments
ggplot( res, aes( method, SE, col = exp_type, pch = adjustment ) ) +
  facet_grid( X_on + W_on ~ exp_type, labeller= label_both ) +
  geom_point( position=position_dodge(width=0.5) ) +
  geom_point( data=raws, color="black" ) +
  geom_hline( yintercept = 0 ) +
  coord_flip() +
  theme_minimal()





ggplot( res, aes( adjustment, SE, col = exp_type ) ) +
  facet_wrap( ~ exp_type ) +
  geom_point( position=position_dodge(width=0.5) ) +
  geom_hline( yintercept = 0 ) +
  coord_flip() +
  theme_minimal()



res <- rps %>%
  group_by( exp_type, method, adjustment ) %>%
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




