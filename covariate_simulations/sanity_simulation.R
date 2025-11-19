
# Sanity check simulations on covariate adjustment being effective


library( tidyverse )
library( clusterRCT )


source( here::here( "covariate_simulations/DGP_code_varying_cluster_sizes.R" ) )
source( here::here( "covariate_simulations/sim_analysis_function.R" ) )


if ( FALSE ) {
  dta = gen_sim_data( exp_type = "quad", K=10 )
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


# Quick single sanity simulation ----





if ( FALSE ) {
  
  rps = map_dfr( 1:100, \(x) {
    dta = gen_sim_data( exp_type = "quad", K=10,
                        X_on = TRUE,
                        W_on = FALSE,
                        R2.1 = 0.9, R2.2 = 0.9 )
    #describe_clusterRCT( Yobs ~ Z | S.id | D.id, data= dta,
    #                     control_formula = ~ X + W )
    
    res = analyze_sim_data( dta, inc_MLM = FALSE )
    res
  }, .progress=TRUE )

  dta = gen_sim_data( exp_type = "quad", K=10,
                      X_on = TRUE,
                      W_on = FALSE,
                      R2.1 = 0.9, R2.2 = 0.9 )
  describe_clusterRCT( Yobs ~ Z | S.id | D.id, data= dta,
                       control_formula = ~ X + W )
  
  SE_calc( params = list( R2.1 = 0.9, R2.2 = 0, ICC.2 = 0.2, ICC.1 = 0.2 ),
           dta = dta )
  SE_calc( params = list( R2.1 = 0, R2.2 = 0, ICC.2 = 0.2, ICC.1 = 0.2 ),
           dta = dta )
  
  head(rps)

  rps %>% 
    group_by( adjustment, method ) %>%
    summarise( SE = sd( ATE_hat, na.rm=TRUE ), .groups = "drop" ) %>%
    pivot_wider( names_from=adjustment, values_from=SE ) %>%
    mutate( red = X / raw )
}

# Run the simulations ----

params = expand.grid(
  exp_type = c( "quad" ),
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


plan( multisession, workers = parallel::detectCores() - 1 )

tictoc::tic()
objs <- mget(ls(), envir = .GlobalEnv)
params$res = future_pmap( params,
                          R2.1 = 0.9,
                          R2.2 = 0.9,
                          ICC.3 = 0,
                          ICC.2 = 0.025,
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
           here::here( "covariate_simulations/results/sanity_simulation_results.rds" ) )






# Analyze the results ----

rps = read_rds( here::here( "covariate_simulations/results/sanity_simulation_results.rds" ) )

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
  facet_grid( X_on  ~ W_on, labeller= label_both ) +
  geom_point( position=position_dodge(width=0.5) ) +
  geom_point( data=raws, color="black" ) +
  geom_hline( yintercept = 0 ) +
  coord_flip() +
  theme_minimal()





