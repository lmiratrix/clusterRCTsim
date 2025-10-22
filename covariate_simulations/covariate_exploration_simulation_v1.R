

# Simulation of matched pairs experiments

# How much do individual level covariates help?


set.seed( 1039 )
library( tidyverse )
library( PUMP )
library( clusterRCT )


ATE = 0.125
tx_het = FALSE


source( "R/data_generation_code.R" )
source( "DGP_code.R" )


if ( FALSE ) {
  mplt = model.params.list
  mplt
  mplt$nbar = 40
  mplt$J = 30
  mplt$K = 20
  ddd = make_covariate_data( mplt )
  ddd
  mplt$ICC.3
  
  mplt$ICC.2 
  mplt$R2.1
  
  clusterRCT::describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=ddd,
                                   control_formula = ~ W )
  clusterRCT::describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=ddd,
                                   control_formula = ~ Xc + Xbar )
  
  
  data=make_covariate_data( model.params.list ) 
}

if ( FALSE ) {
  
  mpl = model.params.list
  mpl$nbar = 100
  mpl$J = 20
  mpl$R2.1 = 0.8
  mpl$R2.2 = 0
  sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", mpl, Tbar = 0.5, include_POs = TRUE )
  
  
  summary( lm( Yobs ~ C.ijk, data=sim.data ) )
  summary( lm( Yobs ~ X.jk, data=sim.data ) )
  
  
  sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- tilt_balanced_PUMP_data(sim.data, 0.5, 0.9 )
  table( sim.data$S.id )
  describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=sim.data,
                       control_formula = ~ C.ijk + X.jk )
  
  class( sim.data )
  sd( sim.data$Yobs )
  mean( sim.data$Y1 - sim.data$Y0 ) / sd( sim.data$Yobs )
  head( sim.data )  
  describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=sim.data,
                       control_formula = ~ C.ijk + X.jk )
  
  
  head( sim.data )
  
  
  compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                   control_formula = ~ C.ijk + X.jk )
  
}

if ( FALSE ) {


  sim.data <- make_covariate_data( model.params.list, scale_W = 0.5 )
  sim.data <- tilt_balanced_PUMP_data(sim.data, 0.3, 0.8 )
  head( sim.data )
  describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=sim.data,
                       control_formula = ~ X + W )  


  model.params.list$R2.2 = 0.9  
  mpl <- model.params.list

  SE = sqrt( mpl$ICC.2 / (mpl$J*mpl$K*0.5*0.5) +
               (1-mpl$ICC.2-mpl$ICC.3) / (mpl$J*mpl$nbar *mpl$K) )
  SE
  
  SE2 = sqrt( (1-mpl$R2.2)*mpl$ICC.2 / (mpl$J*mpl$K*0.5*0.5) +
               (1-mpl$R2.1)*(1-mpl$ICC.2-mpl$ICC.3) / (mpl$J*mpl$nbar *mpl$K) )
  SE2
  
}

# The Experiments ----

## Matched pairs, noise added to level 2 covariate

model.params.list$J = 2


# Simulation functions ----


gen_sim_data <- function( exp_type = c( "matched pairs", "double cluster", "quad" ),
                     X_noise = 0.05,
                     inc_MLM = TRUE ) {
  
 # exp_type = match.arg( exp_type )
  
  model.params.list <- list(
      J = 4                          # number of schools
    , K = 20                          # number of districts
    , nbar = 20                       # number of individuals per school
    , MDES = 0.2            # minimum detectable effect size
    , ICC.3 = 0.2             # district intraclass correlation
    , R2.2 = 0.5              
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
  gen_sim_data( "matched pairs", 0.5 )
}

analyze_sim_data <- function( sim.data ) {

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


if ( FALSE ) {
  # testing
  one_run()
}


params = expand.grid(
  exp_type = c( "matched pairs", "double cluster", "quad" ),
  X_noise = c( 0.05, 0.5 )
)
params$seed = 24244 + 1:nrow(params) * 13
params$reps = 400
params


library( furrr )
cl <- parallel::makeCluster(parallel::detectCores() - 1, type = "FORK")

plan(cluster, workers = cl )
tictoc::tic()
params$res = future_pmap( params,
                         .f = sim_function,
                         .options = furrr_options(seed = NULL), .progress = TRUE )
tictoc::tic()
plan(sequential)

if ( FALSE ) {
R = 400

library( furrr )
plan(multisession, workers = parallel::detectCores() - 1 )
tictoc::tic()
rps = future_map_dfr(1:R, ~ one_run( double_cluster = TRUE ),
                     .id = "runID", .progress=TRUE,
                     .options = furrr_options(seed = TRUE) )
tictoc::tic()
plan(sequential)
}

head( rps )



rps %>%
  group_by( method, adjustment ) %>%
  summarise( SE = sd( ATE_hat ), .groups = "drop" ) %>%
  pivot_wider( names_from=adjustment, values_from=SE )  %>%
  arrange( -X )




res <- rps %>%
  group_by( method, adjustment ) %>%
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






