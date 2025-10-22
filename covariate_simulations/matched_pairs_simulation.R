

# Simulation of matched pairs experiments

# How much do individual level covariates help?


set.seed( 1039 )
library( tidyverse )
library( PUMP )
library( clusterRCT )


ATE = 0.125
tx_het = FALSE

model.params.list <- list(
  M = 1                            # number of outcomes
  , J = 2                          # number of schools
  , K = 20                          # number of districts
  , nbar = 20                       # number of individuals per school
  , S.id = NULL                     # N-length vector of school assignments
  , D.id = NULL                     # N-length vector of district assignments
  , Xi0 = 0                         # scalar grand mean outcome under no treatment
  , MDES = ATE            # minimum detectable effect size
  , R2.3 = 0.0            
  , ICC.3 = 0.2             # district intraclass correlation
  , omega.3 = 0.2*tx_het           # ratio of district effect size variability
  , R2.2 = 0.5              
  , ICC.2 = 0.2             # school intraclass correlation
  , omega.2 = 0.0*tx_het          # ratio of school effect size variability
  , R2.1 = 0.8    
)


source( "R/data_generation_code.R" )
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



one_run <- function( double_cluster = FALSE ) {
  
  sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- tilt_balanced_PUMP_data(sim.data, 0.3, 0.8 )
  
  # Make a double cluster
  if ( double_cluster ) {
    sim.data$D.id[ sim.data$D.id == 1 ] <- 2
  }
  
  sim.data <- sim.data %>% 
    group_by( D.id, S.id ) %>%
    mutate( X = C.ijk + X.jk,
            W = mean( X ) ) %>%
    ungroup()
  
  c1 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                     include_MLM = FALSE,
                                     include_LM = TRUE,
                                     include_agg = FALSE,
                                     include_method_characteristics = FALSE )
  
  c2 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                     control_formula = ~ X,
                                     include_MLM = FALSE,
                                     include_LM = TRUE,
                                     include_agg = FALSE,
                                     include_method_characteristics = FALSE )
  
  c3 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                     control_formula = ~ W,
                                     include_MLM = FALSE,
                                     include_LM = TRUE,
                                     include_agg = FALSE,
                                     include_method_characteristics = FALSE )
  
  c1 = bind_rows( raw=c1, X = c2, W = c3, .id="adjustment" )
  c1$tau = mean( sim.data$Y1 - sim.data$Y0 )
  c1$sdY0 = sd( sim.data$Y0 )
  
  c1
}


if ( FALSE ) {
  # testing
  one_run()
}

R = 100
rps = map_df( 1:R, ~ one_run(), .id = "runID", .progress=TRUE )

head( rps )



rps %>%
  group_by( method, adjustment ) %>%
  summarise( SE = sd( ATE_hat ), .groups = "drop" ) %>%
  pivot_wider( names_from=adjustment, values_from=SE )  %>%
  arrange( -X )


  mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )




res <- rps %>%
  group_by( method, adjustment ) %>%
  summarise( Etau = mean( tau ),
             sdtau = sd( tau ),
             EATE = mean( ATE_hat ),
             SE = sd( ATE_hat ),
             ESE_hat = sqrt( mean( SE_hat^2 ) ),
             bias = mean( EATE - tau ),
             SEb = sd(ATE_hat-tau) / sqrt(n()),
             bias_t = bias / SEb, .groups = "drop" ) %>%
  mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )


res %>%
  arrange( calib ) %>%
  knitr::kable(digits=2)





