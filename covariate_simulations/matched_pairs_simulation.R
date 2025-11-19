

# Simulation of matched pairs experiments

# How much do individual level covariates help?


set.seed( 1039 )
library( tidyverse )
library( PUMP )
library( clusterRCT )




#### Making matched pairs data ####

source( "R/data_generation_code.R" )

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




if ( FALSE ) {

  mpl = model.params.list
  mpl$nbar = 100
  mpl$J = 20
  mpl$R2.1 = 0.8
  mpl$R2.2 = 0
  sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", mpl, Tbar = 0.5, include_POs = TRUE )
  
  # We are explaining variation
  summary( lm( Yobs ~ C.ijk, data=sim.data ) )
  
  # We are explaining nothing as we should
  summary( lm( Yobs ~ X.jk, data=sim.data ) )
  
 
  
  # Check that tilting works (making clusters vary in size, correlated
  # with treatment effect)
  model.params.list$omega.2 = 1
  sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- tilt_balanced_PUMP_data(sim.data, tilt=0.5, corr=0.90 )
  
  sizes = table( sim.data$S.id )
  qplot( sizes )
  
  describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=sim.data,
                       control_formula = ~ C.ijk + X.jk )

  class( sim.data )
  sd( sim.data$Yobs )
  # ATE in effect size units
  mean( sim.data$Y1 - sim.data$Y0 ) / sd( sim.data$Yobs )
  model.params.list$MDES
  head( sim.data )  
  
  agg <- sim.data %>%
    group_by( S.id, D.id ) %>%
    summarise( tau = mean( Y1 - Y0 ),
               Y0 = mean( Y0 ),
               n = n(),
               .groups="drop" ) 
  sd( agg$tau ) / sd( sim.data$Yobs )
  cor( agg$n, agg$tau )
  cor( agg$n, agg$Y0 )
  
  ggplot( agg, aes( n, tau ) ) +
    geom_point() +
    geom_smooth( method="lm", se=FALSE )
  ggplot( agg, aes( n, Y0 ) ) +
    geom_point() +
    geom_smooth( method="lm", se=FALSE )
  
  compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                       control_formula = ~ C.ijk + X.jk )
  
}




#### Simulation code ####



one_run <- function( double_cluster = FALSE ) {
  
  sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
  sim.data <- tilt_balanced_PUMP_data(sim.data, 
                                      tilt = 0.7, 
                                      corr = 0.8 )
  
  # Make a double cluster
  if ( double_cluster ) {
    # Reassign clusters in 1st D to 2nd.
    # (We assume the S.id are all unique)
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





