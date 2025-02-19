
# Illustrate making and using grids in our simulation framework




source( "R/data_generation_code.R" )
source( "R/simulation_support_code_fixed_district.R" )

library( lme4 )
library( lmerTest )
library( future )
library( clusterRCT )

set.seed( 404404 )
grid <- build_districts( K = 20, J.bar = 10, n.bar = 30, 
                         p.bar = 0.5,
                         ATE = 0.2, sd.ATE = 0.1,
                         ICC = 0.2,
                         variable.J = TRUE, 
                         variable.p = TRUE, 
                         variable.nbar = TRUE,
                         num.cov.district = 3, 
                         size.intercept.cor = -0.8,
                         size.impact.cor = -0.75, 
                         proptx.impact.cor = 1, 
                         nbar.J.cor = -0.9, 
                         
                         # school-level params
                         variable.n = TRUE,
                         sd.school = 0.5, ICC.school = 0.2,
                         size.ratio.school = 0.6,
                         size.impact.cor.school = 0.75,
                         empirical = TRUE )

grid


desc <- describe_district_grid( grid )

desc$cv.J
desc$cv.n


# Explore our grid and initial DGP ----
if ( FALSE ) {
  
  
  library(GGally)
  grid %>%
    dplyr::select( district.tau, district.alpha, J, p, n.bar ) %>%
    ggpairs( lower = list(continuous = wrap("smooth", method = "lm", se=FALSE)) ) +
    theme_minimal()
  
  grid
  
  
  # Test data to viz the context
  dat = generate_data_for_districts( grid )
  dat
  
  describe.data(dat)
  
  describe_clusterRCT( Yobs ~ Z | sid | did, data=dat )
  
  dagg = dat %>%
    group_by( sid, did, Z ) %>%
    summarise( tau = mean( Y1 - Y0 ),
               n = n(), .groups = "drop" ) %>%
    group_by( did ) %>%
    mutate( J = n() ) %>%
    ungroup() 
  dagg
  
  dagg$did = fct_reorder( factor(dagg$did), dagg$J )
  ggplot( dagg, aes( n ) ) +
    facet_wrap( ~ did ) +
    geom_dotplot( method="histodot")
  
  ggplot( dagg, aes( J, n ) ) +
    geom_point( aes(  color=as.factor(did) ) )  +
    geom_smooth( se=FALSE, method="lm" )
  
  ggplot( dagg, aes( J, tau ) ) +
    geom_point( aes(  color=as.factor(did) ) )  +
    geom_smooth( se=FALSE, method="lm" )
}


# Run short sim for testing ----

if ( FALSE ) {
  sim <- run.scenario.fixed.district(     
    grid,
    
    # method params
    n.runs = 3,
    R = 2,
    .progress="none",
    parallel = FALSE )
  sim
  
  table( sim$subrun)
}





# running simulation and saving it

if ( FALSE ) {
  
  # running the simulation
  sim <- run_sim( "test", grid=grid, R=100, n.runs=3 )
  sim
  
  # Demo analysis of simulation results
  
  sim <- readRDS( "results/superpop_SE_simulationtest.rds" )
  
  length( table( sim$method ) )
  
  sim
  
  ATE.school = mean( sim$ATE.finite.school )
  ATE.student = mean( sim$ATE.finite.student )
  
  ggplot( sim, aes( method, ATE.hat ) ) +
    geom_hline( yintercept = c( ATE.school, ATE.student ), color="red" ) +
    geom_boxplot() +
    coord_flip() +
    theme_minimal() +
    labs( x = "", title = "Distribution of point estimates" )
  
  
  source( "R/aggregate_simulation_results.R" )
  
  rs <- calculate_performance_statistics( sim, align = TRUE ) %>%
    arrange( weight )
  rs
  
  rsL <- rs %>%
    pivot_longer( c( "SE.s", "bias", "RMSE" ) )
  rsL
  
  ggplot( rsL, aes( x=method, y=value, colour = weight ) ) +
    facet_wrap( ~ name, nrow = 1 ) +
    geom_point() +
    geom_hline( yintercept=0, lty="dashed" ) +
    coord_flip() +
    theme_minimal()
  
  
  
  # Estimated standard errors for individual linear regression (and
  # variants) are terrible in this scenario???
  rs %>%
    group_by( method ) %>%
    summarise( SErat.s = mean.SE.hat / SE.s,
               SErat.f = mean.SE.hat / SE.f,
               SE_fs = SE.s / SE.f ) %>%
    arrange( SErat.s )
  
  
}



