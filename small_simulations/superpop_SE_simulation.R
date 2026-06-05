

# superpop_SE_simulation.R
#
# Compares superpopulation vs. finite-population true standard errors
# across scenarios, holding blocks fixed and varying cluster-level
# sampling. Companion to reports/clusterRCT_SEs.Rmd.
#
# Active simulations: controlled by WHICH_SIMS (currently 3 and 8).
# Output: saved .rds files in results/ + kable performance tables.
# Status: Runs with N_REPS = 5 as of 2026-06-04

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_block.R" ) )

library( lme4 )
library( lmerTest )
library( future )
library( here )

WHICH_SIMS <- c() # Run none of the simulations
WHICH_SIMS <- c( 3, 8 )
N_REPS <- 5  # Set low for testing; increase for a full run

set.seed( 404404 )
grid <- build_blocks( K = 20, J.bar = 10, n.bar = 30,   # build_districts -> build_blocks
                      p.bar = 0.5,
                      ATE = 0.2, sd.block = 0.1,          # sd.ATE -> sd.block
                      ICC.block = 0.2,                     # ICC -> ICC.block
                      cv.block = 0.6,                      # variable.J=TRUE -> cv.block
                      variable.p = TRUE,
                      cv.nbar = 0.33,                      # variable.nbar=TRUE -> cv.nbar
                      num.cov.block = 3,                   # num.cov.district -> num.cov.block
                      size.intercept.cor = -0.8,
                      size.impact.cor = -0.75,
                      proptx.impact.cor = 1,
                      nbar.J.cor = -0.9,
                      sd.cluster = 0.5,                    # sd.school -> sd.cluster
                      ICC.cluster = 0.2,                   # ICC.school -> ICC.cluster
                      cv.cluster = 0.6,                    # size.ratio.school=0.6 -> cv.cluster
                      size.impact.cor.cluster = 0.75,      # size.impact.cor.school
                      empirical = TRUE )


# Explore our grid and initial DGP ----
if ( FALSE ) {
  
  
  library(GGally)
  grid %>%
    dplyr::select( district.tau, district.alpha, J, p, n.bar ) %>%
    ggpairs( lower = list(continuous = wrap("smooth", method = "lm", se=FALSE)) ) +
    theme_minimal()
  
  grid
  
  
  # Test data to viz the context
  dat = generate_data_for_districts( grid,
                                     variable.n = TRUE,
                                     sd.school = 0.5, ICC.school = 0.2,
                                     size.ratio.school = 0.6,
                                     size.impact.cor.school = 0.75 )
  dat
  
  describe.data(dat)
  
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
    ICC.district = 0.5,
    
    # school-level params
    variable.n = TRUE,
    sd.school = 0.5, ICC.school = 0.2,
    size.ratio.school = 0.6,
    size.impact.cor.school = 0.75,
    
    # method params
    n.runs = 3,
    R = 2,
    .progress="none",
    parallel = FALSE )
  sim
}





# Simulation 1! ----


if ( 1 %in% WHICH_SIMS ) {
  sim <- run_sim( "", grid=grid, R=500, n.runs=3 )
}




# Demo analysis of simulation 1 results ----

if ( FALSE ) {
  sim <- readRDS( "results/superpop_SE_simulation.rds" )
  
  length( table( sim$method ) )
  
  #sim <- filter( sim, 
  #               disfavored == 0 ) %>%
  #  dplyr::select( -disfavored, -ATE.finite.dist.school, -ATE.finite.dist.student ) %>%
  #  as_tibble()
  
  sim
  #table( sim$run, sim$subrun )
  
  ATE.school = mean( sim$ATE.finite.school )
  ATE.student = mean( sim$ATE.finite.student )
  
  ggplot( sim, aes( method, ATE.hat ) ) +
    geom_hline( yintercept = c( ATE.school, ATE.student ), color="red" ) +
    geom_boxplot() +
    coord_flip() +
    theme_minimal() +
    labs( x = "", title = "Distribution of point estimates" )
  
  
  source( "aggregate_simulation_results.R" )  # Possible junk code according to Claude's audit

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



# Simulation 2 (No correlation between district-level intercepts and sizes) ----

if ( 2 %in% WHICH_SIMS ) {
  # Possible junk code according to Claude's audit (old build_districts API, sim not in WHICH_SIMS)
  set.seed( 4434343 )
  grid2 = build_districts(  K = 20, J.bar = 10, n.bar = 30,
                            p.bar = 0.5,
                            ATE = 0.2, sd.ATE = 0.1,
                            ICC = 0.2,
                            variable.J = TRUE,
                            variable.p = TRUE,
                            variable.nbar = TRUE,
                            num.cov.district = 3,
                            size.intercept.cor = 0,
                            size.impact.cor = 0,
                            proptx.impact.cor = 0,
                            nbar.J.cor = 0,

                            # school-level params
                            variable.n = TRUE,
                            sd.school = 0.5, ICC.school = 0.2,
                            size.ratio.school = 0.6,
                            size.impact.cor.school = 0,
                            empirical = TRUE )
  run_sim( 2, grid=grid2, R=500, n.runs = 3 )
}


# Simulation 3 (need to get more systematic) ----

if ( 3 %in% WHICH_SIMS ) {
  set.seed( 434343 )
  grid3 = build_blocks( K = 4, J.bar = 7, n.bar = 50,   # build_districts -> build_blocks
                        p.bar = 0.7,
                        ATE = 0.2, sd.block = 0.2,        # sd.ATE -> sd.block
                        sd.cluster = 0.5,                  # sd.school -> sd.cluster
                        ICC.block = 0.2,                   # ICC -> ICC.block
                        ICC.cluster = 0.2,                 # ICC.school -> ICC.cluster
                        cv.block = 0.6,                    # variable.J=TRUE -> cv.block
                        variable.p = TRUE,
                        cv.nbar = 0.33,                    # variable.nbar=TRUE -> cv.nbar
                        cv.cluster = 1.0,                  # size.ratio.school=0.6 -> cv.cluster
                        min.size = 4,
                        size.intercept.cor = -0.7,
                        size.impact.cor = 1,
                        proptx.impact.cor = 1,
                        nbar.J.cor = 1,
                        size.impact.cor.cluster = -0.8,    # size.impact.cor.school
                        empirical = TRUE )
  grid3 <- arrange( grid3, J )
  grid3$J
  grid3$p <- c( 3, 3, 3, 3 ) / grid3$J
  grid3$block.tau                                          # district.tau -> block.tau

  describe_block_grid( grid3 )                             # describe_district_grid -> describe_block_grid

  sim3 <- run_sim( 3, grid = grid3, R = N_REPS, n.runs = 3, parallel = FALSE )
}


## Some testing code to explore simulation 3 ----
if ( FALSE ) {
  
  sim3 %>% group_by( method ) %>%
    filter( disfavored == 0 ) %>%
    summarise( sd_df = sd( df ),
               df = mean( df ) ) 
  
  # some testing to look at df for this scenario
  grid3 = readRDS( "results/superpop_SE_simulation_grid3.rds" )
  describe_district_grid( grid3 )
  
  dt = single.CRT.fixed.district( grid3, just.return.data = TRUE )
  dt
  
  compare_methods( Yobs ~ Z | sid | did, data = dt ) %>%
    filter( disfavored == 0 ) %>%
    knitr::kable( digits = 2 )
  
  debug( linear_model_estimators)
  linear_model_estimators( Yobs ~ Z | sid | did, data = dt ) %>%
    knitr::kable( digits = 2 )
  
  
  debug( MLM_estimators )
  MLM_estimators( Yobs ~ Z | sid | did, data = dt ) %>%
    knitr::kable( digits = 2 )
  
  describe_clusterRCT( Yobs ~ Z | sid | did, data = dt)
  
  ss <- single.CRT.fixed.district( grid3, n.runs = 3 )
  filter( ss, subrun == 1 ) %>%
    dplyr::select( -contains( "dist." ) ) %>%
    filter( disfavored == 0 ) %>%
    knitr::kable( digits=2)
}



# Simulation 4 (no treatment variation other than district level) ----

# This simulation should have finite and superpop standard errors that
# are basically the same

if ( 4 %in% WHICH_SIMS ) {
  # Possible junk code according to Claude's audit (old sd.school API, sim not in WHICH_SIMS)
  set.seed( 34364343 )

  grid4 = readRDS( "results/superpop_SE_simulation_grid3.rds")
  grid4$sd.school = 0

  #options( warn = 2 )
  run_sim( 4, grid = grid4, R = 500, n.runs = 3, parallel=TRUE )
}



# Simulation 5 (massive treatment variation at the school level) ----

# This simulation should have finite and superpop standard errors that
# are basically the same

if ( 5 %in% WHICH_SIMS ) {
  # Possible junk code according to Claude's audit (old sd.school API, sim not in WHICH_SIMS)
  set.seed( 34343143 )
  grid5 = readRDS( "results/superpop_SE_simulation_grid3.rds")
  grid5$sd.school = 1

  run_sim( "5", grid = grid5, 500, 3, parallel=TRUE )
}



# Simulation 6 (nothing unusual at all) ----

if ( 6 %in% WHICH_SIMS ) {
  # Possible junk code according to Claude's audit (old build_districts API, sim not in WHICH_SIMS)
  set.seed( 343431433 )

  grid6 = build_districts(  K = 5, J.bar = 6, n.bar = 40,
                            p.bar = 0.5,
                            ATE = 0.2, sd.ATE = 0.1,
                            sd.school = 0.2,

                            ICC = 0.2,
                            ICC.school = 0.2,

                            variable.J = TRUE,
                            variable.p = FALSE,
                            variable.nbar = TRUE,
                            variable.n = TRUE,

                            size.ratio = 0.4,
                            min.size = 4,

                            size.ratio.school = 0.1,

                            #size.intercept.cor = -0.7,
                            #size.impact.cor = 1,
                            #proptx.impact.cor = 0,
                            #nbar.J.cor = 0,
                            #size.impact.cor.school = -0.8,
                            empirical = TRUE )

  describe_district_grid( grid6 )

  run_sim( "6", grid6, 500, 3, parallel=TRUE )
}


# Simulation 7 (making big superpop SEs vs finite) ----

if ( 7 %in% WHICH_SIMS ) {
  # Possible junk code according to Claude's audit (old build_districts API, sim not in WHICH_SIMS)
  set.seed( 344314333 )

  grid7 = build_districts(  K = 1, J.bar = 9, n.bar = 500,
                            p.bar = 7/9,
                            ATE = 0.2, sd.ATE = 0,
                            sd.school = 0.3,  # tau

                            ICC.school = 0 )

  describe_district_grid( grid7 )
  
  
  set.seed( 344314333 )
  
  
  sim7 <- run_sim( "7", grid = grid7, R=500, n.runs = 3, parallel=TRUE )
  
  
  library(MASS)
  
  
  generate_clusterCRT( 3, 5, 0.2, 0.3^2, -1, 0.5 )  
  
  sim8 <- run_sim( "8", 
                   dgp_function = generate_clusterCRT,
                   dgp_args = list( J=9, n = 500, ICC = 0.0031, tau=0.25, cor.ICC.tau=-1, p=7/9 ),
                   R = 500, n.runs = 3, parallel = TRUE )
  
  
  # One more, for luck, with $ICC = 0.05$, $\tau=0.3^2$, $J = 16$, $n=100$, $p = 3/4$ and $\rho = -0.7$:
  sim9 <- run_sim( "9", 
                   dgp_function = generate_clusterCRT,
                   dgp_args = list( J=16, n = 100, ICC = 0.05, tau=0.3^2, cor.ICC.tau=-0.7, p=3/4 ),
                   R = 400, n.runs = 4, parallel = TRUE )
}

if ( FALSE ) {
  source( here::here( "aggregate_simulation_results.R" ) )
  sim7 = readRDS( "results/superpop_SE_simulation7.rds" )
  sim7
  
  aa <- calculate_performance_statistics( sim7, align = TRUE ) %>%
    arrange( weight ) %>%
    dplyr::select( -bias, -MCSE.bias ) %>%
    mutate( ratio = SE.s / SE.f )
  aa
  
  
  sim7 %>% group_by( method ) %>%
    filter( disfavored == 0 ) %>%
    summarise( sd_df = sd( df ),
               df = mean( df ) ) 
  
  grid7
  
  # Explore data from our grid
  dt = single.CRT.fixed.district( grid7, just.return.data = TRUE )
  dt
  lmer( Y1 ~ Z + (1|sid), data=dt )
  agg <- dt %>% group_by( sid ) %>%
    summarise( Y1 = mean( Y1 ),
               Y0 = mean( Y0 ),
               n = n() )
  sd( agg$Y0 )
  sd( agg$Y1 )
  describe_clusterRCT( Yobs ~ Z | sid | did, data = dt)
  sd( dt$Y0 )
  sd( dt$Yobs )
  sd( dt$Y1 )
  
  compare_methods( Yobs ~ Z | sid, data=dt )
  
  grid7
  ICC = 0.13^2
  sqrt( ICC / (grid$p * (1-grid$p)*grid$J) + (1-ICC) / ( grid7$J * grid7$n.bar * grid7$p * (1-grid7$p) ) )
  
  ICCtau = 0.3^2
  tx_side = ICCtau / (grid$p*grid$J) + (1-ICCtau) / ( grid7$J * grid7$n.bar * grid7$p )
  co_side = 1 / ( grid7$J * grid7$n.bar * (1-grid7$p) )
  sqrt( tx_side + co_side )
  
  
  rns  = single.CRT.fixed.district( grid7,  n.runs = 10, .progress = TRUE )
  rns %>% group_by( method ) %>%
    filter( disfavored == 0 ) %>%
    summarise( SE = sd( ATE.hat ) )
  
  rns$run = 1
  calculate_performance_statistics( rns )
  
  
  compare_methods( Yobs ~ Z | sid, data = dt ) %>%
    filter( disfavored == 0 ) %>%
    knitr::kable( digits = 2 )
  
  
  
  debug( linear_model_estimators)
  linear_model_estimators( Yobs ~ Z | sid | did, data = dt ) %>%
    knitr::kable( digits = 2 )
  
  
  debug( MLM_estimators )
  MLM_estimators( Yobs ~ Z | sid | did, data = dt ) %>%
    knitr::kable( digits = 2 )
  
  describe_clusterRCT( Yobs ~ Z | sid | did, data = dt)
  
  ss <- single.CRT.fixed.district( grid, n.runs = 3 )
  filter( ss, subrun == 1 ) %>%
    dplyr::select( -contains( "dist." ) ) %>%
    filter( disfavored == 0 ) %>%
    knitr::kable( digits=2)
  
  
}



# Simulation 8 (sim 3, but limited cluster variation) ----


if ( 8 %in% WHICH_SIMS ) {
  set.seed( 434343 )
  grid3 = build_blocks( K = 4, J.bar = 7, n.bar = 50,   # build_districts -> build_blocks
                        p.bar = 0.5,
                        ATE = 0.2, sd.block = 0.2,        # sd.ATE -> sd.block
                        sd.cluster = 0.5,                  # sd.school -> sd.cluster
                        ICC.block = 0.2,                   # ICC -> ICC.block
                        ICC.cluster = 0.2,                 # ICC.school -> ICC.cluster
                        cv.block = 0.6,                    # variable.J=TRUE -> cv.block
                        variable.p = TRUE,
                        cv.nbar = 0.33,                    # variable.nbar=TRUE -> cv.nbar
                        cv.cluster = 0.167,                # size.ratio.school=0.1 -> 0.1/0.6
                        min.size = 4,
                        empirical = TRUE )
  describe_block_grid( grid3 )                             # describe_district_grid -> describe_block_grid
  sim8 <- run_sim( 8, grid = grid3, R = N_REPS, n.runs = 3, parallel = FALSE )
}


