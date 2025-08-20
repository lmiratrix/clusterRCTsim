


# Simple simulation with massive amounts of treatment variation.
#
# This sim has equal sized blocks and clusters.
#
# GOAL: The goal is to assess whether any of the SE estimators are
# targeting different things in the face of lots of randomness at the
# block and cluster level.

# This simulation (when varation=FALSE) really looks at how the
# different estimators react to variation differently, when the blocks
# and clusters are all the same size.
#
# In particular, it appears FE methods and FI methods give _different_ SEs, and the FE methods are NOT calibrated to any of the following models:
#
# 1) Fixed dataset, treatment assignment is the only randomness.
# 2) Fixed district, but clusters and individuals sampled randomly
# 3) Full superpopulation with districts, then clusters, then students
#    all sampled randomly.


# Notes
#
# sd.district = sd.school = 5
# FE over inflated by 40%, FI nails it for framing (2), above.
#
# sd.district = sd.school = 0
# All the SEs are about the same, and the calibration 
# is about 1 for everything (around 0.98 or 0.99).
# (There was substantial ICC values.)
#
# sd.district = 0, sd.school = 5
# The estimated SEs _all_ align with the superpopulation
# SE to high accuracy.  The fixed district is also quite close.
#
# sd.district = 5, sd.school = 0
# Here we see large difference in the _estimated_ SEs for FE (being
# much too high) compared to the FI (which are calibrated perfectly to
# the fixed district for (2)).
# calib.s = 0.4
# calib.m = calib.f = 4.2
#
# sd.district = 5, sd.school = 0
# and NO baseline district level variation (ICC.district = 0)
# This was exactly the same as above, with baseline variation.


source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_district.R" ) )


# This will run a single scenario!



if ( FALSE ) {
    grid
    describe_district_grid( grid )
}



one_run <- function( variation = FALSE ) {
    
    ICCd = 0
    ICCs = 0.4
    sdd = 5
    sds = 0
    cv = 0
    if ( variation ) {
        cv = 0.8
    }
    
    grid <- build_districts( K = 10, J.bar = 6, n.bar = 20, p.bar = 0.5,
                             ATE = 0.5, 
                             sd.district = sdd,
                             sd.school = sds,
                             variable.p = FALSE,
                             cv.district = cv,
                             cv.nbar = cv,
                             cv.school = cv,
                             min.size = 4,
                             ICC.district = ICCd,
                             ICC.school = ICCs, 
                             size.intercept.cor = 0,
                             size.impact.cor.school = 0,
                             R2 = NULL, R2.school = NULL, R2.indiv = NULL,
                             num.cov.district = 0,
                             num.cov.school = 0,
                             num.cov.indiv = 0,
                             empirical = FALSE )
    
    
    rps = map_df( 1:3, \(x) {
        rr = single.CRT.fixed.district( grid = grid,
                                        n.runs = 3, ICC.district = ICCd )
        rr
    }, .id = "midrun" )
    
    rps
}

if ( FALSE ) {
    # testing
    one_run()
}


library( furrr )
cli::cli_h1( "Starting multisession plan for simulation" )
plan(multisession, workers = parallel::detectCores() - 1 )
tictoc::tic()
R = 500
rps = future_map_dfr( 1:R, ~ one_run( FALSE ), 
                      .id = "run",
                      .options = furrr_options(seed = TRUE),
                      .progress=TRUE )

tictoc::tic()
plan(sequential)
cli::cli_h1( "Closed multisession plan for simulation" )

# Analyze the simulation results ----

rps
nrow(rps) / length( unique(rps$method ) )

source( here::here( "R/aggregate_simulation_results.R" ) )

#debug( calculate_performance_statistics )
rps$midID <- paste0( rps$run, ".", rps$midrun )
res1 <- calculate_performance_statistics( rps, runID = "midID" )
resB <- calculate_performance_statistics( rps, runID = "run" )
res1$SE.s / resB$SE.s
res1$SE.f / resB$SE.f
res1$SE.s / res1$SE.f
resB$SE.s / resB$SE.f
res1$SE.m = resB$SE.f

res2 <- rps %>%
    group_by( method ) %>%
    summarise( Etau = mean( ATE.finite.student ),
               sdtau = sd( ATE.finite.student ),
               pdef = mean( !is.na( SE.hat ) ) ) 

res <- left_join( res1, res2, by = "method" ) %>%
    rename( EATE = mean.ATE.hat,
            ESE_hat = mean.SE.hat ) %>%
    mutate( calib.s = sqrt( ESE_hat^2 / SE.s^2 ),
            calib.m = sqrt( ESE_hat^2 / SE.m^2 ),
            calib.f = sqrt( ESE_hat^2 / SE.f^2 )
    ) %>%
    dplyr::select( -starts_with( "MCSE" ), -cover )

names(res)            

res$SE.s / res$SE.m
res$SE.s / res$SE.f

res <- res %>%
    mutate( is_FI = str_detect(method,"FI" ) )

res %>%
    dplyr::select( -weight, -pdef, -avg.df, -p.noSE, -is_FI ) %>%
    arrange( ESE_hat ) %>%
    relocate( method, ESE_hat, SE.f, SE.m, SE.s ) %>%
    #    dplyr::select( method:SE.s, starts_with("calib" ) ) %>%
    knitr::kable(digits=2) %>%
    print()



se_low = mean( res$SE.s ) * 0.9
se_high = mean( res$SE.s ) * 1.1
se_low
se_high

ggplot( res, aes( y=ESE_hat, x=SE.s, pch=weight ) ) +
    #    facet_wrap( ~ weight ) +
    facet_wrap( ~ is_FI, labeller = label_both ) +
    geom_point( aes( col="super" ) ) +
    geom_point( aes( x = SE.f, col="finite" ) ) + 
    geom_point( aes( x = SE.m, col="fix district" ) ) + 
    coord_fixed() +
    expand_limits( y=c( se_low, se_high ), 
                   x=c( se_low, se_high ) ) +
    theme_minimal() +
    geom_abline( slope = 1, intercept = 0, col="red" )  +
    labs( x = "True standard error", y = "Average estimated SE" )



