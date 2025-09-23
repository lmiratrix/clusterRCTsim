


# Simple simulation with massive amounts of treatment variation.
# Variant: No blocking, just a single CRT


source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_district.R" ) )


# This will run a single scenario!



if ( FALSE ) {
    grid
    describe_district_grid( grid )
}



one_run <- function( variation = FALSE ) {
    
    ICCd = 0.4
    grid <- build_districts( K = 1, J.bar = 40, n.bar = 20, p.bar = 0.5,
                                 ATE = 0.5, 
                                 sd.district = 0,
                                 sd.school = 2,
                                 variable.p = FALSE,
                                 cv.district = 0,
                                 cv.nbar = 0,
                                 cv.school = 0,
                                 min.size = 4,
                                 ICC.district = 0,
                                 ICC.school = 0.4, 
                                 size.intercept.cor = 0,
                                 size.impact.cor.school = 0,
                                 R2 = NULL,
                                 R2.school = NULL,
                                 R2.indiv = NULL,
                                 num.cov.district = 0,
                                 num.cov.school = 0,
                                 num.cov.indiv = 0,
                                 empirical = FALSE )
    grid
    
    
    rps = map_df( 1:3, \(x) {
        rr = single.CRT.fixed.district( grid = grid,
                                        n.runs = 3, ICC.district = 0,
                                        include_disfavored = TRUE )
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
R = 100
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

ggplot( res, aes( y=ESE_hat, x=SE.s, pch=weight ) ) +
    #    facet_wrap( ~ weight ) +
    facet_wrap( ~ is_FI ) +
    geom_point( aes( col="super" ) ) +
    geom_point( aes( x = SE.f, col="finite" ), alpha=0.5 ) + 
    geom_point( aes( x = SE.m, col="fix district" ), alpha=0.25 ) + 
    coord_fixed() +
    expand_limits( x = c( se_low, se_high ), y = c( se_low, se_high ) ) +
    theme_minimal() +
    geom_abline( slope = 1, intercept = 0, col="red" )  +
    labs( x = "True standard error", y = "Average estimated SE" )



