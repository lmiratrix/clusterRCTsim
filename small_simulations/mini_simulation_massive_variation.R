
# mini_simulation_massive_variation.R
#
# Assesses whether SE estimators target different things when there is
# massive variation at the block and cluster level. Equal-sized blocks
# and clusters reveal how FE vs. FI methods diverge under high variance.
#
# Key finding: Large sd.block with no sd.cluster produces FE SEs ~4x too
# large; FI SEs remain calibrated to the fixed-district model.
# Output: kable table + scatter plot comparing estimated vs. true SEs.
# Status: Runs with N_REPS = 5 as of 2026-06-04

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_block.R" ) )

N_REPS <- 5  # Set low for testing; increase for a full run


if ( FALSE ) {
    # Possible junk code according to Claude's audit
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

    grid <- build_blocks( K = 10, J.bar = 6, n.bar = 20, p.bar = 0.5,
                          ATE = 0.5,
                          sd.block = sdd,              # sd.district -> sd.block
                          sd.cluster = sds,            # sd.school -> sd.cluster
                          variable.p = FALSE,
                          cv.block = cv,               # cv.district -> cv.block
                          cv.nbar = cv,
                          cv.cluster = cv,             # cv.school -> cv.cluster
                          min.size = 4,
                          ICC.block = ICCd,            # ICC.district -> ICC.block
                          ICC.cluster = ICCs,          # ICC.school -> ICC.cluster
                          size.intercept.cor = 0,
                          size.impact.cor.cluster = 0,
                          R2.cluster = NULL,
                          num.cov.block = 0,
                          num.cov.cluster = 0,
                          num.cov.indiv = 0,
                          empirical = FALSE )

    rps = map_df( 1:3, \(x) {
        rr = single.CRT.fixed.block( grid = grid,      # fixed.district -> fixed.block
                                     n.runs = 3,
                                     ICC.block = ICCd ) # ICC.district -> ICC.block
        rr
    }, .id = "midrun" )

    rps
}

if ( FALSE ) {
    # testing
    one_run()
}


rps = map_dfr( 1:N_REPS, ~ one_run( FALSE ),
               .id = "run", .progress = TRUE )

# Analyze the simulation results ----

rps
nrow(rps) / length( unique(rps$method ) )

source( here::here( "R/aggregate_simulation_results.R" ) )

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
    knitr::kable(digits=2) %>%
    print()


se_low  = mean( res$SE.s ) * 0.9
se_high = mean( res$SE.s ) * 1.1

ggplot( res, aes( y=ESE_hat, x=SE.s, pch=weight ) ) +
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
