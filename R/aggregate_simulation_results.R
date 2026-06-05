


# Code to aggregate simulation results to calculate superpop and
# finite sample standard errors and so forth.


library( tidyverse )


source( here::here( "R/MCSE_functions.R" ) )



calculate_range_statistics <- function( simdata ) {
  simdata <- mutate( simdata,
                     weight = str_replace( weight, "Cluster-Cluster", "Cluster" ),
                     weight = str_replace( weight, "Person-Person", "Person" ) )
  
  
  range.dat <- simdata %>% 
    group_by( run, subrun ) %>% 
    mutate(range=diff(range(ATE.hat))) %>%
    ungroup() %>%
    group_by(subrun, weight, range) %>%  
    summarise(range_bygroup=diff(range(ATE.hat)) ) %>%
    pivot_wider( names_from=weight, values_from=range_bygroup, names_prefix="range_" )
  
  range.dat
}




#' This takes a single scenario and calculates performance metrics
#' given a nested simulation
calculate_performance_statistics <- function( simdata, 
                                              ATE.super.student = NULL,
                                              ATE.super.cluster = NULL,
                                              align = TRUE, 
                                              runID = "run" ) {
  
  simdata$run = simdata[[runID]]
  #simdata$subrun = simdata[[subrunID]]
  
  stopifnot( all( !is.na( simdata$method ) ) )
  
  
  ##### Calculate coverage ######
  
  if ( is.null( ATE.super.student ) ) {
    ATE.super.student = mean( simdata$ATE.finite.student )
    ATE.super.cluster = mean( simdata$ATE.finite.cluster )
  }
  
  cover.vector = function( ATE.hat, SE.hat, estimand ) {
    (ATE.hat - 1.96*SE.hat <= estimand) & (estimand <= ATE.hat + 1.96*SE.hat)
  }
  
  
  simdata = mutate( simdata, 
                    cov.f.stu = cover.vector( ATE.hat, SE.hat, ATE.finite.student ),
                    cov.f.clu = cover.vector( ATE.hat, SE.hat, ATE.finite.cluster),
                    #cov.f.dstu = cover.vector( ATE.hat, SE.hat, ATE.finite.dist.student ),
                    #cov.f..clu = cover.vector( ATE.hat, SE.hat, ATE.finite.dist.cluster),
                    
                    cov.s.stu = cover.vector( ATE.hat, SE.hat, ATE.super.student ),
                    cov.s.clu = cover.vector( ATE.hat, SE.hat, ATE.super.cluster ),
                    #cov.s.dstu = cover.vector( ATE.hat, SE.hat, ATE.super.dist.student ),
                    #cov.s..clu = cover.vector( ATE.hat, SE.hat, ATE.super.dist.cluster ))
  )
  
  
  
  
  ##### Aggegate by scenario and calculate summary stats for these scenarios.  #####
  
  sruns = simdata %>% 
    group_by( method, weight ) %>%   # per method... (aggregate over batch, run, subrun)
    summarise( mean.SE.hat = mean( SE.hat, na.rm = TRUE ),
               p.noSE = mean( is.na( SE.hat ) ),
               mean.ATE.hat=mean(ATE.hat),
               
               # bias of estimators for each estimand
               bias.f.stu = mean( ATE.hat - ATE.finite.student ),
               bias.f.clu = mean( ATE.hat - ATE.finite.cluster ),
               #bias.f.dstu = mean( ATE.hat - ATE.finite.dist.student ),
               #bias.f..clu = mean( ATE.hat - ATE.finite.dist.cluster ),
               bias.s.stu = mean( ATE.hat - ATE.super.student ),
               bias.s.clu = mean( ATE.hat - ATE.super.cluster ),
               #bias.s.dstu = mean( ATE.hat - ATE.super.dist.student ),
               #bias.s..clu = mean( ATE.hat - ATE.super.dist.cluster ),
               MCSE.bias = pmax( MCSE_bias( ATE.hat - ATE.finite.student ),
                                 MCSE_bias( ATE.hat - ATE.finite.cluster ) ),
               
               # RMSE of estimators for each estimand
               RMSE.f.stu = sqrt( mean( ( ATE.hat - ATE.finite.student )^2 ) ),
               RMSE.f.clu = sqrt( mean( ( ATE.hat - ATE.finite.cluster   )^2 ) ),
               #RMSE.f.dstu = sqrt( mean( ( ATE.hat - ATE.finite.dist.student )^2 ) ),
               #RMSE.f..clu = sqrt( mean( ( ATE.hat - ATE.finite.dist.cluster   )^2 ) ),
               RMSE.s.stu = sqrt( mean( ( ATE.hat - ATE.super.student  )^2 ) ),
               RMSE.s.clu = sqrt( mean( ( ATE.hat - ATE.super.cluster )^2 ) ),
               #RMSE.s.dstu = sqrt( mean( ( ATE.hat - ATE.super.dist.student  )^2 ) ),
               #RMSE.s..clu = sqrt( mean( ( ATE.hat - ATE.super.dist.cluster )^2 ) ),
               MCSE.RMSE = pmax( MCSE_SE( ATE.hat - ATE.finite.student ),
                                 MCSE_SE( ATE.hat - ATE.finite.cluster ) ),
               
               # true superpopulation sd of estimate
               SE.s = sd( ATE.hat ),
               MCSE.SE = MCSE_SE( ATE.hat ),
               
               # variation in estimated SEs
               sd.SE.hat.s = sd( SE.hat, na.rm=TRUE ),
               
               # coverage rates
               cover.f.stu = mean( cov.f.stu, na.rm=TRUE ),
               cover.f.clu = mean( cov.f.clu, na.rm=TRUE ),
               #cover.f.dstu = mean( cov.f.dstu ),
               #cover.f..clu = mean( cov.f..clu ),
               cover.s.stu = mean( cov.s.stu, na.rm=TRUE ),
               cover.s.clu = mean( cov.s.clu, na.rm=TRUE ),
               #cover.s.dstu = mean( cov.s.dstu ),
               #cover.s..clu = mean( cov.s..clu ),
               MCSE.cover = pmax( MCSE_bias( cov.f.stu ),
                                  MCSE_bias( cov.f.clu ) ),
               
               avg.df = mean( df, na.rm=TRUE ),
               
               power = mean( p_value <= 0.05 ),
               
               # # SD of finite-sample estimands (?)
               # sd.ATE.f.stu = sd( ATE.finite.student ),
               # sd.ATE.f.clu = sd( ATE.finite.cluster ),
               # sd.ATE.f.dstu = sd( ATE.finite.dist.student ),
               # sd.ATE.f..clu = sd( ATE.finite.dist.cluster ),
               
               # var.err2.s.stu = var( ATE.hat - ATE.super.student ),
               # var.err2.s.clu = var( ATE.hat - ATE.super.cluster ),
               # var.err2.s.dstu = var( ATE.hat - ATE.super.dist.student ),
               # var.err2.s..clu = var( ATE.hat - ATE.super.dist.cluster )
               
               .groups = "drop"
    )
  
  
  
  ##### calculate the finite sample SE using our nested simulations  #####
  
  subruns = simdata %>% 
    group_by( method, run ) %>%            # just aggregating over subruns
    summarise( nsr = n(),
               
               # noisy estimate of true finite-sample SE^2:
               SE2 = var( ATE.hat ),           
               # SE.hat.f2 = mean( SE.hat^2 ),
               
               # variation in estimated SE:
               var.SE.hat = var( SE.hat, na.rm=TRUE ), 
               
               .groups="drop" )    
  
  subruns <- subruns %>% 
    group_by( method ) %>%                            # now aggregate over runs
    summarise( 
      # true finite-sample SE^2:
      SE.f = sqrt( mean( SE2 ) ),                
      
      #SE.hat.f = sqrt( mean( SE.hat.f2 ) ),
      
      # average sd of SE estimates:
      sd.SE.hat.f = sqrt( mean( var.SE.hat, na.rm=TRUE ) ) 
    )
  
  # add to sruns
  sruns = merge( sruns, subruns, by=c("method"), all=TRUE )
  
  
  # Return the final results ----
  sruns <- sruns %>%
    as_tibble()
  
  if ( align ) {
    sruns <- align_performance_stats( sruns )
  }
  
  sruns$mean.SE.hat[ is.nan( sruns$mean.SE.hat ) ] <- NA
  sruns$cover[ is.nan( sruns$cover ) ] <- NA
  
  sruns
}



# Align performances with the proper estimand to simplify performance stats ----

align_performance_stats <- function( stats ) {
  
  stats %>%
    mutate( 
      bias = ifelse( grepl( "Student|Person", weight ),
                     bias.f.stu, bias.f.clu ),
      RMSE = ifelse( grepl( "Student|Person", weight ), 
                     RMSE.f.stu, RMSE.f.clu ),
      cover = ifelse( grepl( "Student|Person", weight ),
                      cover.f.stu, cover.f.clu )
    ) %>%
    dplyr::select( method, weight, 
                   mean.SE.hat, mean.ATE.hat, 
                   bias, MCSE.bias, 
                   RMSE, MCSE.RMSE, 
                   cover, MCSE.cover,
                   p.noSE, avg.df,
                   SE.s, SE.f, MCSE.SE )
  
}



# calc_performance: lightweight summary for small-simulation scripts ----
# (Moved here from small_simulations/simulation_DGP.R)

calc_performance <- function( rps ) {
    res <- rps %>%
        group_by( method ) %>%
        summarise( EATE    = mean( ATE_hat ),
                   tau_p   = mean( tau_p ),
                   tau_c   = mean( tau_c ),
                   SE      = sd( ATE_hat ),
                   SE_E    = sd( ATE_hat ) / sqrt(n()),
                   ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
        mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )
    res
}


# Test code ----

if ( FALSE ) {
  
  length( unique( simdata$method ) )
  
  calculate_range_statistics( simdata )
  
  calculate_performance_statistics( simdata )
}


