


# Code to aggregate simulation results to calculate superpop and
# finite sample standard errors and so forth.


library( tidyverse )


source( here::here( "R/MCSE_functions.R" ) )



calculate_range_statistics <- function( simdata ) {
  simdata <- mutate( simdata,
                     weight = str_replace( weight, "Cluster-Cluster", "Cluster" ),
                     weight = str_replace( weight, "Person-Person", "Person" ) )
  
  
  range.dat <- simdata %>% 
    group_by( run,subrun ) %>% 
    mutate(range=diff(range(ATE.hat))) %>%
    ungroup() %>%
    group_by(subrun,weight, range) %>%  
    summarise(range_bygroup=diff(range(ATE.hat)) ) %>%
    pivot_wider( names_from=weight, values_from=range_bygroup, names_prefix="range_" )
  
  range.dat
}




#' This takes a single scenario and calculates performance metrics
#' given a nested simulation
calculate_performance_statistics <- function( simdata, 
                                              ATE.super.student = NULL,
                                              ATE.super.school = NULL,
                                              align = TRUE ) {
  
  stopifnot( all( !is.na( simdata$method ) ) )
  
  
  ##### Calculate coverage ######
  
  if ( is.null( ATE.super.student ) ) {
    ATE.super.student = mean( simdata$ATE.finite.student )
    ATE.super.school = mean( simdata$ATE.finite.school )
  }
  
  cover.vector = function( ATE.hat, SE.hat, estimand ) {
    (ATE.hat - 1.96*SE.hat <= estimand) & (estimand <= ATE.hat + 1.96*SE.hat)
  }
  
  
  simdata = mutate( simdata, 
                    cov.f.stu = cover.vector( ATE.hat, SE.hat, ATE.finite.student ),
                    cov.f.sch = cover.vector( ATE.hat, SE.hat, ATE.finite.school),
                    #cov.f.dstu = cover.vector( ATE.hat, SE.hat, ATE.finite.dist.student ),
                    #cov.f.dsch = cover.vector( ATE.hat, SE.hat, ATE.finite.dist.school),
                    
                    cov.s.stu = cover.vector( ATE.hat, SE.hat, ATE.super.student ),
                    cov.s.sch = cover.vector( ATE.hat, SE.hat, ATE.super.school ),
                    #cov.s.dstu = cover.vector( ATE.hat, SE.hat, ATE.super.dist.student ),
                    #cov.s.dsch = cover.vector( ATE.hat, SE.hat, ATE.super.dist.school ))
  )
  
  
  
  
  ##### Aggegate by scenario and calculate summary stats for these scenarios.  #####
  
  sruns = simdata %>% 
    group_by( method, weight ) %>%   # per method... (aggregate over batch, run, subrun)
    summarise( mean.SE.hat = mean( SE.hat, na.rm = TRUE ),
               p.noSE = mean( is.na( SE.hat ) ),
               mean.ATE.hat=mean(ATE.hat),
               # bias of estimators for each estimand
               bias.f.stu = mean( ATE.hat - ATE.finite.student ),
               bias.f.sch = mean( ATE.hat - ATE.finite.school ),
               #bias.f.dstu = mean( ATE.hat - ATE.finite.dist.student ),
               #bias.f.dsch = mean( ATE.hat - ATE.finite.dist.school ),
               bias.s.stu = mean( ATE.hat - ATE.super.student ),
               bias.s.sch = mean( ATE.hat - ATE.super.school ),
               #bias.s.dstu = mean( ATE.hat - ATE.super.dist.student ),
               #bias.s.dsch = mean( ATE.hat - ATE.super.dist.school ),
               MCSE.bias = pmax( MCSE_bias( ATE.hat - ATE.finite.student ),
                                 MCSE_bias( ATE.hat - ATE.finite.school ) ),
               
               # RMSE of estimators for each estimand
               RMSE.f.stu = sqrt( mean( ( ATE.hat - ATE.finite.student )^2 ) ),
               RMSE.f.sch = sqrt( mean( ( ATE.hat - ATE.finite.school   )^2 ) ),
               #RMSE.f.dstu = sqrt( mean( ( ATE.hat - ATE.finite.dist.student )^2 ) ),
               #RMSE.f.dsch = sqrt( mean( ( ATE.hat - ATE.finite.dist.school   )^2 ) ),
               RMSE.s.stu = sqrt( mean( ( ATE.hat - ATE.super.student  )^2 ) ),
               RMSE.s.sch = sqrt( mean( ( ATE.hat - ATE.super.school )^2 ) ),
               #RMSE.s.dstu = sqrt( mean( ( ATE.hat - ATE.super.dist.student  )^2 ) ),
               #RMSE.s.dsch = sqrt( mean( ( ATE.hat - ATE.super.dist.school )^2 ) ),
               MCSE.RMSE = pmax( MCSE_SE( ATE.hat - ATE.finite.student ),
                                 MCSE_SE( ATE.hat - ATE.finite.school ) ),
               # # SD of finite-sample estimands (?)
               # sd.ATE.f.stu = sd( ATE.finite.student ),
               # sd.ATE.f.sch = sd( ATE.finite.school ),
               # sd.ATE.f.dstu = sd( ATE.finite.dist.student ),
               # sd.ATE.f.dsch = sd( ATE.finite.dist.school ),
               
               # true superpopulation sd of estimate
               SE.s = sd( ATE.hat ),
               MCSE.SE = MCSE_SE( ATE.hat ),
               
               # variation in estimated SEs
               sd.SE.hat.s = sd( SE.hat, na.rm=TRUE ),
               
               # coverage rates
               cover.f.stu = mean( cov.f.stu, na.rm=TRUE ),
               cover.f.sch = mean( cov.f.sch, na.rm=TRUE ),
               #cover.f.dstu = mean( cov.f.dstu ),
               #cover.f.dsch = mean( cov.f.dsch ),
               cover.s.stu = mean( cov.s.stu, na.rm=TRUE ),
               cover.s.sch = mean( cov.s.sch, na.rm=TRUE ),
               #cover.s.dstu = mean( cov.s.dstu ),
               #cover.s.dsch = mean( cov.s.dsch ),
               MCSE.cover = pmax( MCSE_bias( cov.f.stu ),
                                  MCSE_bias( cov.f.sch ) ),
               avg.df = mean( df, na.rm=TRUE ),
               
               # var.err2.s.stu = var( ATE.hat - ATE.super.student ),
               # var.err2.s.sch = var( ATE.hat - ATE.super.school ),
               # var.err2.s.dstu = var( ATE.hat - ATE.super.dist.student ),
               # var.err2.s.dsch = var( ATE.hat - ATE.super.dist.school )
               .groups = "drop"
    )
  
  
  
  ##### calculate the finite sample SE using our nested simulations  #####
  
  subruns = simdata %>% 
    group_by( method, run ) %>%            # just aggregating over subruns
    summarise( nsr = n(),
               SE2 = var( ATE.hat ),           # noisy estimate of true finite-sample SE^2
               # SE.hat.f2 = mean( SE.hat^2 ),
               var.SE.hat = var( SE.hat, na.rm=TRUE ), .groups="drop" )    # variation in estimated SE
  
  subruns <- subruns %>% 
    group_by( method ) %>%                            # now aggregate over runs
    summarise( SE.f = sqrt( mean( SE2 ) ),                # true finite-sample SE^2
               #SE.hat.f = sqrt( mean( SE.hat.f2 ) ),
               sd.SE.hat.f = sqrt( mean( var.SE.hat, na.rm=TRUE ) ) ) # average sd of SE estimates
  
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
                     bias.f.stu, bias.f.sch ),
      RMSE = ifelse( grepl( "Student|Person", weight ), 
                     RMSE.f.stu, RMSE.f.sch ),
      cover = ifelse( grepl( "Student|Person", weight ),
                      cover.f.stu, cover.f.sch )
    ) %>%
    dplyr::select( method, weight, 
                   mean.SE.hat, mean.ATE.hat, 
                   bias, MCSE.bias, 
                   RMSE, MCSE.RMSE, 
                   cover, MCSE.cover,
                   p.noSE, avg.df,
                   SE.s, SE.f, MCSE.SE )
  
}



# Test code ----

if ( FALSE ) {
  
  length( unique( simdata$method ) )
  
  calculate_range_statistics( simdata )
  
  calculate_performance_statistics( simdata )
}


