
# Load the simulation results and get them ready for processing.

# End up with three dataframes:
# scenarios: raw simulation data
# simdata: unnested simulation data
# sruns: aggregate statistics for each method for each simulation run
#
# These will be saved in 'processed_simulation_results.RData'

library( tidyverse )

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd( ".." )

if ( FALSE ) {
  load( "SimulationScenariosCompleteRnd3_Tst.RData" )
  sc2 = scenarios
  load( "SimulationScenariosCompleteRnd3_PtB.RData" )
  nrow( scenarios )
  nrow( sc2 )
  sc2$ID = sc2$ID + 1000
  full_scenarios = bind_rows( scenarios, sc2 )
  rm( sc2, scenarios )
  full_scenarios$variable.n = full_scenarios$variable.p = TRUE
  
  load( "SimulationScenariosCompleteRnd3PartE.RData" )
  nrow( scenarios )
  scenarios$ID = scenarios$ID + 2000
  full_scenarios = bind_rows( full_scenarios, scenarios )
  rm( scenarios )
  scenarios = full_scenarios
}

scenarios = readRDS( "results/complete_odyssey_results.rds" )



head( scenarios )
nrow( scenarios )
scenarios$ID
head( scenarios$run[[1]] )







# # Make the names of the factors better
# scenarios$dependence = ifelse( scenarios$dependence == 1, "corr", "indep" )
# scenarios$proptx.dependence = ifelse( scenarios$proptx.dependence == 1, "pcorr", 
#                                       ifelse( scenarios$proptx.dependence == -1, "ncorr", "indep" ) )
# 
# 
# # Cleverness: Merge the novariation and correlation factors to make single factor!
# table( scenarios$dependence, scenarios$variable.n )
# scenarios$dependence = ifelse( scenarios$variable.n, scenarios$dependence, "const" )
# scenarios$dependence = factor( scenarios$dependence, c( "const", "indep", "corr" ) )
# table( scenarios$dependence )
# 
# table( scenarios$proptx.dependence, scenarios$variable.p )
# scenarios$proptx.dependence = ifelse( scenarios$variable.p, scenarios$proptx.dependence, "const" )
# scenarios$proptx.dependence = factor( scenarios$proptx.dependence, c( "ncorr", "const", "indep", "pcorr" ) )
# table( scenarios$proptx.dependence )
# 
# # Turn to cross site variation in SD
# scenarios = mutate( scenarios, tau = sqrt( tau ) )
# 
# 
# # simple checks on data structure
# head( scenarios$run[[40]] )
# table( table( scenarios$run[[40]]$subrun ) )
# table( table( scenarios$run[[40]]$run ) )
# filter( scenarios$run[[40]], subrun=="1.3.1" )


# Make a dataframe of each simulation
simdata = unnest( scenarios, cols=c(run) )

# Drop simulations from the scenarios frame
scenarios$run = NULL

# # Drop uninteresting estimators
# simdata = filter( simdata, !( method %in% c( "plug_in_big", "hybrid_p", "hybrid_m", "FE-IPTW(n)" )))

simdata


##### drop failed runs due to lmer crashing, etc.  #####

filter( simdata, is.na( method ) )
simdata = filter( simdata, !is.na( method ) )

failed = filter( simdata, is.na( SE.hat ) | is.infinite( SE.hat ) )
nrow( failed )
table( failed$method, useNA ="always" )

if ( nrow( failed ) > 0 ) {
  warning( "Dropping some runs without estimated SEs" )
  simdata = filter( simdata, !(is.na( SE.hat ) | is.infinite( SE.hat )) )
}



##### Calculate coverage ######

cover.vector = function( ATE.hat, SE.hat, estimand ) {
  (ATE.hat - 1.96*SE.hat <= estimand) & (estimand <= ATE.hat + 1.96*SE.hat)
}

simdata = mutate( simdata, 
                  cov.f.stu = cover.vector( ATE.hat, SE.hat, ATE.finite.student ),
                  cov.f.sch = cover.vector( ATE.hat, SE.hat, ATE.finite.school),
                  cov.f.dstu = cover.vector( ATE.hat, SE.hat, ATE.finite.dist.student ),
                  cov.f.dsch = cover.vector( ATE.hat, SE.hat, ATE.finite.dist.school),
                  
                  cov.s.stu = cover.vector( ATE.hat, SE.hat, ATE.super.student ),
                  cov.s.sch = cover.vector( ATE.hat, SE.hat, ATE.super.school ),
                  cov.s.dstu = cover.vector( ATE.hat, SE.hat, ATE.super.dist.student ),
                  cov.s.dsch = cover.vector( ATE.hat, SE.hat, ATE.super.dist.school ))


##### Aggegate by scenario and calculate summary stats for these scenarios.  #####
sruns = simdata %>% group_by( ID, method ) %>%
  summarise( mean.SE.hat = mean( SE.hat ),
             
             bias.f.stu = mean( ATE.hat - ATE.finite.student ),
             bias.f.sch = mean( ATE.hat - ATE.finite.school ),
             bias.f.dstu = mean( ATE.hat - ATE.finite.dist.student ),
             bias.f.dsch = mean( ATE.hat - ATE.finite.dist.school ),
             bias.s.stu = mean( ATE.hat - ATE.super.student ),
             bias.s.sch = mean( ATE.hat - ATE.super.school ),
             bias.s.dstu = mean( ATE.hat - ATE.super.dist.student ),
             bias.s.dsch = mean( ATE.hat - ATE.super.dist.school ),
             
             RMSE.f.stu = sqrt( mean( ( ATE.hat - ATE.finite.student )^2 ) ),
             RMSE.f.sch = sqrt( mean( ( ATE.hat - ATE.finite.school   )^2 ) ),
             RMSE.f.dstu = sqrt( mean( ( ATE.hat - ATE.finite.dist.student )^2 ) ),
             RMSE.f.dsch = sqrt( mean( ( ATE.hat - ATE.finite.dist.school   )^2 ) ),
             RMSE.s.stu = sqrt( mean( ( ATE.hat - ATE.super.student  )^2 ) ),
             RMSE.s.sch = sqrt( mean( ( ATE.hat - ATE.super.school )^2 ) ),
             RMSE.s.dstu = sqrt( mean( ( ATE.hat - ATE.super.dist.student  )^2 ) ),
             RMSE.s.dsch = sqrt( mean( ( ATE.hat - ATE.super.dist.school )^2 ) ),
             
             sd.ATE.f.stu = sd( ATE.finite.student ),
             sd.ATE.f.sch = sd( ATE.finite.school ),
             sd.ATE.f.dstu = sd( ATE.finite.dist.student ),
             sd.ATE.f.dsch = sd( ATE.finite.dist.school ),
             
             SE.s = sd( ATE.hat ),
             sd.SE.hat.s = sd( SE.hat ),
             
             cover.f.stu = mean( cov.f.stu ),
             cover.f.sch = mean( cov.f.sch ),
             cover.f.dstu = mean( cov.f.dstu ),
             cover.f.dsch = mean( cov.f.dsch ),
             cover.s.stu = mean( cov.s.stu ),
             cover.s.sch = mean( cov.s.sch ),
             cover.s.dstu = mean( cov.s.dstu ),
             cover.s.dsch = mean( cov.s.dsch ),
             
             var.err2.s.stu = var( ATE.hat - ATE.super.student ),
             var.err2.s.sch = var( ATE.hat - ATE.super.school ),
             var.err2.s.dstu = var( ATE.hat - ATE.super.dist.student ),
             var.err2.s.dsch = var( ATE.hat - ATE.super.dist.school ))

# Add in the info on the scenario to the results.
sruns = merge( scenarios, sruns, by="ID" )
head( sruns )


##### Ensure no bad sim runs #####
# look for missing SE estimates
badIDs = filter( sruns, is.na( mean.SE.hat ) | is.na( RMSE.f.stu )  )
if ( nrow( badIDs ) > 0 ) {
  warning( "Continuing to drop runs with missing estimated SEs or RMSEs" )
  badIDs = unique( badIDs$ID )
  cat( "Missing SE estimates:\n" )
  print( badIDs )
  sruns = filter( sruns, !is.na( mean.SE.hat ) & !is.na( RMSE.fp ) )  
  stop()
}

# missing IDs from prior filtering
if ( nrow( failed ) > 0 ) {
  
  cat( "Failed runs--looking at where they showed up" )
  badIDs = unique( failed$ID )
  print( badIDs )
  filter( scenarios, ID %in% badIDs )
  stop()
}

# Which IDs are missing stats for some of the methods?
ss = sruns %>% ungroup() %>% dplyr::select( ID, method, mean.SE.hat ) %>%
  spread( method, mean.SE.hat )
head( ss )
not.comp = ss[ !complete.cases(ss), ]
if ( nrow( not.comp ) > 0 ) {
  warning( "Not complete cases found\n" )
  print( not.comp )
  stop()
}


##### calculate the finite sample SE using our nested simulations  #####

subruns = simdata %>% group_by( ID, method, run ) %>% 
    summarise( nsr = n(),
             SE2 = var( ATE.hat ),
             SE.hat.f2 = mean( SE.hat^2 ),
             var.SE.hat = var( SE.hat ) )

head( subruns )
table( subruns$nsr )
subruns = filter( subruns, nsr == 3 )
subruns = subruns %>% group_by( ID, method ) %>%
  summarise( SE.f = sqrt( mean( SE2 ) ),
             sd.SE.hat.f = sqrt( mean( var.SE.hat ) ) )
             #SE.hat.f = sqrt( mean( SE.hat.f2 ) ) )
nrow( sruns )
nrow( subruns )
sruns = merge( sruns, subruns, by=c("ID","method"), all=TRUE )
nrow( sruns )
rm( subruns )


head( sruns )


# Check: How do the true SEs of superpop compare to the finite pop
# (These should generally be above 1.0)
summary( sruns$SE.s / sruns$SE.f )

odds = sruns %>% mutate( ratio = SE.s / SE.f ) %>%
  filter( ratio < 1.0 )
table( odds$ID )
filter( scenarios, ID %in% odds$ID )


# inflate captures how well the estimated standard error predicts the true
# (averaged across finite samples for finite) SE.
sruns = mutate( sruns, 
                inflate.f = mean.SE.hat / SE.f,
                inflate.s = mean.SE.hat / SE.s )
head( sruns )

nrow( sruns )



##### Add in the different method characteristics #####

# # Get characteristics of the different methods
# mc = clusterRCT::method.characteristics()
# mc
# 
# # Add method stats to the simulation results data
# sruns = merge( sruns, mc, by="method" )
# 
# rm( mc )



##### Save #####

save( scenarios, simdata, sruns, file="results/processed_simulation_results.RData" )

