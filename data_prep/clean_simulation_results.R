
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
head( scenarios$run[[10]] )

table( scenarios$J, scenarios$tau )
table( scenarios$proptx.dependence, scenarios$dependence )
table( scenarios$variable.n, scenarios$variable.p )

# Make the names of the factors better
scenarios$dependence = ifelse( scenarios$dependence == 1, "corr", "indep" )
scenarios$proptx.dependence = ifelse( scenarios$proptx.dependence == 1, "pcorr", 
                                      ifelse( scenarios$proptx.dependence == -1, "ncorr", "indep" ) )


# Cleverness: Merge the novariation and correlation factors to make single factor!
table( scenarios$dependence, scenarios$variable.n )
scenarios$dependence = ifelse( scenarios$variable.n, scenarios$dependence, "const" )
scenarios$dependence = factor( scenarios$dependence, c( "const", "indep", "corr" ) )
table( scenarios$dependence )

table( scenarios$proptx.dependence, scenarios$variable.p )
scenarios$proptx.dependence = ifelse( scenarios$variable.p, scenarios$proptx.dependence, "const" )
scenarios$proptx.dependence = factor( scenarios$proptx.dependence, c( "ncorr", "const", "indep", "pcorr" ) )
table( scenarios$proptx.dependence )

# Turn to cross site variation in SD
scenarios = mutate( scenarios, tau = sqrt( tau ) )


# simple checks on data structure
head( scenarios$run[[40]] )
table( table( scenarios$run[[40]]$subrun ) )
table( table( scenarios$run[[40]]$run ) )
filter( scenarios$run[[40]], subrun=="1.3.1" )


# Make a dataframe of each simulation
simdata = unnest( scenarios )

# Drop simulations from the scenarios frame
scenarios$run = NULL

# Drop uninteresting estimators
simdata = filter( simdata, !( method %in% c( "plug_in_big", "hybrid_p", "hybrid_m", "FE-IPTW(n)" )))

simdata


##### drop failed runs due to lmer crashing, etc.  #####

filter( simdata, is.na( method ) )
simdata = filter( simdata, !is.na( method ) )

failed = filter( simdata, is.na( SE.hat ) )
nrow( failed )
table( failed$method, useNA ="always" )

if ( nrow( failed ) > 0 ) {
  warning( "Dropping some runs without estimated SEs" )
  simdata = filter( simdata, !is.na( SE.hat ) )
}



##### Calculate coverage ######

cover.vector = function( ATE.hat, SE.hat, estimand ) {
  (ATE.hat - 1.96*SE.hat <= estimand) & (estimand <= ATE.hat + 1.96*SE.hat)
}

simdata = mutate( simdata, 
                  cov.fp = cover.vector( ATE.hat, SE.hat, ATE.finite.person ),
                  cov.fs = cover.vector( ATE.hat, SE.hat, ATE.finite.site),
                  cov.sp = cover.vector( ATE.hat, SE.hat, ATE.super.person ),
                  cov.ss = cover.vector( ATE.hat, SE.hat, ATE.super.site ) )



##### Aggegate by scenario and calculate summary stats for these scenarios.  #####
sruns = simdata %>% group_by( ID, method ) %>%
  summarise( mean.SE.hat = mean( SE.hat ),
             bias.fp = mean( ATE.hat - ATE.finite.person ),
             bias.fs = mean( ATE.hat - ATE.finite.site ),
             bias.sp = mean( ATE.hat - ATE.super.person ),
             bias.ss = mean( ATE.hat - ATE.super.site ),
             RMSE.fp = sqrt( mean( ( ATE.hat - ATE.finite.person )^2 ) ),
             RMSE.fs = sqrt( mean( ( ATE.hat - ATE.finite.site   )^2 ) ),
             RMSE.sp = sqrt( mean( ( ATE.hat - ATE.super.person  )^2 ) ),
             RMSE.ss = sqrt( mean( ( ATE.hat - ATE.super.site    )^2 ) ),
             sd.ATE.fp = sd( ATE.finite.person ),
             sd.ATE.fs = sd( ATE.finite.site ),
             SE.s = sd( ATE.hat ),
             sd.SE.hat.s = sd( SE.hat ),
             cover.fp = mean( cov.fp ),
             cover.fs = mean( cov.fs ),
             cover.sp = mean( cov.sp ),
             cover.ss = mean( cov.ss ),
             var.err2.sp = var( ATE.hat - ATE.super.person ),
             var.err2.ss = var( ATE.hat - ATE.super.site ) )

# Add in the info on the scenario to the results.
sruns = merge( scenarios, sruns, by="ID" )
head( sruns )


##### Ensure no bad sim runs #####
# look for missing SE estimates
badIDs = filter( sruns, is.na( mean.SE.hat ) | is.na( RMSE.fp )  )
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

# Get characteristics of the different methods
mc = blkvar::method.characteristics()
mc

# Add method stats to the simulation results data
sruns = merge( sruns, mc, by="method" )

rm( mc )



##### Save #####

save( scenarios, simdata, sruns, file="results/processed_simulation_results.RData" )

