
# Load the simulation results and get them ready for processing.

# End up with three dataframes:
# scenarios: raw simulation data, one row per simulation scenario
# simdata: unnested simulation data
# sruns: aggregate statistics for each method for each simulation run
#
# These will be saved in 'processed_simulation_results.RData'


rm(list=ls())
library( tidyverse )

scenarios = readRDS( "results/complete_odyssey_results.rds" )

simdata <- scenarios %>% 
  unnest(run)
table(simdata$method)
nrow(simdata)

#rm(list=ls())
#library( tidyverse )
#load("results_sav/processed_simulation_results_estimators_estimand_2024-09-30.RData")
#load("results_sav/processed_simulation_results_sim_base_1_2024-09-30.RData")
#load("results_sav/processed_simulation_results_estimators_estimand_2_2024-09-30.RData")

# Drop simulations from the scenarios frame
scenarios$run = NULL
nrow(simdata)

# # Drop uninteresting estimators
# simdata = filter( simdata, !( method %in% c( "plug_in_big", "hybrid_p", "hybrid_m", "FE-IPTW(n)" )))


##### drop failed runs due to lmer crashing, etc.  #####

# TODO: Q: how is method ever NA?
failed <- simdata %>% 
  filter(is.na(method) | is.na(SE.hat) | is.infinite(SE.hat))

#don't drop the runs without the estimated SEs 
# table( failed$method, useNA ="always" )
#if ( nrow( failed ) > 0 ) {
 # warning( "Dropping some runs without estimated SEs" )
  #simdata = filter( simdata, !(is.na( SE.hat ) | is.infinite( SE.hat )) )
#}

nrow(simdata)
with(simdata,table(is.na(method)))
with(simdata,table(is.na(ATE.hat)))
with(simdata,table(is.na(SE.hat)))

##### Calculate coverage ######

cover.vector = function( ATE.hat, SE.hat, df, estimand ) {
  t=qt(0.975,df)
  (ATE.hat - t*SE.hat <= estimand) & (estimand <= ATE.hat + t*SE.hat)
}

CI.length = function( ATE.hat, SE.hat, df, estimand ) {
  t=qt(0.975,df)
  2*t*SE.hat
}


simdata = mutate(simdata, 
                 # Coverage for finite populations
                 cov.f.stu = cover.vector(ATE.hat, SE.hat, df, ATE.finite.student),
                 cov.f.sch = cover.vector(ATE.hat, SE.hat, df, ATE.finite.school),
                 cov.f.dstu = cover.vector(ATE.hat, SE.hat, df, ATE.finite.dist.student),
                 cov.f.dsch = cover.vector(ATE.hat, SE.hat, df, ATE.finite.dist.school),
                 
                 # Coverage for super populations
                 cov.s.stu = cover.vector(ATE.hat, SE.hat, df, ATE.super.student),
                 cov.s.sch = cover.vector(ATE.hat, SE.hat, df, ATE.super.school),
                 cov.s.dstu = cover.vector(ATE.hat, SE.hat, df, ATE.super.dist.student),
                 cov.s.dsch = cover.vector(ATE.hat, SE.hat, df, ATE.super.dist.school),
                 
                 # CI length
                 CI.length = CI.length(ATE.hat, SE.hat, df))

nrow(simdata)
table(simdata$method)

# Finite population coverage
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.f.stu)))
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.f.sch)))
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.f.dstu)))
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.f.dsch)))

# Super population coverage
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.s.stu)))
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.s.sch)))
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.s.dstu)))
with(simdata,table((is.na(SE.hat)|df==0),is.na(cov.s.dsch)))

# CI length
with(simdata,table((is.na(SE.hat)|df==0),is.na(CI.length)))

failed.se <- simdata %>% filter(!(is.na(SE.hat)) & is.na(cov.f.stu))
table(failed.se$df)
with(simdata,table(round(J.bar,0),round(df,0)))

#add in ranges across methods 
#range.dat <-simdata %>% 
 #                 filter(!(weight %in% c("cluster/block","person/block")) & !(method %in% c("DB_HT","DB_Raj"))) %>% 
  #                group_by(ID,subrun ) %>% 
   #               mutate(range=diff(range(ATE.hat))) %>%
    #              ungroup() %>%
     #            group_by(ID,subrun,weight,population) %>%  
#                summarise(range_bygroup=diff(range(ATE.hat)),range=first(range))

#range.dat = merge( scenarios, range.dat, by="ID" )

##### Aggegate by scenario and calculate summary stats for these scenarios.  #####

#confirm weight, bias and disfavored are unique by method



# Check if variables are unique by ID and method

duplicates <- simdata %>%
  group_by(method) %>%
  summarize(
    unique_biased_values = n_distinct(biased),
    unique_disfavored_values = n_distinct(disfavored),
    unique_weight_values = n_distinct(weight),
    .groups = "drop"
  ) %>%
  filter(unique_biased_values > 1 | unique_disfavored_values > 1 | unique_weight_values > 1)


# If empty, variables are unique by ID and method
if(nrow(duplicates) == 0) {
  print("Variables are unique by ID and method")
} else {
  print(paste("Found", nrow(duplicates), "non-unique combinations"))
  print(duplicates)
}

nrow(simdata)
table(simdata$method)

sruns = simdata %>% 
  group_by( ID, method, weight, disfavored, biased) %>%   # per scenario and method... (aggregate over batch, run, subrun)
  summarise( mean.SE.hat = mean( SE.hat ),
             mean.ATE.hat=mean(ATE.hat),
             mean.df=mean(df),
             # bias of estimators for each estimand
             bias.f.stu = mean( ATE.hat - ATE.finite.student ),
             bias.f.sch = mean( ATE.hat - ATE.finite.school ),
             bias.f.dstu = mean( ATE.hat - ATE.finite.dist.student ),
             bias.f.dsch = mean( ATE.hat - ATE.finite.dist.school ),
             bias.s.stu = mean( ATE.hat - ATE.super.student ),
             bias.s.sch = mean( ATE.hat - ATE.super.school ),
             bias.s.dstu = mean( ATE.hat - ATE.super.dist.student ),
             bias.s.dsch = mean( ATE.hat - ATE.super.dist.school ),
             
             # RMSE of estimators for each estimand
             RMSE.f.stu = sqrt( mean( ( ATE.hat - ATE.finite.student )^2 ) ),
             RMSE.f.sch = sqrt( mean( ( ATE.hat - ATE.finite.school   )^2 ) ),
             RMSE.f.dstu = sqrt( mean( ( ATE.hat - ATE.finite.dist.student )^2 ) ),
             RMSE.f.dsch = sqrt( mean( ( ATE.hat - ATE.finite.dist.school   )^2 ) ),
             RMSE.s.stu = sqrt( mean( ( ATE.hat - ATE.super.student  )^2 ) ),
             RMSE.s.sch = sqrt( mean( ( ATE.hat - ATE.super.school )^2 ) ),
             RMSE.s.dstu = sqrt( mean( ( ATE.hat - ATE.super.dist.student  )^2 ) ),
             RMSE.s.dsch = sqrt( mean( ( ATE.hat - ATE.super.dist.school )^2 ) ),
             
             # # SD of finite-sample estimands
              #sd.ATE.f.stu = sd( ATE.finite.student ),
              #sd.ATE.f.sch = sd( ATE.finite.school ),
              #sd.ATE.f.dstu = sd( ATE.finite.dist.student ),
              #sd.ATE.f.dsch = sd( ATE.finite.dist.school ),
             
             # true superpopulation sd of estimate
             SE.s = sd( ATE.hat ),
             
             # variation in estimated SEs
             sd.SE.hat.s = sd( SE.hat ),
             
             # coverage rates
             cover.f.stu = mean( cov.f.stu ),
             cover.f.sch = mean( cov.f.sch ),
             cover.f.dstu = mean( cov.f.dstu ),
             cover.f.dsch = mean( cov.f.dsch ),
             cover.s.stu = mean( cov.s.stu ),
             cover.s.sch = mean( cov.s.sch ),
             cover.s.dstu = mean( cov.s.dstu ),
             cover.s.dsch = mean( cov.s.dsch ),
             CI.length=mean(CI.length))
             
             
             # var.err2.s.stu = var( ATE.hat - ATE.super.student ),
             # var.err2.s.sch = var( ATE.hat - ATE.super.school ),
             # var.err2.s.dstu = var( ATE.hat - ATE.super.dist.student ),
             # var.err2.s.dsch = var( ATE.hat - ATE.super.dist.school )

# Add in the info on the scenario to the results.
nrow(sruns)
table(sruns$method)

sruns = merge( scenarios, sruns, by="ID" )
head( sruns )
nrow(sruns)
table(sruns$method)

##### Ensure no bad sim runs #####

# look for missing SE estimates
badIDs = filter( sruns, is.na( mean.SE.hat ) | is.na( RMSE.f.stu )  )
if ( nrow( badIDs ) > 0 ) {
  #warning( "Continuing to drop runs with missing estimated SEs or RMSEs" )
  badIDs = unique( badIDs$ID )
  cat( "Missing SE estimates:\n" )
  print( badIDs )
  #sruns = filter( sruns, !is.na( mean.SE.hat ) & !is.na( RMSE.f.stu ) )  
  stop()
}
nrow(sruns)
table(sruns$method)


# missing IDs from prior filtering
if ( nrow( failed ) > 0 ) {
  cat( "Failed runs--looking at where they showed up" )
  badIDs = unique( failed$ID )
  print( badIDs )
  filter( scenarios, ID %in% badIDs )
  stop()
}

#which are the scenarios where the FI SE ran 
filter(scenarios, !(ID %in% badIDs ))
#in the lowdf sims J.Bar has to be at least 4, anyway this is a thing to come back to in the memo 

# Which IDs are missing stats for some of the methods?
ss = sruns %>% ungroup() %>% dplyr::select( ID, method, mean.SE.hat ) %>%
  spread( method, mean.SE.hat )
head( ss )
not.comp = ss[ !complete.cases(ss), ]

temp = merge( scenarios, not.comp, by="ID" )
if ( nrow( not.comp ) > 0 ) {
  warning( "Not complete cases found\n" )
  print( not.comp )
  stop()
}

nrow(sruns)
table(sruns$method)

##### calculate the finite sample SE using our nested simulations  #####

subruns = simdata %>% 
  group_by( ID, method, run ) %>%            # just aggregating over subruns
  summarise( nsr = n(),
             SE2 = var( ATE.hat ),           # noisy estimate of true finite-sample SE^2
             # SE.hat.f2 = mean( SE.hat^2 ),
             var.SE.hat = var( SE.hat ) )    # variation in estimated SE
subruns <- subruns %>% 
  filter(nsr == 3) %>%                                  # need all three subruns
  group_by( ID, method ) %>%                            # now aggregate over runs
  summarise( SE.f = sqrt( mean( SE2 ) ),                # true finite-sample SE^2
             #SE.hat.f = sqrt( mean( SE.hat.f2 ) ),
             sd.SE.hat.f = sqrt( mean( var.SE.hat ) ) ) # average sd of SE estimates

# add to sruns
sruns = merge( sruns, subruns, by=c("ID","method"), all=TRUE )
rm( subruns )


head( sruns )
nrow(sruns)
table(sruns$method)


# Check: How do the true SEs of superpop compare to the finite pop
# (These should generally be above 1.0)
summary( sruns$SE.s / sruns$SE.f )
odds = sruns %>% mutate( ratio = SE.s / SE.f ) %>%
  filter( ratio < .99 )
table( odds$ID )
filter( scenarios, ID %in% odds$ID )
head(odds) %>% as.data.frame()

# inflate captures how well the estimated standard error predicts the true
# (averaged across finite samples for finite) SE.
sruns = mutate( sruns, 
                inflate.f = mean.SE.hat / SE.f,
                inflate.s = mean.SE.hat / SE.s )


head( sruns )
nrow(sruns)
table(sruns$method)

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
#manually fix group because I didn't in the simulation program, this doesn't effect anything but the label 
sruns$group<-"lowdf"
scenarios$group<-"lowdf"

group=unique(sruns$group)
#group=unique(simdata$group)

#save the results 
#for now just save the aggregate results in an agg folder

file_name=paste("results_agg/processed_simulation_results_agg_",group,"_",Sys.Date(),".RData",sep="")


#rename sruns so that it can be combined with other sruns data 
sruns.name<-paste("sruns",group,sep=".")
assign(sruns.name,sruns)

save(scenarios,list=c(sruns.name), file=file_name )


#old save code 
  #file_name=paste("results/processed_simulation_results_",group,"_",Sys.Date(),".RData",sep="")
  #file_name=paste("processed_simulation_results_",group,"_",Sys.Date(),".RData",sep="")
  #rename the range data 
  #range.dat.name<-paste("range.dat",group,sep=".")
  #assign(range.dat.name,range.dat)

  #save(scenarios, simdata, list=c(sruns.name,range.dat.name), file=file_name )
  
    

