
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

#rm(list=ls())
#library( tidyverse )
#load("results_sav/processed_simulation_results_estimators_estimand_2024-09-30.RData")
#load("results_sav/processed_simulation_results_sim_base_1_2024-09-30.RData")
#load("results_sav/processed_simulation_results_estimators_estimand_2_2024-09-30.RData")

# Drop simulations from the scenarios frame
scenarios$run = NULL


## sometimes there is no SE.hat becuase SE.hat is not defined 

#failed <- simdata %>% filter(is.na(method))
failed <- simdata %>% filter (is.na(SE.hat))
#failed <- simdata %>% filter(is.infinite(SE.hat))

#if ( nrow( failed ) > 0 ) {
 # warning( "Dropping some runs without estimated SEs" )
  #simdata = filter( simdata, !(is.na( SE.hat ) | is.infinite( SE.hat )) )
#}


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
                 cov.f.clu = cover.vector(ATE.hat, SE.hat, df, ATE.finite.cluster),
                 cov.f.dstu = cover.vector(ATE.hat, SE.hat, df, ATE.finite.dist.student),
                 cov.f.dclu = cover.vector(ATE.hat, SE.hat, df, ATE.finite.dist.cluster),
                 
                 # Coverage for super populations
                 cov.s.stu = cover.vector(ATE.hat, SE.hat, df, ATE.super.student),
                 cov.s.clu = cover.vector(ATE.hat, SE.hat, df, ATE.super.cluster),
                 cov.s.dstu = cover.vector(ATE.hat, SE.hat, df, ATE.super.dist.student),
                 cov.s.dclu = cover.vector(ATE.hat, SE.hat, df, ATE.super.dist.cluster),
                 
                 # CI length
                 CI.length = CI.length(ATE.hat, SE.hat, df))


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

sruns = simdata %>% 
  group_by( ID, method, weight, disfavored, biased) %>%   # per scenario and method... (aggregate over batch, run, subrun)
  summarise( mean.SE.hat = mean( SE.hat ),
             mean.ATE.hat=mean(ATE.hat),
             mean.df=mean(df),
             # bias of estimators for each estimand
             bias.f.stu = mean( ATE.hat - ATE.finite.student ),
             bias.f.clu = mean( ATE.hat - ATE.finite.cluster ),
             bias.f.dstu = mean( ATE.hat - ATE.finite.dist.student ),
             bias.f.dclu = mean( ATE.hat - ATE.finite.dist.cluster ),
             bias.s.stu = mean( ATE.hat - ATE.super.student ),
             bias.s.clu = mean( ATE.hat - ATE.super.cluster ),
             bias.s.dstu = mean( ATE.hat - ATE.super.dist.student ),
             bias.s.dclu = mean( ATE.hat - ATE.super.dist.cluster ),
             
             # RMSE of estimators for each estimand
             RMSE.f.stu = sqrt( mean( ( ATE.hat - ATE.finite.student )^2 ) ),
             RMSE.f.clu = sqrt( mean( ( ATE.hat - ATE.finite.cluster   )^2 ) ),
             RMSE.f.dstu = sqrt( mean( ( ATE.hat - ATE.finite.dist.student )^2 ) ),
             RMSE.f.dclu = sqrt( mean( ( ATE.hat - ATE.finite.dist.cluster   )^2 ) ),
             RMSE.s.stu = sqrt( mean( ( ATE.hat - ATE.super.student  )^2 ) ),
             RMSE.s.clu = sqrt( mean( ( ATE.hat - ATE.super.cluster )^2 ) ),
             RMSE.s.dstu = sqrt( mean( ( ATE.hat - ATE.super.dist.student  )^2 ) ),
             RMSE.s.dclu = sqrt( mean( ( ATE.hat - ATE.super.dist.cluster )^2 ) ),
             
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
             cover.f.clu = mean( cov.f.clu ),
             cover.f.dstu = mean( cov.f.dstu ),
             cover.f.dclu = mean( cov.f.dclu ),
             cover.s.stu = mean( cov.s.stu ),
             cover.s.clu = mean( cov.s.clu ),
             cover.s.dstu = mean( cov.s.dstu ),
             cover.s.dclu = mean( cov.s.dclu ),
             CI.length=mean(CI.length))
             
             
             # var.err2.s.stu = var( ATE.hat - ATE.super.student ),
             # var.err2.s.sch = var( ATE.hat - ATE.super.school ),
             # var.err2.s.dstu = var( ATE.hat - ATE.super.dist.student ),
             # var.err2.s.dsch = var( ATE.hat - ATE.super.dist.school )

# Add in the info on the scenario to the results.
sruns = merge( scenarios, sruns, by="ID" )
head( sruns )


##### Ensure no bad sim runs #####

# look for missing SE estimates
badIDs = filter( sruns, is.na( mean.SE.hat )) 
#badIDs = filter( sruns, is.na( RMSE.f.stu  ))

#we don't want to drop these, because sometimes ATE is defined and SE isn't 
#if ( nrow( badIDs ) > 0 ) {
 # warning( "Continuing to drop runs with missing estimated SEs or RMSEs" )
  #badIDs = unique( badIDs$ID )
  #cat( "Missing SE estimates:\n" )
  #print( badIDs )
  #sruns = filter( sruns, !is.na( mean.SE.hat ) & !is.na( RMSE.fp ) )  
  #stop()
#}


# missing IDs from prior filtering, but I don't think we actually care about this
#if ( nrow( failed ) > 0 ) {
 # cat( "Failed runs--looking at where they showed up" )
  #badIDs = unique( failed$ID )
  #print( badIDs )
  #temp<-filter( scenarios, ID %in% badIDs )
  #print(temp)
  #stop()
#}

# Which IDs are missing stats for some of the methods?
#again sometimes SE's are missing for the DB methods, which is just a known problem
ss = sruns %>% ungroup() %>% dplyr::select( ID, method, mean.SE.hat ) %>%
  spread( method, mean.SE.hat )
head( ss )
not.comp = ss[ !complete.cases(ss), ]

temp = merge( scenarios, not.comp, by="ID" )
if ( nrow( not.comp ) > 0 ) {
  warning( "Not complete cases found\n" )
  print( not.comp )
  # stop()  # disabled 2026-06-10: the incomplete cases are the known DB-method gap
  #         # (some estimators don't return an SE on harder designs). Per the comments
  #         # above we keep these rather than drop them, so this is expected, not an
  #         # error. The warning() and print() above still log it for the record.
}


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
  
    

