
# Take all the RDS files from Odyssey and bundle them up into a single results file
#
# These will be saved in 'complete_odyssey_results.rds'
#
# Also runs checks on how long things took to run.

library( tidyverse )


# Load a collection of files and put them all together into a single, 
# very large dataframe
load.all.sims = function( filehead="results/" ) {
  
  files = list.files( filehead, full.names=TRUE)
  
  res = map_df( files, function( fname ) {
    cat( "Reading results from ", fname, "\n" )
    rs = readRDS( file = fname )
    rs$filename = fname
    rs
  })
  res
}


results = load.all.sims( filehead="results/" )


##### Merge all reps of a scenario together ######

cat( "Merging results by ID\n" )

# there are a bunch of ID columns, here's a summary:
#  - ID: which simulation scenario (1 to nrow(make.scenario.list()))
#  - batch: which slurm job (1 to a bajillion)
#  - run: which run of the scenario (1 to NUM_TRIALS_PER_SCENARIO = 4)
#  - subrun: which rerandomization of the run (1 to 3)

results <- results %>% 
  mutate(batch = as.numeric( as.factor( results$filename ))) %>% 
  unnest(cols = run) %>% 
  mutate(run = paste( batch , run, sep="." ),
         subrun = paste( batch, subrun, sep="." )) %>% 
  dplyr::select( -filename, -batch )

results <- results %>%
  group_by(ID) %>%
  mutate(across(starts_with("ATE.super"), mean),   # true superpop ATEs are scenario-level averages
         sim.per.min = mean( sim.per.min ),        # timing of each scenario
         R = length( unique( run ) )) %>%          # ??? used later though
  ungroup() %>%
  
  # regroup data by scenario
  group_by( ID, 
            # keep all scenario-level data-generating parameters
            K, J.bar, p.bar, sd.district, ICC.district, variable.J, variable.p, 
            n.bar, sd.school, ICC.school, variable.n,
            size.impact.cor.district, size.impact.cor.school, proptx.impact.cor.district,
            # superpop ATEs are scenario-level attributes, keep them
            ATE.super.district, ATE.super.school, ATE.super.student,
            ATE.super.dist.school, ATE.super.dist.student,
            # scenario timing is scenario-level attribute, keep it
            sim.per.min, R) %>% 
  nest() %>% 
  rename(run = data)

saveRDS( results, file="results/complete_odyssey_results.rds" )




# TODO: update below, e.g., instead of J*n.bar it's now K*J.bar*n.bar, etc.

# ##### Check that we got all scenarios and how many reps for each scenario #####
# 
# source( "simulation_support_code_crt.R" )
# 
# lst = make.scenario.list(group="initial")
# lst$ID = 1:nrow(lst)
# setdiff( lst$ID, results$ID )
# setdiff( results$ID, lst$ID )
# nrow( results )
# nrow( lst )
# 
# head( results$run[[1]] )
# 
# results$R = map_dbl( results$run, function( dd ) { length( unique( dd$run ) ) } )
# 
# results = arrange( results, ID )
# results
# 
# table( results$ATE.super.site )
# 
# summary( results$R )
# 
# 
# ##### Make overall runtime stats for each scenario ######
# 
# simstats = results %>% dplyr::select( -run )
# simstats
# 
# #ggplot( simstats, aes( R ) ) +
# #  facet_wrap( ~ J ) +
# #  geom_histogram()
# 
# quantile( simstats$sim.per.min, 0.05 )
# slow = filter( simstats, sim.per.min < quantile( sim.per.min, 0.05 ) )
# slow
# slow$J * slow$n.bar
# 
# summary( simstats$sim.per.min )
# 
# simstats$sec.per.sim = 60 / simstats$sim.per.min
# 
# # Calculate how long 1000 trials would run, expected
# timings = 1000 * simstats$sec.per.sim / (60 * 60)
# summary( timings )
# 
# simstats$n = simstats$n.bar * simstats$J
# M0 = lm( log(sec.per.sim) ~ I(J/5)  + I(n/100) + dependence + proptx.dependence + tau, data=simstats )
# summary( M0 )
# exp( coef( M0 ) )
# 
# timings = simstats %>% group_by( J, n ) %>% summarise( low.time = quantile( sec.per.sim, 0.90 ),
#                                             median.time = median( sec.per.sim ) )
# ggplot( timings, aes( n, low.time, col=J, group=J ) ) + geom_line() + geom_point()
# 
# 
# simstats = mutate( simstats,
#                          e.time = 1.13^(J/5) +	1.0183^(n/100) )
# summary( simstats$e.time )
# 
# qplot( e.time, sec.per.sim, data=simstats ) +
#   geom_abline(intercept = 0, slope = 1 )



