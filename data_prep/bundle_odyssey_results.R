
# 
# Take all the RDS files from Odyssey and bundle them up into a single results file
#
# These will be saved in 'complete_odyssey_results.rds'
#
# Also runs checks on how long things took to run.

library( tidyverse )


# Load a collection of files and put them all together into a single, very large
# dataframe
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


results = load.all.sims( filehead="raw_results/" )
#res2 = load.all.sims( filehead = "od_res_run1/" )
#results
#results = bind_rows( results, res2, .id="batch" )
head( results )
nrow( results )
length( unique( results$filename ) )
#table( results$batch )

#toplevel.sim.per.min = results$sim.per.min 
#results$sim.per.min = NULL
#summary( toplevel.sim.per.min )
#qplot( toplevel.sim.per.min )

results$batch = as.numeric( as.factor( results$filename ) )
table( table( results$ID ) )


##### Merge all reps of a scenario together ######

cat( "Merging results by ID\n" )

head( results$run[[1]] )
results$run[[4]]$run
results$batch[[4]]

head( results )

if ( FALSE ) {
  r.100 = filter( results, ID == 200 )
  summary( r.100$ATE.super.person )
  qplot( r.100$ATE.super.person )
  filter( r.100, ATE.super.person < 0.202 )$filename
  head( r.100 )
  head( r.100$run[[1]] )
  r.100$filename[[1]]
  head( r.100$run[[2]] )
  r.100$filename[[2]]
  r.100$run[[1]]$ATE.hat == r.100$run[[2]]$ATE.hat
  r.100$run[[1]]$ATE.hat == r.100$run[[50]]$ATE.hat
}

#results %>% group_by( ID ) %>% 
#  summarise( ATE.super.person = )

#pb<- txtProgressBar(min=0,max=nrow(results),style=3)
#for ( r in 1:nrow( results ) ) {
#  results$run[[r]]$run = paste( results$batch[[r]], results$run, sep="." )
#  results$subrun[[r]]$run = paste( results$batch[[r]], results$subrun, sep="." )
#  setTxtProgressBar(pb,r)
#}
#close(pb)

results = unnest( results )
nrow( results )
results = mutate( results, run = paste( batch , run, sep="." ),
                  subrun = paste( batch, subrun, sep="." ) )
results
head( results$subrun )
head( results$run )

results = results %>% dplyr::select( -filename, -batch ) %>%
  group_by( ID ) %>%
  mutate( ATE.super.person = mean( ATE.super.person ),
          sim.per.min = mean( sim.per.min ),
          R = length( unique( run ) ) ) %>% ungroup() %>%
  group_by( ID, J, n.bar, tau, dependence, proptx.dependence, variable.n, variable.p, ATE.super.site, ATE.super.person, sim.per.min, R ) %>%
    nest()

results = rename( results, run = data )
results

filter( results, duplicated( results$ID ) )

saveRDS( results, file="results/complete_odyssey_results.rds" )



##### Check that we got all scenarios and how many reps for each scenario #####

source( "simulation_support_code.R" )

lst = make.scenario.list()
lst$ID = 1:nrow(lst)
setdiff( lst$ID, results$ID )
setdiff( results$ID, lst$ID )
nrow( results )
nrow( lst )

head( results$run[[1]] )

results$R = map_dbl( results$run, function( dd ) { length( unique( dd$run ) ) } )

results = arrange( results, ID )
results

table( results$ATE.super.site )

summary( results$R )


##### Make overall runtime stats for each scenario ######

simstats = results %>% dplyr::select( -run )
simstats

#ggplot( simstats, aes( R ) ) +
#  facet_wrap( ~ J ) +
#  geom_histogram()

quantile( simstats$sim.per.min, 0.05 )
slow = filter( simstats, sim.per.min < quantile( sim.per.min, 0.05 ) )
slow
slow$J * slow$n.bar

summary( simstats$sim.per.min )

simstats$sec.per.sim = 60 / simstats$sim.per.min

# Calculate how long 1000 trials would run, expected
timings = 1000 * simstats$sec.per.sim / (60 * 60)
summary( timings )

simstats$n = simstats$n.bar * simstats$J
M0 = lm( log(sec.per.sim) ~ I(J/5)  + I(n/100) + dependence + proptx.dependence + tau, data=simstats )
summary( M0 )
exp( coef( M0 ) )

timings = simstats %>% group_by( J, n ) %>% summarise( low.time = quantile( sec.per.sim, 0.90 ),
                                            median.time = median( sec.per.sim ) )
ggplot( timings, aes( n, low.time, col=J, group=J ) ) + geom_line() + geom_point()


simstats = mutate( simstats,
                         e.time = 1.13^(J/5) +	1.0183^(n/100) )
summary( simstats$e.time )

qplot( e.time, sec.per.sim, data=simstats ) +
  geom_abline(intercept = 0, slope = 1 )



