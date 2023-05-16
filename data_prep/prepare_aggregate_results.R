

##
## Calculate summary stats of the scenarios for the markdown file.
##
## Takes in processed sim results after they have been aggregated and cleaned up
## a bit.  Saves a much smaller datafile with aggregate data.
##
## These will be saved in 'aggregated_simulation_results.RData'
##
## This script also has many side calculations and investigations that are used
## to enrich the Rmd report.
##

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd( ".." )

library( tidyverse )

library( ggthemes )

my_t = theme_tufte() + theme( legend.position="bottom", 
                              legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                              panel.border = element_blank() )
theme_set( my_t )

cat( "Loading the big data file\n" )
load( "results/processed_simulation_results.RData" )


cat( "Now processing\n" )


# Get characteristics of the different methods
mc = blkvar::method.characteristics()
mc

table( sruns$method )

# Select methods to get the unique point estimates
selected = c( "DB-FP-Persons", "DB-FP-Sites", "FE", "RICC", "FIRC", "RIRC" )
core.selected = c( "DB-Persons", "DB-Sites", "FE", "FIRC" )
cat( "Selecting:\n" )
print( selected )
cat( "Core Selecting:\n" )
print( core.selected )


# NOTE ON VARIABLES:
# sruns - dataframe of aggregated results.  These were calculated in prepare_aggregate_results
head( sruns )
nrow( sruns ) / length( unique( sruns$method ) )



# For each scenario, calculate relative performance rates using the design based
# estimators as reference categories
calc.stats = function( runs.m ) {
  m.fp = filter( runs.m, method=="DB-Persons" )$RMSE.fp
  m.sp = filter( runs.m, method=="DB-Persons" )$RMSE.sp
  m.fs = filter( runs.m, method=="DB-Sites" )$RMSE.fs
  m.ss = filter( runs.m, method=="DB-Sites" )$RMSE.ss
  #fbase = filter( runs.m, method=="DB-FP-Persons" )
  fbase = filter( runs.m, method=="DB-Persons" )
  Var.s.base = fbase$SE.s^2 / ( 2 * (fbase$R-1) )
  
  #mns = filter( mns, method %in% selected )
  #  mns = dplyr::select( runs.m, method, tau, dependence, proptx.dependence, weight, population, 
  #                       SE.f, SE.s, RMSE.fp, RMSE.sp, RMSE.fs, RMSE.ss, bias.sp, bias.ss ) %>%
  mns = mutate( runs.m, 
                bias.p.pct = 100 * bias.sp^2 / RMSE.sp^2,
                bias.s.pct = 100 * bias.ss^2 / RMSE.ss^2,
                RMSE.fp = 100 * RMSE.fp / m.fp,
                RMSE.sp = 100 * RMSE.sp / m.sp,
                RMSE.fs = 100 * RMSE.fs / m.fs,
                RMSE.ss = 100 * RMSE.ss / m.ss,
                SE.s.inf = 100 * SE.s / SE.f,
                SE.f = 100 * SE.f / fbase$SE.f,
                SE.SE.s = 100 * sqrt( (1/fbase$SE.s^2) * (SE.s^2 / (2*(R-1))) + (SE.s^2/fbase$SE.s^4) * Var.s.base ),
                SE.s = 100 * SE.s / fbase$SE.s
  )
  mns = rename( mns, 
                bias.per = bias.sp,
                bias.site = bias.ss )
  mns
}


# Change method names to names without reference to superpop vs. finite pop
refactor_methods = function( dat ) {
  mutate( dat, method = factor( method, 
                                levels=c( "DB-FP-Persons", "DB-FP-Sites", "FE", "RICC", "FIRC", "RIRC" ),
                                labels=c( "DB-Persons", "DB-Sites", "FE", "RICC", "FIRC", "RIRC" ) ) )
}


# For making tables with intervals
get.interval = function( vals, digits=0 ) {
  qt = quantile( vals, c( 0.05, 0.95 ) )
  if ( diff( range( qt, 100 ) ) < 0.00001 ) {
    "Pig"
  } else {
    #sprintf( "%.0f (%.0f-%.0f)", mean( vals ), qt[[1]], qt[[2]] )
    sprintf( paste0( "(%.", digits, "f-%.", digits, "f)" ), qt[[1]], qt[[2]] )
  }
}


#### Calculate performance statistics of the main estimators for all scenarios ####


sruns.core = sruns %>% 
  filter( method %in% selected ) %>%
  refactor_methods()

mns = sruns.core %>%
  group_by( ID ) %>% do( calc.stats(.) ) %>% ungroup()

head( mns )
nrow( mns ) / length( unique( mns$method ) )
table( mns$method )

mns.full = mns
mns = dplyr::filter( mns, method %in% core.selected )
mns$method = droplevels( mns$method )
table( mns$method )


#### Divide scenarios into variable (main) scenarios and the rest  #####

mns.novar = filter( mns, variable.p == FALSE | variable.n == FALSE )
#mns = filter( mns, variable.p & variable.n )

sruns.novar = filter( sruns, variable.p == FALSE | variable.n == FALSE )
#sruns = filter( sruns, variable.p & variable.n )



####  True SE performance metrics ####


cat( "Now processing SE performance metrics\n" )

table( mns$method )
SE.table = mns %>% filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( finite = mean( SE.f ),
             range.f = get.interval( SE.f ),
             super = mean( SE.s ),
             range.s = get.interval( SE.s ),
             inflation = mean( SE.s.inf ),
             range.i = get.interval( SE.s.inf ) )

SE.table[ SE.table=="Pig" ] = NA
SE.table


filter( mns, method=="FIRC", SE.f > 115 )
mns$n = mns$n.bar * mns$J 
ggplot( filter( mns, method=="FIRC" ), aes( n, SE.f, pch=dependence, col=proptx.dependence ) ) +
  facet_grid( J ~ tau ) +
  geom_point()


filter( mns, method=="FIRC", tau==0, J == 20 )

M0 = lm( SE.f ~ dependence + proptx.dependence + J + n.bar + n + tau, data=filter( mns, method=="FIRC" ) )
summary( M0 )
round( coef( M0 ), digits=1 )


# Looking at simulation uncertainty for the table
head( mns )
summary( mns$SE.SE.s )

mns %>% filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( finite = mean( SE.f ),
             SE = sd( SE.f ) / sqrt( n() ),
             super = mean( SE.s ),
             SE.s = sd( SE.s ) / sqrt( n() ) ) %>%
  mutate_at( c("SE","SE.s"), round, digits=3 )


#### Looking at scenarios with no proportion or site size variation ####


# Sanity check that no variation means point estimates are all the same.  They are
simdata.const = filter( simdata, variable.p == FALSE, variable.n == FALSE ) %>%
          filter( method %in% selected )
head( simdata.const )
dbs = simdata.const %>% filter( J == 40, n.bar == 200 ) %>%
  dplyr::select( ID, subrun, method, ATE.hat ) %>% spread( key=method, value=ATE.hat )
head( dbs )


mns.const = sruns.novar %>%   filter( method %in% selected ) %>%
  refactor_methods() %>%
  group_by( ID ) %>% do( calc.stats(.) ) %>% ungroup()
mns.const
nrow( mns.const )

table( mns.const$J, mns.const$n.bar )
table( mns.const$tau )
table( mns.const$ID )
length( unique( mns.const$ID ) )

# Should all be near 0.20
summary( mns.const$ATE.super.person )

head( mns.const )

SE.table.const = mns.const %>% group_by( method ) %>%
  summarise( finite = mean( SE.f ),
             range.f = get.interval( SE.f ),
             super = mean( SE.s ),
             range.s = get.interval( SE.s ),
             inflation = mean( SE.s.inf ),
             range.i = get.interval( SE.s.inf ) )

SE.table.const[ SE.table.const=="Pig" ] = NA
SE.table.const
  
if ( FALSE ) {
  tt = filter( mns, method == "DB-Sites" )
  nrow( tt )
  head( tt )
  aa = dplyr::select( tt, method, tau, dependence, proptx.dependence )
  aa = unique( aa )
  aa
  nrow( aa )
}



#### RMSE ####
cat( "Now processing RMSE performance metrics\n" )

RMSE.table = mns %>% filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( mn.fp = mean( RMSE.fp ),
             rng.fp = get.interval( RMSE.fp ),
             #sd.RMSE.fp = sd( RMSE.fp ),
             mn.sp = mean( RMSE.sp ),
             rng.sp = get.interval( RMSE.sp ),
             #sd.RMSE.sp = sd( RMSE.sp ),
             mn.fs = mean( RMSE.fs ),
             rng.fs = get.interval( RMSE.fs ),
             #    rng.fs = get.interval( RMSE.fs ),
             #sd.RMSE.fs = sd( RMSE.fs ),
             mn.ss = mean( RMSE.ss ),
             rng.ss = get.interval( RMSE.ss ) )
#sd.RMSE.ss = sd( RMSE.ss ) )

RMSE.table[ RMSE.table == "Pig" ] = NA
RMSE.table
#for ( i in 2:ncol(mns.agg) ) {
#  mns.agg[ mns.agg[[i]] < 0.0001, i ] = NA
#}

if (FALSE) {
# Looking at ranges
mns2 = gather(mns, RMSE.fp, RMSE.sp, RMSE.fs, RMSE.ss, key="estimand", value="RMSE" )
head( mns2 )
mns2 = mns2 %>% group_by( method,estimand ) %>% 
  summarise( q20 = quantile( RMSE, 0.20 ),
             q80 = quantile( RMSE, 0.80 ) )
head( mns2 )
filter( mns2, method=="FIRC" )
}


#names(mns.agg) = gsub( "mean.", "", names(mns.agg ) )
#for ( i in c(8,6,4,2 ) ) {
#  mns.agg[,i] = paste0( "$", round( mns.agg[,i,drop=TRUE] ), " \\pm ", round( 2 * mns.agg[,i+1, drop=TRUE] ), "$"  )
#}
#mns.agg = mns.agg[ c(1,2,4,6,8) ]
#mns.agg[ mns.agg == "$100 \\pm NA$" ] = "100"


# SE for the mean RMSE estimates
mns %>% filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( mn.fp = mean( RMSE.fp ),
             SE.fp = sd( RMSE.fp ) / sqrt(n()),
             mn.sp = mean( RMSE.sp ),
             SE.sp = sd( RMSE.sp ) / sqrt(n()),
             mn.fs = mean( RMSE.fs ),
             SE.fs = sd( RMSE.fs ) / sqrt(n()),
             mn.ss = mean( RMSE.ss ),
             SE.ss = sd( RMSE.ss ) / sqrt(n()) ) %>%
  mutate_if( is.numeric, round, digits=2 )


#### Bias ####

cat( "Now processing bias performance metrics\n" )

#mns = sruns %>% filter( method %in% selected ) %>% refactor_methods() %>%
#  group_by( method, weight, tau, proptx.dependence, dependence ) %>%
#  summarise( mean.bias.fp = mean( bias.fp ),
#             sd.bias.fp = sd( bias.fp ),
#             mean.bias.fs = mean( bias.fs ),
#             sd.bias.fs = sd( bias.fs ) )

get.interval.bias = function( vals ) {
  qt = quantile( vals, c( 0.1, 0.9 ) )
  if ( diff( range( qt, 100 ) ) < 0.00001 ) {
    "Pig"
  } else {
    #sprintf( "%.0f (%.0f-%.0f)", mean( vals ), qt[[1]], qt[[2]] )
    sprintf( "(%.2f-%.2f)", qt[[1]], qt[[2]] )
  }
}

bias.table = mns %>% filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( mean.bias.per = mean( bias.per ),
             rng.bias.per = get.interval.bias( bias.per ),
             sd.bias.per = sd( bias.per ),
             mean.bias.site= mean( bias.site ),
             rng.bias.site = get.interval.bias( bias.site ),
             sd.bias.site = sd( bias.site ) )

bias.table$sd.bias.per = bias.table$sd.bias.site = NULL

bias.table


summary( mns$bias.p.pct )

bias.table.pct = mns %>% group_by( method ) %>%
  summarise( mean.bias.per = mean( bias.p.pct ),
             rng.bias.per = get.interval( bias.p.pct ),
             mean.bias.site = mean( bias.s.pct ),
             rng.bias.site = get.interval( bias.s.pct ) )

bias.table.pct


#### Examining RICC with respect to FE   and RIRC with respect to FIRC ####

head( simdata )
table( simdata$method )
csim = filter( simdata, method %in% c( "RICC", "FE", "RIRC", "FIRC" ), variable.n == TRUE, variable.p == TRUE )

head( csim )
csim$SE.hat = NULL
csim = spread( csim, method, ATE.hat )
head( csim )

cor( csim$FE, csim$RICC )
round( summary( csim$FE - csim$RICC ), digits= 3 )
mean( abs( csim$FE - csim$RICC ) > 0.01 )

cor( csim$FE, csim$FIRC )
round( summary( csim$FE - csim$FIRC ), digits= 3 )

cor( csim$RIRC, csim$FIRC )
round( summary( csim$RIRC - csim$FIRC ), digits= 3 )
sd( csim$RIRC - csim$FIRC )
quantile( csim$RIRC - csim$FIRC, c( 0.05, 0.95 ) )
mean( abs( csim$RIRC - csim$FIRC ) >= 0.01 )
mean( abs( csim$RIRC - csim$FIRC ) >= 0.02 )

# How do RMSEs compare?
cmns = filter( mns.full, method %in% c( "RICC", "FE", "RIRC", "FIRC" ) ) %>%
  dplyr::select( method, ID, RMSE.fp ) %>%
  spread( method, RMSE.fp )
head( cmns )
nrow( cmns )
summary( cmns$RICC / cmns$FE )
t.test( cmns$RICC / cmns$FE )
qplot( cmns$RICC / cmns$FE )

nrow( cmns )
summary( cmns$RIRC / cmns$FIRC )
t.test( cmns$RIRC / cmns$FIRC )
qplot( cmns$RIRC / cmns$FIRC )



#### Sanity checking: RMSE^2 = bias^2 + SE^2 ####

head( sruns )
ss = sruns %>% mutate( tt.f = bias.fp^2 + SE.f^2,
                       tt.s = bias.sp^2 + SE.s^2, 
                       tt.ss = bias.ss^2 + SE.s^2,
                       RMSE.check.f = tt.f - RMSE.fp^2,
                       RMSE.check.ss = tt.ss / RMSE.ss^2,
                       birat = bias.sp^2 / RMSE.sp^2 )

summary( ss$RMSE.check.ss )
filter( ss, RMSE.check.ss > 8 ) %>% dplyr::select( method, ID, bias.ss, RMSE.ss, SE.s )

if (FALSE) {
  sss = filter( simdata, ID == 1041, method=="FE" )
  nrow( sss )
  SE2.true = var( sss$ATE.hat )
  unique( sss$ATE.super.site )
  bias.true = mean( sss$ATE.hat ) - 0.2  
  bias.true
  RMSE.true = 4
}

wack = filter( ss, birat > 1 )
nrow( wack )
summary( wack$birat )
head( wack )

round( summary( ss$tt.f / ss$RMSE.fp^2 ), digits=3 )

round( summary( ss$tt.s / ss$RMSE.fs^2 ), digits=3 )

# Check results: Close, but not perfect.



# Looking at SE vs RMSE
head( sruns )

summary( 100 * sruns$SE.f^2 / sruns$RMSE.fp^2)

summary( 100 * sruns$SE.s^2 / sruns$RMSE.sp^2)






##### Comparing the superpop SE.hat to finite pop SE.hat #####

cat( "Now processing SE.hat comparison metrics\n" )

# How often are superpop SEs lower than finite pop SEs?
# Comparing FIRC to DB
head( simdata )
unique( simdata$method )
ssub = filter( simdata, method %in% c( "FIRC", "DB-FP-Persons", "DB-FP-Sites", "FE-CR" ) ) %>% 
  dplyr::select( ID, J, n.bar, tau, dependence, subrun, method, SE.hat )
head( ssub )
ssub = mutate( ssub, subrun = paste( ID, subrun, sep="." ) )
ssub = spread( ssub, key=method, value=SE.hat )

head( ssub )
nrow( ssub )
# Are SE estimates correlated?

# This divided by dependance---but that seems unneeded.
#ssub = ssub %>% group_by( J, n.bar, tau, dependence ) %>% 
#  summarise( FIRC.less = mean( FIRC < `DB-FP-Persons` ),
#             FE.CR.less = mean( `FE-CR` < `DB-FP-Persons` ),
#             FIRC.way.less = mean( FIRC < `DB-FP-Persons` - 0.025 ),
#             FE.CR.way.less = mean( `FE-CR` < `DB-FP-Persons` - 0.025 ),
#             )
#ggplot( ssub, aes( tau, percent, col=method, lty=dependence ) ) +
#  facet_grid( J ~ n.bar, labeller = label_both ) +
#  geom_point() + geom_line()
ssub = ssub %>% group_by( J, n.bar, tau ) %>% 
  summarise( FIRC.less = mean( FIRC < `DB-FP-Persons` ),
             FE_CR.less = mean( `FE-CR` < `DB-FP-Persons` ),
             FIRC.way_less = mean( FIRC < `DB-FP-Persons`*0.90 ),
             FE_CR.way_less = mean( `FE-CR` < `DB-FP-Persons`*0.90 ),
             FIRC.less_site = mean( FIRC < `DB-FP-Sites` ),
             FIRC.way_less_site = mean( FIRC < `DB-FP-Sites`*0.90 )   )

head( ssub )
ssub2 = gather( ssub, FIRC.less, FE_CR.less, FIRC.way_less, FE_CR.way_less, FIRC.less_site, FIRC.way_less_site, key="method", value="percent" )
head( ssub2 )
table( ssub2$method )
ssub2 = separate( ssub2, method, into=c("method","threshold"), sep="\\." )
table(ssub2$method )
table(ssub2$threshold )
ssub2
SE.comp.results = ssub2



##### Looking at correlation of the SE hats ######

cat( "Now processing correlation of SE hats for specific scenario\n" )

head( simdata )
unique( simdata$method )
ssub = simdata %>% filter( J %in% c( 20, 40 ), n.bar == 200, dependence == "indep", proptx.dependence == "indep" ) %>%
  filter( method %in% c( "FIRC", "DB-FP-Persons", "FE-CR", "FE-Club", "DB-SP-Sites" ) ) %>% 
  dplyr::select( ID, J, n.bar, tau, dependence, proptx.dependence, run, subrun, method, SE.hat )
head( ssub )
ssub = spread( ssub, key=method, value=SE.hat )
head( ssub )
ssub = separate(ssub, subrun, into=c("batch", "rr", "scount" ) )
head( ssub )

#s2 = filter( ssub, J  %in% c(20,80), n.bar == 100, tau %in% c(0, 0.2), dependence == "corr", proptx.dependence == 0 )
s2 = filter( ssub, scount == 1)
table( s2$ID )
head( s2 )
nrow( s2 )
s2 = mutate( s2, smaller.FIRC = ifelse( 0.99*`DB-FP-Persons` > `FIRC`, "flip", "ok" ),
             smaller.CR = ifelse( 0.99*`DB-FP-Persons` > `FE-CR`, "flip", "ok" ),
             small = paste0( smaller.CR, "-", smaller.FIRC ) )

SE.corr.results = s2
head( SE.corr.results )
#rm( simdata, s2, ssub2, ssub )
#rm( scenarios )

cat( "Saving all the processed results\n" )

save( bias.table, bias.table.pct, RMSE.table, SE.table, SE.comp.results, SE.corr.results, mns, sruns, sruns.core,
      get.interval, get.interval.bias, file = "results/aggregated_simulation_results.RData" )




