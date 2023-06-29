

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


library( tidyverse )

library( ggthemes )

my_t = theme_tufte() + theme( legend.position="bottom", 
                              legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                              panel.border = element_blank() )
theme_set( my_t )

cat( "Loading the big data file\n" )
load( "results/processed_simulation_results.RData" )

cat( "Now processing\n" )



# selected is just all methods
table( sruns$method )
simdata %>% ungroup() %>% 
  select(ID, subrun, method, ATE.hat)
selected = c( "DB_clust (FE)", "DB_clust (int)", 
              "DB_indiv (FE)", "DB_indiv (int)",
              "LR (agg)", "LR (agg, wt)", "LR w/ CRSE", 
              "MLM" )
# core.selected have unique point estimates
core.selected = c( "DB_clust (int)", 
                   "DB_indiv (int)",
                   "LR (agg)", "LR (agg, wt)", 
                   "MLM" )
cat( "Selecting:\n" )
print( selected )
cat( "Core Selecting:\n" )
print( core.selected )




# For each scenario, calculate relative performance rates using the design based
# estimators as reference categories
calc.stats = function( runs.m ) {
  m.f.stu = filter( runs.m, method=="DB_indiv (FE)" )$RMSE.f.stu
  m.s.stu = filter( runs.m, method=="DB_indiv (FE)" )$RMSE.s.stu
  m.f.dstu = filter( runs.m, method=="DB_indiv (FE)" )$RMSE.f.dstu
  m.s.dstu = filter( runs.m, method=="DB_indiv (FE)" )$RMSE.s.dstu
  
  m.f.sch = filter( runs.m, method=="DB_clust (FE)" )$RMSE.f.sch
  m.s.sch = filter( runs.m, method=="DB_clust (FE)" )$RMSE.s.sch
  m.f.dsch = filter( runs.m, method=="DB_clust (FE)" )$RMSE.f.dsch
  m.s.dsch = filter( runs.m, method=="DB_clust (FE)" )$RMSE.s.dsch
  
  fbase = filter( runs.m, method=="DB_indiv (FE)" )
  Var.s.base = fbase$SE.s^2 / ( 2 * (fbase$R-1) )     # WHAT IS R???
  
  #mns = filter( mns, method %in% selected )
  #  mns = dplyr::select( runs.m, method, tau, dependence, proptx.dependence, weight, population, 
  #                       SE.f, SE.s, RMSE.fp, RMSE.sp, RMSE.fs, RMSE.ss, bias.sp, bias.ss ) %>%
  mns = mutate( runs.m, 
                
                # proportion of RMSE that is bias, for superpopulation quantities only (?)
                bias.stu.pct = 100 * bias.s.stu^2 / RMSE.s.stu^2,
                bias.sch.pct = 100 * bias.s.sch^2 / RMSE.s.sch^2,
                bias.dstu.pct = 100 * bias.s.dstu^2 / RMSE.s.dstu^2,
                bias.dsch.pct = 100 * bias.s.dsch^2 / RMSE.s.dsch^2,
                
                RMSE.f.stu = 100 * RMSE.f.stu / m.f.stu,
                RMSE.s.stu = 100 * RMSE.s.stu / m.s.stu,
                RMSE.f.sch = 100 * RMSE.f.sch / m.f.sch,
                RMSE.s.sch = 100 * RMSE.s.sch / m.s.sch,
                RMSE.f.dstu = 100 * RMSE.f.dstu / m.f.dstu,
                RMSE.s.dstu = 100 * RMSE.s.dstu / m.s.dstu,
                RMSE.f.dsch = 100 * RMSE.f.dsch / m.f.dsch,
                RMSE.s.dsch = 100 * RMSE.s.dsch / m.s.dsch,
                
                SE.s.inf = 100 * SE.s / SE.f,
                SE.f = 100 * SE.f / fbase$SE.f,
                SE.SE.s = 100 * sqrt( (1/fbase$SE.s^2) * (SE.s^2 / (2*(R-1))) + (SE.s^2/fbase$SE.s^4) * Var.s.base ),
                SE.s = 100 * SE.s / fbase$SE.s
  )
  mns = rename( mns, 
                bias.stu = bias.s.stu,
                bias.sch = bias.s.sch )
  mns
}


# # Change method names to names without reference to superpop vs. finite pop
# refactor_methods = function( dat ) {
#   mutate( dat, method = factor( method, 
#                                 levels=c( "DB-FP-Persons", "DB-FP-Sites", "FE", "RICC", "FIRC", "RIRC" ),
#                                 labels=c( "DB-Persons", "DB-Sites", "FE", "RICC", "FIRC", "RIRC" ) ) )
# }


# For making tables with intervals
#  - computes (0.05, 0.95) quantile interval
get.interval = function( vals, digits=0 ) {
  qt = quantile( vals, c( 0.05, 0.95 ) )
  if ( diff( range( qt, 100 ) ) < 0.00001 ) {   # what is this???
    "Pig"
  } else {
    #sprintf( "%.0f (%.0f-%.0f)", mean( vals ), qt[[1]], qt[[2]] )
    sprintf( paste0( "(%.", digits, "f-%.", digits, "f)" ), qt[[1]], qt[[2]] )
  }
}


#### Calculate performance statistics of the main estimators for all scenarios ####


sruns.core = sruns %>% 
  filter( method %in% selected )
  # refactor_methods()

# runs.m <- sruns.core %>% filter(ID == 1)
mns.full = sruns.core %>%
  group_by( ID ) %>% 
  do( calc.stats(.) ) %>% 
  ungroup()
mns = mns.full %>% 
  dplyr::filter(method %in% core.selected )
# mns$method = droplevels( mns$method )


#### Divide scenarios into variable (main) scenarios and the rest  #####

mns.novar = filter( mns, variable.p == FALSE | variable.n == FALSE )
#mns = filter( mns, variable.p & variable.n )

sruns.novar = filter( sruns, variable.p == FALSE | variable.n == FALSE )
#sruns = filter( sruns, variable.p & variable.n )



####  True SE performance metrics ####


cat( "Now processing SE performance metrics\n" )

table( mns$method )
SE.table = mns %>% 
  filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( finite = mean( SE.f ),
             range.f = get.interval( SE.f ),
             super = mean( SE.s ),
             range.s = get.interval( SE.s ),
             inflation = mean( SE.s.inf ),
             range.i = get.interval( SE.s.inf ) )

SE.table[ SE.table=="Pig" ] = NA
SE.table


# filter( mns, method=="FIRC", SE.f > 115 )
# mns$n = mns$n.bar * mns$J 
# ggplot( filter( mns, method=="FIRC" ), aes( n, SE.f, pch=dependence, col=proptx.dependence ) ) +
#   facet_grid( J ~ tau ) +
#   geom_point()
# 
# 
# filter( mns, method=="FIRC", tau==0, J == 20 )
# M0 = lm( SE.f ~ dependence + proptx.dependence + J + n.bar + n + tau, data=filter( mns, method=="FIRC" ) )
# summary( M0 )
# round( coef( M0 ), digits=1 )


# Looking at simulation uncertainty for the table
head( mns )
summary( mns$SE.SE.s )

mns %>% 
  filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( finite = mean( SE.f ),
             SE = sd( SE.f ) / sqrt( n() ),
             super = mean( SE.s ),
             SE.s = sd( SE.s ) / sqrt( n() ) ) %>%
  mutate_at( c("SE","SE.s"), round, digits=3 )


#### Looking at scenarios with no proportion or site size variation ####


# # Sanity check that no variation means point estimates are all the same.  They are
# simdata.const = filter( simdata, variable.p == FALSE, variable.n == FALSE ) %>%
#           filter( method %in% selected )
# head( simdata.const )
# dbs = simdata.const %>% filter( J == 40, n.bar == 200 ) %>%
#   dplyr::select( ID, subrun, method, ATE.hat ) %>% spread( key=method, value=ATE.hat )
# head( dbs )
# 
# 
# mns.const = sruns.novar %>%   filter( method %in% selected ) %>%
#   refactor_methods() %>%
#   group_by( ID ) %>% do( calc.stats(.) ) %>% ungroup()
# mns.const
# nrow( mns.const )
# 
# table( mns.const$J, mns.const$n.bar )
# table( mns.const$tau )
# table( mns.const$ID )
# length( unique( mns.const$ID ) )
# 
# # Should all be near 0.20
# summary( mns.const$ATE.super.person )
# 
# head( mns.const )
# 
# SE.table.const = mns.const %>% group_by( method ) %>%
#   summarise( finite = mean( SE.f ),
#              range.f = get.interval( SE.f ),
#              super = mean( SE.s ),
#              range.s = get.interval( SE.s ),
#              inflation = mean( SE.s.inf ),
#              range.i = get.interval( SE.s.inf ) )
# 
# SE.table.const[ SE.table.const=="Pig" ] = NA
# SE.table.const
#   
# if ( FALSE ) {
#   tt = filter( mns, method == "DB-Sites" )
#   nrow( tt )
#   head( tt )
#   aa = dplyr::select( tt, method, tau, dependence, proptx.dependence )
#   aa = unique( aa )
#   aa
#   nrow( aa )
# }



#### RMSE ####
cat( "Now processing RMSE performance metrics\n" )

RMSE.table = mns %>% 
  filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( mn.f.stu = mean( RMSE.f.stu ),
             rng.f.stu = get.interval( RMSE.f.stu ),
             mn.s.stu = mean( RMSE.s.stu ),
             rng.s.stu = get.interval( RMSE.s.stu ),
             mn.f.dstu = mean( RMSE.f.dstu ),
             rng.f.dstu = get.interval( RMSE.f.dstu ),
             mn.s.dstu = mean( RMSE.s.dstu ),
             rng.s.dstu = get.interval( RMSE.s.dstu ),
             
             mn.f.sch = mean( RMSE.f.sch ),
             rng.f.sch = get.interval( RMSE.f.sch ),
             mn.s.sch = mean( RMSE.s.sch ),
             rng.s.sch = get.interval( RMSE.s.sch ),
             mn.f.dsch = mean( RMSE.f.dsch ),
             rng.f.dsch = get.interval( RMSE.f.dsch ),
             mn.s.dsch = mean( RMSE.s.dsch ),
             rng.s.dsch = get.interval( RMSE.s.dsch ) )

RMSE.table[ RMSE.table == "Pig" ] = NA
RMSE.table
#for ( i in 2:ncol(mns.agg) ) {
#  mns.agg[ mns.agg[[i]] < 0.0001, i ] = NA
#}

# if (FALSE) {
# # Looking at ranges
# mns2 = gather(mns, RMSE.fp, RMSE.sp, RMSE.fs, RMSE.ss, key="estimand", value="RMSE" )
# head( mns2 )
# mns2 = mns2 %>% group_by( method,estimand ) %>% 
#   summarise( q20 = quantile( RMSE, 0.20 ),
#              q80 = quantile( RMSE, 0.80 ) )
# head( mns2 )
# filter( mns2, method=="FIRC" )
# }


#names(mns.agg) = gsub( "mean.", "", names(mns.agg ) )
#for ( i in c(8,6,4,2 ) ) {
#  mns.agg[,i] = paste0( "$", round( mns.agg[,i,drop=TRUE] ), " \\pm ", round( 2 * mns.agg[,i+1, drop=TRUE] ), "$"  )
#}
#mns.agg = mns.agg[ c(1,2,4,6,8) ]
#mns.agg[ mns.agg == "$100 \\pm NA$" ] = "100"


# SE for the mean RMSE estimates
mns %>% 
  filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( mn.f.stu = mean( RMSE.f.stu ),
             SE.f.stu = sd( RMSE.f.stu ) / sqrt(n()),
             mn.s.stu = mean( RMSE.s.stu ),
             SE.s.stu = sd( RMSE.s.stu ) / sqrt(n()),
             mn.f.dstu = mean( RMSE.f.dstu ),
             SE.f.dstu = sd( RMSE.f.dstu ) / sqrt(n()),
             mn.s.dstu = mean( RMSE.s.dstu ),
             SE.s.dstu = sd( RMSE.s.dstu ) / sqrt(n()),
             
             mn.f.sch = mean( RMSE.f.sch ),
             SE.f.sch = sd( RMSE.f.sch ) / sqrt(n()),
             mn.s.sch = mean( RMSE.s.sch ),
             SE.s.sch = sd( RMSE.s.sch ) / sqrt(n()),
             mn.f.dsch = mean( RMSE.f.dsch ),
             SE.f.dsch = sd( RMSE.f.dsch ) / sqrt(n()),
             mn.s.dsch = mean( RMSE.s.dsch ),
             SE.s.dsch = sd( RMSE.s.dsch ) / sqrt(n()) ) %>%
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

bias.table = mns %>% 
  filter( method %in% core.selected ) %>%
  group_by( method ) %>%
  summarise( mean.bias.stu = mean( bias.stu ),
             rng.bias.stu = get.interval.bias( bias.stu ),
             sd.bias.stu = sd( bias.stu ),
             mean.bias.sch= mean( bias.sch ),
             rng.bias.sch = get.interval.bias( bias.sch ),
             sd.bias.sch = sd( bias.sch ) )

bias.table$sd.bias.stu = bias.table$sd.bias.sch = NULL
bias.table


summary( mns$bias.stu.pct )

bias.table.pct = mns %>% 
  group_by( method ) %>%
  summarise( mean.bias.stu = mean( bias.stu.pct ),
             rng.bias.stu = get.interval( bias.stu.pct ),
             mean.bias.sch = mean( bias.sch.pct ),
             rng.bias.sch = get.interval( bias.sch.pct ) )
bias.table.pct



# #### Examining RICC with respect to FE   and RIRC with respect to FIRC ####
# 
# head( simdata )
# table( simdata$method )
# csim = filter( simdata, method %in% c( "RICC", "FE", "RIRC", "FIRC" ), variable.n == TRUE, variable.p == TRUE )
# 
# head( csim )
# csim$SE.hat = NULL
# csim = spread( csim, method, ATE.hat )
# head( csim )
# 
# cor( csim$FE, csim$RICC )
# round( summary( csim$FE - csim$RICC ), digits= 3 )
# mean( abs( csim$FE - csim$RICC ) > 0.01 )
# 
# cor( csim$FE, csim$FIRC )
# round( summary( csim$FE - csim$FIRC ), digits= 3 )
# 
# cor( csim$RIRC, csim$FIRC )
# round( summary( csim$RIRC - csim$FIRC ), digits= 3 )
# sd( csim$RIRC - csim$FIRC )
# quantile( csim$RIRC - csim$FIRC, c( 0.05, 0.95 ) )
# mean( abs( csim$RIRC - csim$FIRC ) >= 0.01 )
# mean( abs( csim$RIRC - csim$FIRC ) >= 0.02 )
# 
# # How do RMSEs compare?
# cmns = filter( mns.full, method %in% c( "RICC", "FE", "RIRC", "FIRC" ) ) %>%
#   dplyr::select( method, ID, RMSE.fp ) %>%
#   spread( method, RMSE.fp )
# head( cmns )
# nrow( cmns )
# summary( cmns$RICC / cmns$FE )
# t.test( cmns$RICC / cmns$FE )
# qplot( cmns$RICC / cmns$FE )
# 
# nrow( cmns )
# summary( cmns$RIRC / cmns$FIRC )
# t.test( cmns$RIRC / cmns$FIRC )
# qplot( cmns$RIRC / cmns$FIRC )



#### Sanity checking: RMSE^2 = bias^2 + SE^2 ####

head( sruns )
ss = sruns %>% 
  mutate( tt.f = bias.f.stu^2 + SE.f^2,
          tt.s = bias.s.stu^2 + SE.s^2, 
          tt.s.sch = bias.s.sch^2 + SE.s^2,
          
          RMSE.check.f = tt.f - RMSE.f.stu^2,
          RMSE.check.s = tt.s - RMSE.s.stu^2,
          RMSE.check.s.sch = tt.s.sch / RMSE.s.sch^2,
          birat = bias.s.stu^2 / RMSE.s.stu^2 )
ss$RMSE.check.f %>% summary()
ss$RMSE.check.s %>% summary()
ss$RMSE.check.s.sch %>% summary()

# summary( ss$RMSE.check.s.sch )
# ss %>% 
#   filter( RMSE.check.s.sch > 8 ) %>% 
#   dplyr::select( method, ID, bias.s.sch, RMSE.s.sch, SE.s )
# 
# if (FALSE) {
#   sss = filter( simdata, ID == 1041, method=="FE" )
#   nrow( sss )
#   SE2.true = var( sss$ATE.hat )
#   unique( sss$ATE.super.site )
#   bias.true = mean( sss$ATE.hat ) - 0.2  
#   bias.true
#   RMSE.true = 4
# }

wack = filter( ss, birat > 1 )
nrow( wack )
summary( wack$birat )
head( wack )

round( summary( ss$tt.f / ss$RMSE.fp^2 ), digits=3 )

round( summary( ss$tt.s / ss$RMSE.fs^2 ), digits=3 )

# Check results: Close, but not perfect.



# Looking at SE vs RMSE
head( sruns )

summary( 100 * sruns$SE.f^2 / sruns$RMSE.f.stu^2)

summary( 100 * sruns$SE.s^2 / sruns$RMSE.s.stu^2)






# ##### Comparing the superpop SE.hat to finite pop SE.hat #####
# 
# cat( "Now processing SE.hat comparison metrics\n" )
# 
# # How often are superpop SEs lower than finite pop SEs?
# # Comparing FIRC to DB
# head( simdata )
# unique( simdata$method )
# ssub = filter( simdata, method %in% c( "FIRC", "DB-FP-Persons", "DB-FP-Sites", "FE-CR" ) ) %>% 
#   dplyr::select( ID, J, n.bar, tau, dependence, subrun, method, SE.hat )
# head( ssub )
# ssub = mutate( ssub, subrun = paste( ID, subrun, sep="." ) )
# ssub = spread( ssub, key=method, value=SE.hat )
# 
# head( ssub )
# nrow( ssub )
# # Are SE estimates correlated?
# 
# # This divided by dependance---but that seems unneeded.
# #ssub = ssub %>% group_by( J, n.bar, tau, dependence ) %>% 
# #  summarise( FIRC.less = mean( FIRC < `DB-FP-Persons` ),
# #             FE.CR.less = mean( `FE-CR` < `DB-FP-Persons` ),
# #             FIRC.way.less = mean( FIRC < `DB-FP-Persons` - 0.025 ),
# #             FE.CR.way.less = mean( `FE-CR` < `DB-FP-Persons` - 0.025 ),
# #             )
# #ggplot( ssub, aes( tau, percent, col=method, lty=dependence ) ) +
# #  facet_grid( J ~ n.bar, labeller = label_both ) +
# #  geom_point() + geom_line()
# ssub = ssub %>% 
#   group_by( J, n.bar, tau ) %>% 
#   summarise( FIRC.less = mean( FIRC < `DB-FP-Persons` ),
#              FE_CR.less = mean( `FE-CR` < `DB-FP-Persons` ),
#              FIRC.way_less = mean( FIRC < `DB-FP-Persons`*0.90 ),
#              FE_CR.way_less = mean( `FE-CR` < `DB-FP-Persons`*0.90 ),
#              FIRC.less_site = mean( FIRC < `DB-FP-Sites` ),
#              FIRC.way_less_site = mean( FIRC < `DB-FP-Sites`*0.90 )   )
# 
# head( ssub )
# ssub2 = gather( ssub, FIRC.less, FE_CR.less, FIRC.way_less, FE_CR.way_less, FIRC.less_site, FIRC.way_less_site, key="method", value="percent" )
# head( ssub2 )
# table( ssub2$method )
# ssub2 = separate( ssub2, method, into=c("method","threshold"), sep="\\." )
# table(ssub2$method )
# table(ssub2$threshold )
# ssub2
# SE.comp.results = ssub2



# ##### Looking at correlation of the SE hats ######
# 
# cat( "Now processing correlation of SE hats for specific scenario\n" )
# 
# head( simdata )
# unique( simdata$method )
# ssub = simdata %>% filter( J %in% c( 20, 40 ), n.bar == 200, dependence == "indep", proptx.dependence == "indep" ) %>%
#   filter( method %in% c( "FIRC", "DB-FP-Persons", "FE-CR", "FE-Club", "DB-SP-Sites" ) ) %>% 
#   dplyr::select( ID, J, n.bar, tau, dependence, proptx.dependence, run, subrun, method, SE.hat )
# head( ssub )
# ssub = spread( ssub, key=method, value=SE.hat )
# head( ssub )
# ssub = separate(ssub, subrun, into=c("batch", "rr", "scount" ) )
# head( ssub )
# 
# #s2 = filter( ssub, J  %in% c(20,80), n.bar == 100, tau %in% c(0, 0.2), dependence == "corr", proptx.dependence == 0 )
# s2 = filter( ssub, scount == 1)
# table( s2$ID )
# head( s2 )
# nrow( s2 )
# s2 = mutate( s2, smaller.FIRC = ifelse( 0.99*`DB-FP-Persons` > `FIRC`, "flip", "ok" ),
#              smaller.CR = ifelse( 0.99*`DB-FP-Persons` > `FE-CR`, "flip", "ok" ),
#              small = paste0( smaller.CR, "-", smaller.FIRC ) )
# 
# SE.corr.results = s2
# head( SE.corr.results )
# #rm( simdata, s2, ssub2, ssub )
# #rm( scenarios )

cat( "Saving all the processed results\n" )

save( bias.table, bias.table.pct, 
      RMSE.table, SE.table, 
      # SE.comp.results, SE.corr.results, 
      mns, sruns, sruns.core,
      get.interval, get.interval.bias, file = "results/aggregated_simulation_results.RData" )




