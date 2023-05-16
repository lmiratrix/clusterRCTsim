# ---
# title: "Calculate paper results"
# author: "Miratrix"
# date: "`r Sys.Date()`"
# ---

# This script generates the core figures and tables (pretaining to the
# simulations) found in the paper and supplements.
#
# (C) 2020 Miratrix
# 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library( tidyverse )
options(list(dplyr.summarise.inform = FALSE))

# Load the prepared data.
# See 'data_prep' for how this file is generated from the raw results.
load( "../results/aggregated_simulation_results.RData" )


library( ggthemes )
my_theme = theme_tufte() + theme( legend.position="bottom", 
                                  legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                  legend.box = "vertical", 
                                  panel.border=element_rect(size=0.3, fill=NA))
theme_set( my_theme )



#### Relative SE plot ####

mns = mutate( mns, 
              biased = ifelse( method %in% c("FIRC","FE"), "biased", "DB" ) )


mns2 = mns %>% filter( proptx.dependence == "indep", dependence=="indep" ) %>% 
  gather( SE.f, SE.s, key="estimand", value="SE" )

s2 = mns2 %>% ungroup( ) %>% 
  group_by( method, weight, population, tau, proptx.dependence, dependence, estimand ) %>% 
  summarise( mean.SE = mean( SE ) )

table( s2$estimand )

SEnames = c( `SE.f`="SE (finite sample)", `SE.s`="SE (superpopulation)" )
ggplot( s2, aes( x=tau, y=mean.SE, col=method, fill=method, pch = method, lty=method ) ) +
  facet_grid( . ~ estimand, labeller = as_labeller(SEnames)  ) +
  geom_line() + geom_point(size=2) +
  labs( y = "Avg Percent Increase of SE", x = expression( tau ) ) +
  scale_shape_manual("", values=c(`DB-Persons`=21,`DB-Sites`=22,FE=23,FIRC=24)) +
  scale_fill_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_color_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_linetype_manual("", values = c(`DB-Persons` = 1, `DB-Sites`=1, FE=2, FIRC=2) ) 




## ----RMSE table, echo=FALSE------------------------------
knitr::kable( RMSE.table, digits=0 )



## ----RMSE_finite_person, echo=FALSE, results="hide", fig.height=5.5----
s2 = mns %>% 
  #filter( proptx.dependence == "indep" ) %>%
  ungroup( ) %>% group_by( method, weight, population, tau, proptx.dependence, dependence, biased ) %>% 
  summarise( mean.RMSE.fp = mean( RMSE.fp ),
             mean.RMSE.fs = mean( RMSE.fs ),
             mean.RMSE.sp = mean( RMSE.sp ),
             mean.RMSE.ss = mean( RMSE.ss ) )
s2 = gather( s2, mean.RMSE.fp, mean.RMSE.fs, mean.RMSE.sp, mean.RMSE.ss,
             key="weighting", value="mean.RMSE" )

s2 = mutate( s2, weighting = fct_recode(weighting, 
                                        "person"="mean.RMSE.fp",
                                        "site"="mean.RMSE.fs",
                                        "sup person" = "mean.RMSE.sp",
                                        "sup site" = "mean.RMSE.ss" ))
s2 = rename( s2, ptx.dep = proptx.dependence,
             size.dep = dependence )



s3 = filter( s2,  weighting=="person",
             method !="DB-Sites")

ptx.names = c( ncorr="Neg Corr p", const="Constant p", indep="Indep p", pcorr="Pos Corr p")
size.names = c( corr="Corr N", const="Constant N", indep="Indep N" )

ggplot( s3, aes( tau, mean.RMSE, col=method, fill=method, pch = method, lty=method  ) ) +
  facet_grid( ptx.dep ~ size.dep, labeller = labeller( ptx.dep=as_labeller(ptx.names), size.dep=as_labeller(size.names) ) ) +
  geom_hline( yintercept = c(90,100, 110), size=0.5, col="grey") +
  geom_point(size=3) + geom_line( ) +
  labs( title= "Average Relative RMSEs for person-weighted, finite population estimands",
        y="Percent increase in RMSE", x=expression(tau)) +
  coord_cartesian( ylim=c(80,120) ) + 
  scale_shape_manual("", values=c(`DB-Persons`=21,`DB-Sites`=22,FE=23,FIRC=24)) +
  scale_fill_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_color_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_linetype_manual("", values = c(`DB-Persons` = 1, `DB-Sites`=1, FE=2, FIRC=2) ) 




## ----RMSE_super_person, echo=FALSE, results="hide", fig.height=5.5----

s3 = filter( s2,  weighting=="sup person",
             method !="DB-Sites")

ggplot( s3, aes( tau, mean.RMSE, lty=method, col=method, fill=method, pch=method  ) ) +
  facet_grid( ptx.dep ~ size.dep, labeller = labeller( ptx.dep=as_labeller(ptx.names), size.dep=as_labeller(size.names) )  ) +
  geom_hline( yintercept = c(90,100, 110), size=0.5, col="grey") +
  geom_point(size=3) + geom_line( ) +
  labs( title= "Average Relative RMSEs for person-weighted, superpopulation estimands",
        y="Percent increase in RMSE", x=expression(tau) ) +
  coord_cartesian( ylim=c(80,120) ) + 
  scale_shape_manual("", values=c(`DB-Persons`=21,`DB-Sites`=22,FE=23,FIRC=24)) +
  scale_fill_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_color_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_linetype_manual("", values = c(`DB-Persons` = 1, `DB-Sites`=1, FE=2, FIRC=2) ) 



## ----RMSE_finite_site, echo=FALSE, results="hide", fig.height=5.5----

s4 = filter( s2, weighting=="site" )
ggplot( s4, aes( tau, mean.RMSE,  lty=method, col=method, fill=method, pch=method  ) ) +
  facet_grid( ptx.dep ~ size.dep, labeller = labeller( ptx.dep=as_labeller(ptx.names), size.dep=as_labeller(size.names) ) ) +
  geom_hline( yintercept = c(90,100, 110), size=0.5, col="grey") +
  geom_point(size=2) + geom_line( ) +
  labs( title= "Average Relative RMSEs for site-weighted, finite population estimands",
        y="Percent increase in RMSE", x=expression(tau) ) + 
  scale_shape_manual("", values=c(`DB-Persons`=21,`DB-Sites`=22,FE=23,FIRC=24)) +
  scale_fill_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_color_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_linetype_manual("", values = c(`DB-Persons` = 1, `DB-Sites`=1, FE=2, FIRC=2) ) 




## ----RMSE_super_site, echo=FALSE, results="hide", fig.height=5.5----

s4 = filter( s2,  weighting=="sup site" )

ggplot( s4, aes( tau, mean.RMSE, lty=method, col=method, fill=method, pch=method  ) ) +
  facet_grid( ptx.dep ~ size.dep, labeller = labeller( ptx.dep=as_labeller(ptx.names), size.dep=as_labeller(size.names) ) ) +
  geom_hline( yintercept = c(90,100, 110), size=0.5, col="grey") +
  geom_point(size=3) + geom_line( ) +
  labs( title= "Average Relative RMSEs for site-weighted, superpopulation estimands",
        y="Percent increase in RMSE", x=expression(tau)) + 
  scale_shape_manual("", values=c(`DB-Persons`=21,`DB-Sites`=22,FE=23,FIRC=24)) +
  scale_fill_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_color_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_linetype_manual("", values = c(`DB-Persons` = 1, `DB-Sites`=1, FE=2, FIRC=2) ) 


ss = s4 %>% filter( weighting =="sup site" ) %>% 
  ungroup() %>% dplyr::select( ptx.dep, size.dep, tau, method, mean.RMSE ) %>%
  spread( method, mean.RMSE )
head( ss )
ss = mutate( ss, ratio = FIRC / `DB-Sites` )
summary( ss$ratio )



#### Bias Exploration ####

## Bias Table ---- echo=FALSE-----------------------------------------
knitr::kable( bias.table, digits=2 )


mns2 = gather( mns, bias.fp, bias.fs, key="estimand", value="bias" )
head( mns2 )
s2 = mns2 %>% ungroup( ) %>% 
  group_by( method, weight, population, tau, proptx.dependence, dependence, estimand ) %>% 
  summarise( mean.bias = mean( bias ) )
s2 = mutate( s2, estimand = fct_recode( estimand, "person" = "bias.fp", "site"="bias.fs" ) )
#mns2 = filter( mns2, dependence =="corr" )
s2 = rename( s2, ptx.dep = proptx.dependence,
             size.dep = dependence )


## ----mean_bias, echo=FALSE, results="hide", fig.height=4.5, fig.width = 5----
s3 = filter( s2,  estimand=="person", size.dep != "const", ptx.dep != "const" )


ggplot( s3, aes( x=tau, y=mean.bias, col=method, fill=method, pch = method, lty=method ) ) +
  facet_grid( ptx.dep ~ size.dep, labeller = labeller( ptx.dep=as_labeller(ptx.names), size.dep=as_labeller(size.names) ) ) +
  geom_line() + geom_point() +
  geom_hline( yintercept = 0 )+ 
  scale_shape_manual("", values=c(`DB-Persons`=21,`DB-Sites`=22,FE=23,FIRC=24)) +
  scale_fill_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_color_manual("", values = c(`DB-Persons` = "red", FE="red", `DB-Sites`="blue", FIRC = "blue")) +
  scale_linetype_manual("", values = c(`DB-Persons` = 1, `DB-Sites`=1, FE=2, FIRC=2) ) +
  labs( y="Mean Bias", x=expression(tau) )




## ---- echo=FALSE-----------------------------------------
names(bias.table.pct) = c( "Method", "Mean Percent Bias (Person)", "90% range", "Mean Percent Bias (Site)", "90% Range" )
knitr::kable( bias.table.pct, digits=0 )


#### Estimated Standard Errors ####

## ----calibratedSEs, echo=FALSE, results="hide"-----------
inf.names = c( inflate.f = "Inflation (Finite)", inflate.s = "Inflation (Superpop)" )
srun2 = dplyr::select( sruns, method, population, inflate.f, inflate.s, tau, population, J, n.bar ) %>%
  gather( inflate.f, inflate.s, key="estimand", value="inflation" )

ggplot( srun2, aes( x=method, y=inflation ) ) +
  facet_grid( ~ estimand, scales = "free", labeller = as_labeller(inf.names) ) +
  geom_boxplot() +
  coord_flip() + 
  geom_hline( yintercept = 1, lwd=1,col="red" ) +
  labs( y = "Ratios of average SE.hat to true SE across scenarios", x="")



## ----calibratedSEs3, echo= FALSE, results="hide"---------
ggplot( filter( srun2, population == "superpop", estimand=="inflate.s"), aes( x = as.factor(tau), y = inflation ) ) +
  facet_wrap( ~ method, nrow=1 ) +
  geom_boxplot( fill="grey" ) + geom_hline( yintercept = 1 ) +
  labs( title="Estimated to true SE for superpopulation targeting estimators", x=expression(tau), y="Relative Inflation" )



## ----calibratedSEs4, echo= FALSE, results="hide", include=FALSE----
ggplot( filter( srun2, population == "finite", estimand=="inflate.f"), aes( x = as.factor(J), y = inflation ) ) +
  facet_wrap( ~ method, nrow=1 ) +
  geom_boxplot( fill="grey" ) + geom_hline( yintercept = 1 ) +
  labs( title="Estimated to true SE for superpopulation targeting estimators", x="J (number sites)", y="Relative Inflation" )




M0 = lm( inflation ~ method + tau + J + n.bar, data=srun2 )
summary( M0 )



#### Are the estimated standard errors diﬀerently precise? #####

srun2 = dplyr::select( sruns, method, weight, population, 
                       ID, J, n.bar, dependence, proptx.dependence, tau, 
                       sd.SE.hat.s, sd.SE.hat.f, SE.s, SE.f ) %>%
  mutate( sd.SE.hat.s = 100 * sd.SE.hat.s / SE.s,
          sd.SE.hat.f = 100 * sd.SE.hat.f / SE.f )
srun2 = mutate( srun2, sd.SE = ifelse( population=="finite", sd.SE.hat.f, sd.SE.hat.s ) )

ggplot( srun2, aes( x=method, y=sd.SE, fill=population ) ) +
  #  facet_grid( ~ population ) +
  geom_boxplot() +
  coord_flip()+ 
  labs( y = "100*Ratio of Standard deviation of the standard errors to actual SEs", x="" )



# Filter to mild cross site variation
# (To see how FIRC compares to DB-SP-Persons, reviewer request.)
head( srun2 )
sr2 = filter( srun2, method %in% c( "DB-SP-Persons", "FIRC", "DB-SP-Sites", "FE-CR" ) )
ggplot( sr2, aes( x=method, y=sd.SE, fill=population ) ) +
  facet_grid( ~ tau ) +
  geom_boxplot() +
  coord_flip()+ 
  labs( y = "100*Ratio of Standard deviation of the standard errors to actual SEs", x="" )


## ---- echo=FALSE-----------------------------------------
summ<- srun2 %>% group_by( method, weight, population ) %>% 
  summarise( mean.sd = mean(  sd.SE ),
             range.sd = get.interval( sd.SE, digits = 2 ) ) %>%
  arrange( population, weight )
knitr::kable(summ, digits=1 )


## ----SP_Site_v_FIRC_SE, echo=FALSE, results="hide", include=FALSE----
head( srun2 )
m2 = srun2 %>% dplyr::select( method, ID, J, n.bar, dependence, proptx.dependence, tau, sd.SE.hat.s ) %>%
  spread( method, sd.SE.hat.s )
head( m2 )
m2 = mutate( m2, ratio = `DB-SP-Sites` / FIRC )
ggplot( m2, aes( tau, ratio ) ) +
  facet_grid( dependence ~ proptx.dependence ) +
  geom_point() +
  geom_hline( yintercept = 1 )
quantile( m2$ratio, c(0.05, 0.5, 0.95 ) )
mean( m2$ratio )



#### How often are superpop SEs estimated less than finite? #####

rss = SE.comp.results %>% ungroup() %>%
  filter( threshold %in% c("less", "way_less" ) )
ggplot( rss,
        aes( tau, I(100*percent), col=method, lty=threshold ) ) +
  facet_grid( J ~ n.bar, labeller = label_both ) +
  geom_point() + geom_line() +
  coord_cartesian( ylim=c(0,100) ) +
  geom_hline( yintercept=c(0,0.25), col="grey", size=0.3 )  +
  labs( x=expression(tau), y="Percent Chance") +
  scale_x_continuous(breaks=c(0,0.1, 0.2) ) +
  guides( lty=FALSE )

SE.comp.results %>% ungroup() %>%
  filter( tau == 0.1 ) %>% group_by( method, J, n.bar ) %>%
  summarise( per.less = mean( percent ) ) %>%
  spread( method, per.less )





#####  Calculate inflation of DB-Site over FIRC true SEs    ######

unique( sruns.core$method )
ss = filter( sruns.core, method %in% c( "FIRC", "DB-Sites" ) )
head( ss )

ss = dplyr::select( ss, c( method, ID, SE.s, sd.SE.hat.s ) )

ss = mutate( ss, frac.unc = sd.SE.hat.s / SE.s )

# Table of averages across scenarios
ss %>% group_by( method ) %>% 
  summarise( mean.SE.s = mean( SE.s ),
             mean.sd.SE.hat.s = mean( sd.SE.hat.s ),
             frac.unc = mean( frac.unc ))

head( ss )
ss2 = ss %>% dplyr::select( ID, method, sd.SE.hat.s ) %>%
  spread( method, sd.SE.hat.s )
head( ss2 )
summary( 100 * ss2$`DB-Sites` / ss2$FIRC )


##### Look at ratio of true SE of site vs FIRC  ######
ss = dplyr::select( ss, c( method, ID, SE.s ) )
ss2 = spread( ss, method, SE.s )
head( ss2 )
ss = spread( ss, method, SE.s )
head( ss )
ss = mutate( ss, inflate = 100 * `DB-Sites` / FIRC )
summary( ss$inflate )


head( ss )



##### Look at how scenario true ATEs vary #####

ss.scen = filter( sruns.core, method=="DB-Persons" )
dim( ss.scen )
head( ss.scen )

ss.scen = mutate( ss.scen, del.ATE = ATE.super.person - ATE.super.site )
qplot( ss.scen$del.ATE )

qplot( ss.scen$tau, ss.scen$del.ATE ) +
  labs( x = "Cross site impact heterogeniety", 
        y = "Differnence between Site and Person" )





#### Some regressions to see how simulation factors predict true standard error #####

head( sruns.core )

sruns.core$proptx.dependence = factor( sruns.core$proptx.dependence, levels=c(0, -1, 1), labels=c("indep","neg","pos") )
sruns.core$dependence = factor( sruns.core$dependence, levels=c("indep","corr" ) )
sruns.core$tau.f = as.factor( sruns.core$tau )

table( sruns$proptx.dependence, sruns$dependence )
str( sruns )
sruns = mutate( sruns,
                log_SE.s = log( SE.s ),
                J.f = factor( J ),
                n.bar.f = factor( n.bar ),
                tau.f = factor( tau ) )
M0 = lm( log_SE.s ~ method + J.f + n.bar.f + dependence + proptx.dependence + tau.f, 
          data=sruns )

summary( M0 )
exp( coef( M0 ) ) 



M0 = lm(  log(RMSE.fp) ~ method + as.factor(J) + as.factor(n.bar) + dependence + proptx.dependence + tau.f, data=sruns )

summary( M0 )
exp( coef( M0 ) ) 


