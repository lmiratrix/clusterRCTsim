# ---
# title: "Focused Comparison of DB Site to FIRC"
# author: "Miratrix, Weiss, Henderson"
# date: "`r Sys.Date()`"
# ---

# This script looks at DB Site vs FIRC and assess how they perform relative to
# each other.
#
# Written in response to reviewer queries.


#### Load data and libraries ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library( tidyverse )
load( "../results/aggregated_simulation_results.RData" )

library( ggthemes )
my_theme = theme_tufte() + theme( legend.position="bottom", 
                                  legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                  legend.box = "vertical", 
                                  panel.border=element_rect(size=0.3, fill=NA))
theme_set( my_theme )


# Get just DB-Sites and FIRC
sruns.core = filter( sruns.core, method %in% c(" DB-Sites", "FIRC" ) )
nrow( sruns.core )

mns = filter( mns, method %in% c( "DB-Sites", "FIRC" ) )
levels( mns$method )
table( mns$method )


mns
table( mns$dependence )
table( mns$dependence, mns$proptx.dependence )
table( mns$dependence, mns$proptx.dependence, mns$tau )



##### Look at true Superpopulation SE  #####

m2 = mns %>% dplyr::select( method, ID, J, n.bar, 
                            dependence, proptx.dependence, tau, SE.s, RMSE.ss ) %>%
  pivot_wider( names_from=method, values_from=c( SE.s, RMSE.ss ) )
head( m2 )

m2 = mutate( m2, ratio_SE = `SE.s_DB-Sites` / SE.s_FIRC,
             ratio_RMSE = `RMSE.ss_DB-Sites` / RMSE.ss_FIRC )
qplot( m2$ratio_SE )

ggplot( m2, aes( tau, ratio_SE, col=J ) ) +
  facet_grid( dependence ~ . ) +
  geom_point() +
  geom_hline( yintercept = 1 ) 

quantile( m2$ratio_SE, c(0.05, 0.5, 0.95 ) )
mean( m2$ratio_SE )


####  How do RMSEs compare?  #####
ggplot( m2, aes( tau, ratio_RMSE, col=J ) ) +
  facet_grid( dependence ~ proptx.dependence ) +
  geom_jitter( width=0.01) +
  geom_hline( yintercept = 1 ) 


quantile( m2$ratio_RMSE, c(0.05, 0.5, 0.95 ) )
mean( m2$ratio_RMSE )

m2
ggplot( m2, #filter( m2 , dependence != "const" ), 
        aes( `ratio_SE`, `ratio_RMSE`, col=tau ) ) +
  facet_grid( dependence ~ J ) +
  geom_point( alpha=0.75 ) +
#  geom_hline( yintercept = 1 ) +
#  geom_vline( xintercept = 1 ) +
  labs( y = "Relative Superpopulation Site RMSE of DB-Site / FIRC",
        x = "Relative Superpopulation variability in estimating SE of DB-Site / FIRC") +
#  coord_fixed() +
  theme( legend.position="right", 
         legend.direction="vertical", legend.key.width=unit(1,"cm"),
         legend.box = "vertical", 
         panel.border=element_rect(size=0.3, fill=NA))





##### When does DB win?  #####

# We know for constant site size point estimates are identical, so we drop
# those.
m2_noconst = filter( m2, dependence != "const" )

names( m2_noconst )
#m2_noconst$dependence = relevel( m2_noconst$dependence, ref = "const" )
m2_noconst$proptx.dependence = relevel( m2_noconst$proptx.dependence, ref = "const" )
M0 = lm( log(ratio_RMSE) ~ (as.factor(tau) + as.factor(n.bar) + as.factor(J) + dependence + proptx.dependence)^2, data=m2_noconst )
summary( M0 )

m2_noconst$predRat = exp( predict( M0 ) )
mean( m2_noconst$predRat < 0.99 )
mean( m2_noconst$predRat < 1 )
nrow( m2_noconst )
ff =  m2_noconst %>% filter( predRat < 1 ) %>%
  select( -`SE.s_DB-Sites`, -SE.s_FIRC, -ID )
ff
mean( ff$predRat )
table( m2_noconst$J )

ggplot( m2_noconst, aes( ratio_RMSE, predRat ) ) +
  geom_point() +
  geom_abline( slope=1, intercept = 0 ) +
  geom_hline( yintercept = 1 )

ggplot( m2_noconst, aes( n.bar, predRat, col=tau ) ) +
  facet_grid( dependence ~ J ) +
  geom_point( alpha=0.75 ) +
  geom_hline( yintercept = 1 ) +
  theme( legend.position="right", 
         legend.direction="vertical", legend.key.width=unit(1,"cm"),
         legend.box = "vertical", 
         panel.border=element_rect(size=0.3, fill=NA))


# Percent of time DB wins
mean( m2_noconst$ratio_RMSE < 0.99 )
mean( m2_noconst$ratio_RMSE < 0.97 )
min( m2_noconst$ratio_RMSE )

mgd = filter( m2_noconst, ratio_RMSE < 0.99 )
table( mgd$J ) 
table( mgd$tau ) 
#filter( mgd, tau==0 ) %>% select( -`SE.s_DB-Sites`, -SE.s_FIRC )




#### When site sizes are identical (no site variation) ####

# If proportion vary in site, then we still have differences, but these are small

# Ratio of true SEs  SE[DB] / SE[FIRC]
m3 = filter( m2, dependence == "const" )
head( m3 )
summary( m3$ratio_SE )
# so DB is _slightly_ higher.

# Check: equal prop tx then all are same
m4 = filter( m2, dependence == "const", proptx.dependence=="const" )
head( m4 )
summary( m4$ratio_SE )


m3 = filter( m2, dependence == "const" )
head( m3 )
summary( m3$ratio_RMSE )




##### Standard Error Estimation Comparison #####

srun2 = dplyr::select( sruns, method, weight, population, 
                       ID, J, n.bar, dependence, proptx.dependence, tau, 
                       sd.SE.hat.s, sd.SE.hat.f, SE.s, SE.f ) %>%
  mutate( sd.SE.hat.s = 100 * sd.SE.hat.s / SE.s,
          sd.SE.hat.f = 100 * sd.SE.hat.f / SE.f )
srun2 = mutate( srun2, sd.SE = ifelse( population=="finite", sd.SE.hat.f, sd.SE.hat.s ) )


m2_SE = srun2 %>% dplyr::select( method, ID, J, n.bar, dependence, proptx.dependence, 
                                 tau, sd.SE.hat.s ) %>%
  spread( method, sd.SE.hat.s )
head( m2_SE )
m2_SE = mutate( m2_SE, ratio_SE_hat = `DB-SP-Sites` / FIRC )
ggplot( m2_SE, aes( tau, ratio_SE_hat ) ) +
  facet_grid( dependence ~ proptx.dependence ) +
  geom_point() +
  geom_hline( yintercept = 1 )

quantile( m2_SE$ratio_SE_hat, c(0.05, 0.5, 0.95 ) )
mean( m2_SE$ratio_SE_hat )



##### SE estimation vs. RMSE ######

all( m2$ID  == m2_SE$ID )
m2$ratio_SE_hat = m2_SE$ratio_SE_hat

ggplot( m2, #filter( m2 , dependence != "const" ), 
        aes( ratio_SE_hat, ratio_RMSE, col=tau ) ) +
  facet_grid( dependence ~ J ) +
  geom_point( alpha=0.75 ) +
  geom_hline( yintercept = 1 ) +
  geom_vline( xintercept = 1 ) +
  labs( y = "Relative superpop site RMSE\nof DB-Site vs FIRC",
        x = "Relative variability in estimating SE of DB-Site vs FIRC") +
  coord_fixed() +
  theme( legend.position="right", 
                                    legend.direction="vertical", legend.key.width=unit(1,"cm"),
                                    legend.box = "vertical", 
                                    panel.border=element_rect(size=0.3, fill=NA))

ggsave( filename="../figures/compare_FIRC_v_DB.pdf", width=8, height=5)








