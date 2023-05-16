# ---
# title: "Explore a single simulation scenario"
# author: "Miratrix, Weiss, Henderson"
# date: "`r Sys.Date()`"
# ---

# Exploration script to evaluate results from a single scenario of the SE
# assessment simulation.
#
# In particular we can see how estimates of standard errors are correlated
# across datasets.

library( tidyverse )

source( "../simulation_support_code.R" )
#source( "prepare simulation results/load_simulation_results.R" )
load( "../results/processed_simulation_results.RData" )


#### Describe realized experimental data across scenarios #####

nrow( simdata )
scen = filter( simdata, !duplicated( ID ) )

nrow( scen )
head( scen )
summary( scen$R )
table( scen$J )
table( scen$n.bar )
table( scen$n.bar * scen$J )
table( scen$tau )


##### Prep data by selecting single scenario to explore #####
scenarios = make.scenario.list()

scenarios[3,]

# Subset to targeted simulation
simdata = filter( simdata, ID == 3 )

head( simdata )
dim( simdata )

# Number of simulation runs
table( simdata$method )
length( unique( simdata$subrun ) )
sample( simdata$subrun, 20 )
sample( simdata$run, 20 )
length( unique( simdata$run ) )

# Drop repeated permutations of finite population inference within dataset.  Keep 1 row of evaluation per dataset (and method)
simdata = filter( simdata, endsWith( subrun, "1" ) )
nrow( simdata )

ggplot( simdata, aes( x=method, y=SE.hat ) ) +
  geom_boxplot() +
  coord_flip()


##### Looking at correlation of SE estimates #######

# Look at pairwise SE estimates in the simulation
s2 = simdata %>% dplyr::select( method, SE.hat, ID, run, subrun ) %>% 
  spread( method, SE.hat ) 
head( s2 )
s2 = s2 %>% dplyr::select( `DB-FP-Persons`,FIRC, RIRC, `FE-CR`)
head( s2 )


# Make plots with best fit lines
panel.lm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) {   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)
  abline(0, 1, col = col.lm, ...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y[ok]~ x[ok]), 
           col = col.lm, lty=2, ...)
}

# Look at how SE estimates co-vary.
# Solid line is 45', dashed is best fit
pairs(s2, panel=panel.lm, asp=1)


# Individual pairwise plots
ggplot( s2, aes( `DB-FP-Persons`, `FE-CR` ) ) +
  geom_point() +
  geom_abline( intercept=0, slope=1, col="red" )

ggplot( s2, aes( FIRC, `FE-CR` ) ) +
  geom_point() +
  geom_abline( intercept=0, slope=1, col="red" )

ggplot( s2, aes( FIRC, RIRC ) ) +
  geom_point() +
  geom_abline( intercept=0, slope=1, col="red" )



