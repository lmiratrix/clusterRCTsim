
##########################################
# Author: Sophie Litschwartz
# Date Created:11/27/2024
# Project Name: Cluster RCT Estimators
# File Name: generate_grids.R
# File Description: 
#generates blocks to be used across parallel iterations of the simulation
#this creates finite block sample simulations 
################################
rm(list=ls())
library( tidyverse )
library( lme4 )
library( lmerTest )
library( future )
library("furrr")


set.seed( 40444 )

group="fullfactor"


#call the source code

source( "R/data_generation_code.R" )
source( "R/simulation_support_code_fixed_district.R" )


#Blocks all the same: 
btype1<-data.frame( variable.J=FALSE, 
                    size.ratio=0,
                    sd.ATE=0,
                    size.impact.cor = 0, 
                    nbar.J.cor = 0,   
                    proptx.impact.cor = 0,
                    size.intercept.cor = 0, 
                    variable.nbar=F, 
                    variable.p=F )

#Blocks all messy: 
#blocktype=c("messy","same","reasonable")
btype2<-data.frame( variable.J=TRUE,
                    size.ratio=1.3, 
                    sd.ATE=.6, 
                    size.impact.cor = .8, 
                    nbar.J.cor = .8,   
                    proptx.impact.cor = .8, 
                    size.intercept.cor = .8, 
                    variable.nbar=T, 
                    variable.p=T )

#variable.nbar=T, variable.p=F
btype3<-data.frame( variable.J=TRUE, 
                    size.ratio=.4,
                    sd.ATE=.1, 
                    size.impact.cor = .2, 
                    nbar.J.cor = .1,   
                    proptx.impact.cor = 0,
                    size.intercept.cor = .1,
                    variable.nbar=T, 
                    variable.p=F )

#variables for the block scenarios
blockvars <- bind_rows( same = btype1, messy = btype2, reasonable = btype3,
                        .id = "blocktype" )

#This the parameters that need to go in it
scenarios = expand.grid( 
  K=5, 
  J.bar=c(7,10,16),
  n.bar=20, 
  p.bar=0.5,
  ATE = .2,
  ICC = .2,
  num.cov.district = 0, 
  blocktype=c("messy","same","reasonable"),
  
  # school-level params
  sd.school =  c(0,.3), 
  ICC.school = c(0,0.20),
  variable.n = c(TRUE,FALSE),
  R2.school = .5,
  num.cov.school = c(0,1,5),
  
  # correlation params at school level
  size.ratio.school =  c(0,.4,1.3),
  size.impact.cor.school = c(0,.8),
  
  # Individual-level params
  num.cov.indiv = 0,
  
  # Make districts as close as possible to desired specified variation
  empirical = TRUE
  )


#fix the R2 
scenarios <- scenarios %>% 
  mutate(R2.school = ifelse(num.cov.school==0, 0, R2.school) )

scenarios$ID = 1:nrow(scenarios)

scenarios <- merge(scenarios, blockvars, by="blocktype", all = T) %>% 
  select(-blocktype) %>% 
  as_tibble()


#make the grid 
scenarios$grid <- pmap( dplyr::select( scenarios, -ID ), build_districts)


#scenarios$run = map(scenarios$grid, run.scenario.fixed.district, R=1 , parallel=F)

saveRDS( scenarios, "scenarios.rds" )

