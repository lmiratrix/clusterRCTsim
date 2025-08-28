
##########################################
#
# Author: Sophie Litschwartz
# Date Created:11/27/2024
# Project Name: Cluster RCT Estimators
# File Name: generate_grids.R
# File Description: 
#
# generates full set of scenarios, including the statistics on the
# districts for each scenario, for the main cloud-based simulation.
#
################################

rm(list=ls())
library( tidyverse )
library( lme4 )
library( lmerTest )
library( future )
library(furrr)


set.seed( 40444 )

group="fullfactor"


#call the source code

source( "R/data_generation_code.R" )
source( "R/simulation_support_code_fixed_district.R" )

#blocktype=c("messy","same","reasonable")

#Blocks all the same: 
btype1<-data.frame( cv.district=0,
                    sd.ATE=0,
                    size.impact.cor = 0, 
                    nbar.J.cor = 0,   
                    proptx.impact.cor = 0,
                    size.intercept.cor = 0, 
                    variable.p=F,
                    cv.nbar= 0)

#Blocks all messy: 
btype2<-data.frame( cv.district=.85, 
                    sd.ATE=.6, 
                    size.impact.cor = .8, 
                    nbar.J.cor = .8,   
                    proptx.impact.cor = .8, 
                    size.intercept.cor = .8, 
                    variable.p=T,
                    cv.nbar= .33)

                
#Blocks vaguely reasonable
btype3<-data.frame( cv.district=.33,
                    sd.ATE=0.2, 
                    size.impact.cor = .2, 
                    nbar.J.cor = .3,   
                    proptx.impact.cor = 0,
                    size.intercept.cor = .3,
                    variable.p=F,
                    cv.nbar= .33)


#variables for the block scenarios
blockvars <- bind_rows( same = btype1, 
                        reasonable = btype3,
                        messy = btype2, 
                        .id = "blocktype" )

# How do district characteristics vary?
blockvars

#This the parameters that need to go in it
scenarios = expand.grid( 
  K=5, 
  J.bar=c(5,10,20),
  n.bar=20, 
  p.bar=c(0.5,0.7),
  ATE = .2,
  ICC.district = .2,
  num.cov.district = 0,
  blocktype=c("messy","same","reasonable"),
  
  # school-level params
  sd.school =  c(0,.3), 
  ICC.school = c(0,0.20),
  R2.school = .5,
  num.cov.school = c(0,1,5),
  
  # correlation params at school level
  cv.school =  c(0,.33,.6),
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

scenarios <- left_join(scenarios, blockvars, by="blocktype") %>% 
#  select(-blocktype) %>% 
  relocate( ID ) %>%
  as_tibble()


#make the grid 
scenarios$grid <- pmap( dplyr::select( scenarios, -ID, -blocktype ),
                        build_districts )

scenarios

scenarios[131,] %>%
  dplyr::select( blocktype, grid ) %>%
  unnest( cols = grid ) %>%
  dplyr::select( -did, -R2.school )

#scenarios$run = map(scenarios$grid, run.scenario.fixed.district, R=1 , parallel=F)

saveRDS( scenarios, "scenarios.rds" )

