
##########################################
#
# Author: Sophie Litschwartz
# Date Created:11/27/2024
# Project Name: Cluster RCT Estimators
# File Name: generate_grids.R
# File Description: 
#
# generates full set of scenarios, including the statistics on the
# blocks for each scenario, for the main cloud-based simulation.
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

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_block.R" ) )

#blocktype=c("messy","same","reasonable")

#Blocks all the same: 
btype1<-data.frame( cv.block=0, #checked
                    sd.block=0, #checked
                    size.impact.cor = 0, #checked 
                    nbar.J.cor = 0, #checked   
                    proptx.impact.cor = 0, #checked
                    size.intercept.cor = 0, #checked
                    variable.p=F, #checked
                    cv.nbar= 0 ) #checked

#Blocks all messy: 
btype2<-data.frame( cv.block=.85, #checked
                    sd.block=.6, #checked
                    size.impact.cor = .8, #checked  
                    nbar.J.cor = .8,   #checked
                    proptx.impact.cor = .8, #checked
                    size.intercept.cor = .8, #checked
                    variable.p=T, #checked
                    cv.nbar= .33) #checked

                
#Blocks vaguely reasonable
btype3<-data.frame( cv.block=.33,#checked
                    sd.block=0.2, #checked
                    size.impact.cor = .3,#checked 
                    nbar.J.cor = .3,#checked   
                    proptx.impact.cor = .3,#checked
                    size.intercept.cor = .3,#checked
                    variable.p=T,#checked 
                    cv.nbar= .33) #checked


#variables for the block scenarios
blockvars <- bind_rows( same = btype1, 
                        reasonable = btype3,
                        messy = btype2, 
                        .id = "blocktype" )

# How do block characteristics vary?
blockvars

#This the parameters that need to go in it
scenarios = expand.grid( 
  K=5, #checked
  J.bar=c(5,10,20), #checked
  n.bar=20, #checked
  p.bar=c(0.5,0.7), #checked
  ATE = .2, #checked
  ICC.block = .2, #checked
  num.cov.block = 0, #checked
  blocktype=c("messy","same","reasonable"),
  
  # cluster-level params
  sd.cluster =  c(0,.3), #checked
  ICC.cluster = c(0,0.20), #checked
  R2.cluster = .5, #checked
  num.cov.cluster = c(0,1,5), #checked
  
  # correlation params at cluster level
  cv.cluster =  c(0,.33,.6), #checked
  size.impact.cor.cluster = c(0,.8), #checked
  
  # Individual-level params
  num.cov.indiv = 0,
  
  # Make blocks as close as possible to desired specified variation
  empirical = TRUE
  )


#fix the R2 
scenarios <- scenarios %>% 
  mutate(R2.cluster = ifelse(num.cov.cluster==0, 0, R2.cluster) )

scenarios$ID = 1:nrow(scenarios)

scenarios <- left_join(scenarios, blockvars, by="blocktype") %>% 
#  select(-blocktype) %>% 
  relocate( ID ) %>%
  as_tibble()

if ( FALSE ) {
  # some testing code
  key_ID <- 1067
  set.seed( 44040 )
  debug( build_blocks )
  grid <- do.call( build_blocks, 
           scenarios[ scenarios$ID == key_ID, ] %>%
             dplyr::select( -ID, -blocktype ) %>% 
             as.list() )
  
  describe_block_grid(grid)
  
}

#make the grid 
scenarios$grid <- pmap( dplyr::select( scenarios, -ID, -blocktype ),
                        build_blocks )

scenarios

# Check a specific scenario
scenarios[131,] %>%
  dplyr::select( blocktype, grid ) %>%
  unnest( cols = grid ) %>%
  dplyr::select( -did, -R2.cluster )

#scenarios$run = map(scenarios$grid, run.scenario.fixed.block, R=1 , parallel=F)

scenarios$group = group
scenarios$date.created=format(Sys.Date(), "%m/%d/%Y")
saveRDS( scenarios,
         here::here( "results/scenarios.rds" ) )

