##########################################
#
# Author: Sophie Li.cluwartz
# Date Created:06/13/2025
# Project Name: Cluster RCT Estimators
# File Name: generate_grids_lowdf.R
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

group="lowdf"


#call the source code

source( here::here( "R/data_generation_code.R" ) )
source( here::here( "R/simulation_support_code_fixed_block.R" ) )



#blocktype=c("messy","same","reasonable")

#Blocks all the same: 
btype1<-data.frame( sd.block=0,
                    size.impact.cor = 0, 
                    size.intercept.cor = 0)

#Blocks all messy: 
btype2<-data.frame( sd.block=.6, 
                    size.impact.cor = .8, 
                    size.intercept.cor = .8)

#Blocks vaguely reasonable
btype3<-data.frame( sd.block=0.2, 
                    size.impact.cor = .2, 
                    size.intercept.cor = .3)

#variables for the block scenarios
blockvars <- bind_rows( same = btype1, 
                        reasonable = btype3,
                        messy = btype2, 
                        .id = "blocktype" )

# How do block characteristics vary?
blockvars

# Define the four specific scenario sets
scenario_sets <- data.frame(
  set_id = 1:4,
  K = c(20, 10, 7, 20),
  J.bar = c(2.4, 4, 4, 2),
  num.cov.cluster = c(5, 6, 6, 0),
  cv.block = c(0.33, 0.33, 0, 0),
  # For sets 3 and 4, need special handling of size correlations
  special_zero_cors = c(FALSE, FALSE, TRUE, TRUE)
)

#This the parameters that need to go in
scenarios = expand.grid( 
  set_id = 1:4,

  # Fixed block-level factors (same for all scenarios)
  n.bar=20, 
  ATE = .2,
  ICC.block = .2,
  num.cov.block = 0,
  nbar.J.cor = 0,   
  proptx.impact.cor = 0,
  variable.p = FALSE,
  cv.nbar = 0,
  p.bar=c(0.5),
  
  blocktype=c("messy","same","reasonable"),
  
  # cluster-level params
  sd.cluster =  c(0,.3), 
  ICC.cluster = c(0,0.20),
  R2.cluster = .5,
  
  # correlation params at cluster level
  cv.cluster =  c(0,.33,.6),
  size.impact.cor.cluster = c(0,.8),
  
  # Individual-level params
  num.cov.indiv = 0,
  
  # Make blocks as close as possible to desired specified variation
  empirical = TRUE
)

# Join with scenario sets to get K, J.bar, num.cov.cluster
scenarios <- scenarios %>%
  left_join(scenario_sets %>% select(set_id, K, J.bar, num.cov.cluster, cv.block, special_zero_cors), by = "set_id")


# Filter out invalid combinations where size.impact.cor.cluster > 0 but sd.cluster = 0 or cv.cluster = 0
scenarios <- scenarios %>%
  filter(!(size.impact.cor.cluster > 0 & (sd.cluster == 0 | cv.cluster == 0)))

# Fix the R2.cluster - set to NULL (0) if num.cov.cluster == 0
scenarios <- scenarios %>% 
  mutate(R2.cluster = ifelse(num.cov.cluster == 0, 0, R2.cluster))

scenarios$ID = 1:nrow(scenarios)

# Join with block variables and handle special cases for sets 3 and 4
scenarios <- scenarios %>%
  left_join(blockvars, by="blocktype") %>% 
  # For sets 3 and 4, override size correlations to be 0 for "messy" and "reasonable" blocks
  mutate(
    size.impact.cor = ifelse(special_zero_cors & blocktype %in% c("messy", "reasonable"), 0, size.impact.cor),
    size.intercept.cor = ifelse(special_zero_cors & blocktype %in% c("messy", "reasonable"), 0, size.intercept.cor)
  ) %>%
  select(-special_zero_cors, -set_id) %>%
  relocate( ID ) %>%
  as_tibble()

#make the grid 
scenarios$grid <- pmap( dplyr::select( scenarios, -ID, -blocktype ),
                        build_blocks )

#save the true J.bar
scenarios$J.bar.true <- map_dbl(scenarios$grid, ~ mean(.x$J))

scenarios

scenarios[131,] %>%
  dplyr::select( blocktype, grid ) %>%
  unnest( cols = grid ) %>%
  dplyr::select( -did, -R2.cluster )

#scenarios$run = map(scenarios$grid, run.scenario.fixed.block, R=1 , parallel=F)

scenarios$group = group
saveRDS( scenarios, 
         here::here( "scenarios_lowdf.rds" ) )