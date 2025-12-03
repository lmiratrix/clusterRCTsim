
##########################################
# Author: Sophie Litschwartz
# Date Created:10/16/2024
# Project Name: Cluster RCT Estimators
# File Name: Combine_Results_Across_Simulations.R
# File Description: 
#Pulls simulation results from each group of data generating processes 
#Combines them into two master files: one that's an aggregation of range results 
#and one that's aggregated data collapsed by estimator 
#Master file is then saved as combined_simulation_results.RData in a results_sav folder 
################################


###Combining Results Across Simulations
rm(list=ls())
library(tidyverse)
load("processed_simulation_results_agg_fullfactor_2025-11-25.RData")
scenarios.fullfactor<-scenarios
load("processed_simulation_results_agg_lowdf_2025-07-01.RData")
scenarios.lowdf<-scenarios
sruns<-bind_rows(sruns.fullfactor, sruns.lowdf)
scenarios<-bind_rows(scenarios.fullfactor, scenarios.lowdf)
save(scenarios,sruns, file= "combined_simulation_results.RData")

