
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
load("processed_simulation_results_estimators_estimand_2024-10-02.RData")
load("processed_simulation_results_sim_base_1_2024-10-02.RData")
load("processed_simulation_results_estimators_estimand_2_2024-10-02.RData")
load("processed_simulation_results_performance_1_2024-10-02.RData")
sruns<-bind_rows(sruns.sim_base_1,sruns.estimators_estimand,sruns.estimators_estimand_2,sruns.performance_1)
range.dat<-bind_rows(range.dat.sim_base_1,range.dat.estimators_estimand,range.dat.estimators_estimand_2,range.dat.performance_1)
save(sruns,range.dat, file= "combined_simulation_results.RData")

