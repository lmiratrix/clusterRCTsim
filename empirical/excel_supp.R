# SETUP ----

# load libraries
library(tidyverse)
library(openxlsx)

# load data
source( here::here( "empirical/load_empirical_results.R" ) )

# FORMAT FILES ----

study_char <- combined %>% 
  filter(Covariates == 'best') %>%
  select(Collection:Pretest) %>%
  distinct() %>%
  select(-Covariates)

estimates <- combined %>%
  filter(Covariates == 'best') %>%
  select(UniqueID, Estimator, ATE_hat_std:Estimand)

# SAVE AS EXCEL ----
wb <- createWorkbook() 
addWorksheet(wb, 'Study characteristics')
addWorksheet(wb, 'Estimates')
writeData(wb, 'Study characteristics', study_char)
writeData(wb, 'Estimates', estimates)

saveWorkbook(wb, here::here('empirical/tables/excel_supp.xlsx'))
