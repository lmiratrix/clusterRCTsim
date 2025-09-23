# SETUP ----

# load libraries
library(tidyverse)
library(openxlsx)

# load data
source( here::here( "empirical/load_empirical_results.R" ) )

# # export study characteristics file with fields to fill out for curated datasets that are missing it
# combined %>% 
#   select(Collection:Trt_arm, lead_author, year, journal, article_url) %>%
#   distinct() %>%
#   writexl::write_xlsx('/data/share/pse/IES CRT/doc/study_research_info.xlsx')

# FORMAT FILES ----

study_char <- combined %>% 
  select(Collection:Pretest, lead_author, year, journal, article_url) %>%
  distinct()

estimates <- combined %>%
  select(UniqueID, Covariates, Estimator, ATE_hat_std:Estimand)

# SAVE AS EXCEL ----
wb <- createWorkbook() 
addWorksheet(wb, 'Study characteristics')
addWorksheet(wb, 'Estimates')
writeData(wb, 'Study characteristics', study_char)
writeData(wb, 'Estimates', estimates)

saveWorkbook(wb, here::here('empirical/tables/excel_supp.xlsx'), overwrite = T)
