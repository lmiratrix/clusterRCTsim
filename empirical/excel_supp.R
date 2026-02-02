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
        select(Collection,	UniqueID,	Project,	Outcome,	Outcome_desc,
               Trt_arm,	Covariates,	n,	J,	nbar,	ncv,	n.25,	n.75,	n.IQR,
               p.tx,	K,	Jbar,	Jcv,	J.25,	J.75,	J.IQR,	n_block,	n_block_cv,
               tx.avg,	tx.cv,	tx.25,	tx.75,	tx.IQR,	cluster_ICC,	block_ICC,
               sdY0,	matched_pairs,	num_singletons,	num_doubletons,	R2.1,	SER2.1,
               ncov.1,	R2.2,	SER2.2,	ncov.2, cov.1,	cov.2,	Pretest,
               lead_author,	year,	journal,	article_url) %>%
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
