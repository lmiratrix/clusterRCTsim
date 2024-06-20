


# Load libraries and data ----

library(tidyverse)
library(readxl)
library(gridExtra)
library(openxlsx)

combined <-  readRDS('/data/share/pse/IES CRT/sd7/master240528_combined.RDS') 

# Research question 1 ----
# To what extent can the choice of estimator of the overall average treatment effect (β) result in a different impact estimate?
#    a. To what extent do effect estimates vary across all estimators? 
#    b. For different types of estimands (e.g., site average vs individual, or finite population vs not),
#       to what extent do effect estimates vary across estimators that target that type of estimand? 

## 1a. Range of effect estimates across all estimands ----

ate_summary <- combined %>%
  group_by(Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T)) %>%
  mutate(Range = Max - Min)

write.xlsx(ate_summary, '/data/share/pse/IES CRT/tables/ate_summary.xlsx')

makeDotplot <- function(df, var, 
                        xlab       = NULL, 
                        ylab       = NULL, 
                        title      = NULL, 
                        subtitle   = NULL, 
                        binw       = NULL,
                        filename   = NULL,
                        coordratio = NULL,
                        byest      = F) {
  p <- ggplot(df, aes(x = !!sym(var))) +
    geom_dotplot(binwidth   = binw,
                 dotsize    = 0.8,
                 stackratio = 1.2,
                 method     = 'histodot') +
    theme_minimal() +
    labs(x        = xlab,
         y        = ylab,
         title    = title,
         subtitle = subtitle) +
    theme(axis.text.y        = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())
  
  if (!is.null(coordratio)) {
    p <- p +
      coord_fixed(ratio = coordratio)
  }
  
  if (byest == T) {
    p <- p +
      facet_wrap(~ Estimand) +
      theme(strip.background = element_rect(fill = "lightgray", 
                                            linewidth = 0))
  }
  
  if (!is.null(filename)) {
    ggsave(paste0('/data/share/pse/IES CRT/tables/', filename, '.png'), plot = p)
  }
  
  return(p)
}

### With covariates ----
ate_summary %>%
  filter(Covariates == 1) %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              subtitle = 'With covariates',
              filename = 'ate_range_covar_240220') 

ate_summary %>%
  filter(Covariates == 1) %>%
  makeDotplot(var        = 'Range',
              xlab       = 'Range of Point Estimates of β',
              ylab       = 'Number of Outcomes',
              filename   = 'paper/ate_range_covar_paper_240513',
              coordratio = .15)

### Without covariates ----
ate_summary %>%
  filter(Covariates == 0) %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'ate_range_240220')

## 1b. Effect estimates by estimand ----

ate_summary_estimand <- combined %>%
  # # combining finit vs super population estimands 
  # mutate(Estimand = str_replace_all(Estimand, 'SP\\*', 'SP')) %>%
  group_by(Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            WhichMin = Estimator[which.min(SE_hat_std)],
            WhichMax = Estimator[which.max(SE_hat_std)]) %>%
  mutate(Range = Max - Min)

write.xlsx(ate_summary_estimand, '/data/share/pse/IES CRT/tables/ate_summary_estimand.xlsx')

ate_summary_estimand %>%
  filter(Covariates == 1) %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP(\\*|)_')) %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β',
              ylab     = 'Number of Outcomes',
              filename = 'paper/ate_range_covar_byestimand_paper_240513',
              coordratio = .15,
              byest    = T)

# Research question 2 ----
# To what extent can the choice of estimator of the standard error of the overall 
# average effect (SEˆ(β)) result in a different standard error? 
#   a. To what extent do the estimated standard errors vary across all estimators? 
#   b. For different types of estimands, to what extent do estimated standard errors 
#      vary across estimators that target that type of estimand? 

se_summary <- combined %>%
  # replace values of Inf with NA
  mutate(SE_hat_std = ifelse(is.infinite(SE_hat_std), NA, SE_hat_std)) %>%
  group_by(Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = Max / Min)

write.xlsx(se_summary, '/data/share/pse/IES CRT/tables/se_summary.xlsx')

## 2a. Range of SE estimates across all estimators ----

### With covariates ----
se_summary %>%
  filter(Covariates == 1) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'With covariates',
              filename = 'se_ratio_covar_240220')

se_summary %>%
  filter(Covariates == 1) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'paper/se_ratio_covar_paper_240513') 

# remove outliers above 50
se_summary %>%
  filter(Covariates == 1) %>%
  filter(Ratio < 50) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_rmoutlier_240513') 

### Without covariates ----
se_summary %>%
  filter(Covariates == 0) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'se_ratio_240220')

# making version without outlier
se_summary %>%
  filter(Covariates == 0) %>%
  filter(Ratio < 100) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'se_ratio_rmoutlier_240220')

## 2b. SE estimates by estimand ----

se_summary_estimand <- combined %>%
  # # combining finit vs super population estimands 
  # mutate(Estimand = str_replace_all(Estimand, 'SP\\*', 'SP')) %>%
  group_by(Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T),
            WhichMin = Estimator[which.min(SE_hat_std)],
            WhichMax = Estimator[which.max(SE_hat_std)]) %>%
  mutate(Ratio = Max / Min) 

write.xlsx(se_summary_estimand, '/data/share/pse/IES CRT/tables/se_summary_estimand.xlsx')

se_summary_estimand %>%
  filter(Covariates == 1) %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP(\\*|)_')) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'paper/se_ratio_covar_byestimand_paper_240513',
              byest    = T,
              coordratio = 60) 

# remove outliers above 50
se_summary_estimand %>%
  filter(Covariates == 1) %>%
  filter(Ratio < 50) %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP_')) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_byestimand_rmoutlier_240513',
              byest    = T,
              coordratio = 10) 

# Other checks ----

## QC effect size estimates ----

# read in Excel file with reported effect sizes
report_es <- read_excel('/data/share/pse/IES CRT/tables/Report vs Our Outcomes 240429.xlsx') %>%
  select(Project = `Study Name`,
         Outcome_desc = Outcome,
         ATE_hat_std_qc = `Report Effect size`) %>%
  mutate(Trt_arm = case_when(Project == 'EMERGE - BB' ~ 'BB',
                             Project == 'EMERGE - BBSSR' ~ 'BBSSR',
                             Project == 'CARES - Incredible Years' ~ 'Incredible Years',
                             Project == 'CARES - Preschool PATHS' ~ 'Preschool PATHS',
                             Project == 'CARES - Tools of the Mind' ~ 'Tools of the Mind',
                             TRUE ~ NA_character_),
         Outcome_desc = case_when(grepl('F1_TR_WRS', Outcome_desc) ~ combined %>% 
                                    filter(grepl('F1_TR_WRS', Outcome_desc)) %>% 
                                    pull(Outcome_desc) %>% 
                                    unique(),
                                  grepl('Cumulative credits', Outcome_desc) ~ combined %>% 
                                    filter(grepl('Cumulative credits', Outcome_desc)) %>% 
                                    pull(Outcome_desc) %>% 
                                    unique(),
                                  TRUE ~ Outcome_desc),
         Grade_lvl = ifelse(grepl('Content Literacy Continuum', Project),
                            str_remove(Project, 'Content Literacy Continuum - '),
                            NA),
         Project = case_when(Project == 'Making pre-K Count' ~ 'MPC',
                             Project == 'Building Blocks - TRIAD' ~ 'TRIAD',
                             Project == 'Building Blocks - TRIAD Scale-Up' ~ 'TRIAD Scale-up',
                             grepl('EMERGE', Project) ~ 'EMERGE',
                             grepl('CARES', Project) ~ 'CARES',
                             Project == 'Tennessee Class Size Experiment (STAR)' ~ 'Tennessee Star',
                             grepl('Content Literacy Continuum', Project) ~ 'CLC',
                             TRUE ~ Project))

# merge onto combined results (using covars)
esqc <- combined %>%
  filter(Covariates == 1) %>%
  select(Project, Trt_arm, Outcome, Outcome_desc, ATE_hat_std, Estimator) %>%
  full_join(report_es,
            relationship = 'many-to-many')

# plot reported effect size against all of ours
qcplot_list <- list()
list_index <- 1
for (i in 1:nrow(report_es)) {
  if (report_es$Grade_lvl[i] == 'Grade 9' | is.na(report_es$Grade_lvl[i])) {
    proj <- report_es$Project[i]
    trt  <- report_es$Trt_arm[i]
    y    <- report_es$Outcome_desc[i]
    var  <- esqc %>%
      filter(Project == proj) %>%
      filter(Outcome_desc == y) %>%
      pull(Outcome) %>%
      unique()
    if (proj == 'CLC') {
      qc <- report_es$ATE_hat_std_qc[report_es$Outcome_desc == y]
      qc_color <- c('blue', 'red')
    } else {
      qc <- report_es$ATE_hat_std_qc[i]
      qc_color <- 'red'
    }
    
    
    plot_title <- paste0('Study: ', proj, 
                         ifelse(is.na(trt), '', paste0(' (Arm: ', trt, ')')),
                         ifelse(proj == 'CLC', ' (blue = grade 9, red = grade 10)', ''))
    plot_subtitle <- paste0('Outcome: ', var)
    # plot_object <- paste0('ggplot_', i)
    
    qcplot_list[[list_index]] <- esqc %>%
             rowwise() %>%
             filter(Project == proj) %>%
             filter(ifelse(is.na(trt), TRUE, Trt_arm == trt)) %>%
             filter(Outcome_desc == y) %>% 
             filter(ifelse(proj == 'CLC', Grade_lvl == 'Grade 9', TRUE)) %>%
             ggplot(aes(x = ATE_hat_std)) +
             geom_dotplot(binwidth   = NULL,
                          dotsize    = 0.8,
                          stackratio = 1.2,
                          method     = 'histodot') +
             theme_minimal() +
             labs(x        = 'Effect size estimate',
                  y        = NULL,
                  title    = str_wrap(plot_title, 25),
                  subtitle = plot_subtitle) +
             theme(axis.text.y  = element_blank()) +
             geom_vline(xintercept = qc,
                        color = qc_color)
    
    list_index <- list_index + 1
  }
}

plot_grid <- do.call('grid.arrange', qcplot_list)

ggsave(filename = '/data/share/pse/IES CRT/tables/QC/QC_EffectSize_240429.png',
       plot = plot_grid,
       width = 15,
       height = 15)

## Outlier investigation ----

# modeling approach
# m_ate <- lm(ATE_hat_std ~ 1 + Outcome,
#             data      = combined,
#             na.action = na.exclude)
# 
# m_se <- lm(SE_hat_std ~ 1 + Outcome,
#            data      = combined,
#            na.action = na.exclude)
# 
# combined$mresid_ate <- resid(m_ate)
# combined$mresid_se  <- resid(m_se)
# 
# ggplot(combined %>% filter(Covariates == 1), aes(mresid_ate, fill = Outcome)) +
#   facet_wrap(~ Estimator) +
#   geom_dotplot()
# 
# ggplot(combined %>% filter(Covariates == 1), aes(mresid_se, fill = Outcome)) +
#   facet_wrap(~ Estimator) +
#   geom_dotplot()

# median approach
combined2 <- combined %>%
  # filter(Covariates == 1) %>%
  group_by(UniqueID, Covariates) %>%
  mutate(median_ate   = median(ATE_hat_std, na.rm = T),
         medresid_ate = ATE_hat_std - median_ate)

(ggplot(combined2 %>% filter(Covariates == 1), aes(medresid_ate)) +
  facet_wrap(~ Estimator) +
  geom_dotplot() +
  theme(legend.position = 'none') +
  labs(title    = 'Variation in ATE by estimator',
       subtitle = 'With covariates',
       x        = 'Difference from median ATE') +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())) %>%
  ggsave(paste0('/data/share/pse/IES CRT/tables/ATE_median_byestimator_covar_240220.png'), plot = .,
         width  = 9,
         height = 7)

(ggplot(combined2 %>% filter(Covariates == 1) %>% mutate(Trt_arm = ifelse(is.na(Trt_arm), '', Trt_arm)), aes(medresid_ate)) +
  facet_wrap(vars(Project, Outcome, Trt_arm)) +
  geom_dotplot() +
  theme(legend.position = 'none') +
  labs(title    = 'Variation in ATE by UniqueID',
       subtitle = 'With covariates',
       x        = 'Difference from median ATE') +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text   = element_text(size = 6, margin = margin(.1, .1, .1, .1, 'mm')))) %>%
  ggsave(paste0('/data/share/pse/IES CRT/tables/ATE_median_byuniqueID_covar_240220.png'), plot = .,
         width  = 11,
         height = 8) 

## Assessing impact of covariates ----

ggplot(combined2 %>% 
         mutate(Trt_arm = ifelse(is.na(Trt_arm), '', Trt_arm)), 
       aes(x = medresid_ate,
           fill = factor(Covariates), group = factor(Covariates),
           alpha = factor(Covariates))) +
  facet_wrap(vars(Project, Outcome, Trt_arm)) +
  geom_dotplot() +
  labs(title    = 'Variation in ATE by UniqueID',
       subtitle = 'With and without covariates',
       x        = 'Difference from median ATE') +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text   = element_text(size = 6, margin = margin(.1, .1, .1, .1, 'mm')),
        legend.position = 'bottom')+
  scale_fill_manual(values = c('black', 'black'), name = 'Covariates') +
  scale_alpha_manual(values = c(.2, 1), name = 'Covariates')

# SREE proposal dotplot ----

ate_summary_sree <- combined %>%
  filter(Estimator != 'LR_FI_CRVE_Cluster') %>%
  group_by(Project, Trt_arm, Outcome, Covariates) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T)) %>%
  mutate(Range = Max - Min)

ate_summary_sree %>%
  filter(Covariates == 1) %>%
  makeDotplot(var      = 'Range',
              xlab     = paste0('Range across ', max(.$N), ' Point Estimates'), 
              ylab     = 'Count', 
              title    = 'Range of ATE Estimates across all Estimators',
              subtitle = 'Each point represents one study-treatment arm-outcome combination',
              filename = 'ate_range_covar_sree2024')
