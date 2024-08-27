

source( here::here( "empirical/setup.R" ) )


# Research question 1 ----
# To what extent can the choice of estimator of the overall average treatment effect (β) result in a different impact estimate?
#    a. To what extent do effect estimates vary across all estimators? 
#    b. For different types of estimands (e.g., site average vs individual, or finite population vs not),
#       to what extent do effect estimates vary across estimators that target that type of estimand? 

## 1a. Range of effect estimates across all estimands ----

ate_summary <- combined %>%
  group_by(Collection,Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = "drop") %>%
  mutate(Range = Max - Min)
ate_summary

write.xlsx(ate_summary, plot_path( 'ate_summary.xlsx') )


#' Make plot of differences and save to file
#' 
makeDotplot <- function(df, var, 
                        xlab       = NULL, 
                        ylab       = NULL, 
                        title      = NULL, 
                        subtitle   = NULL, 
                        binw       = NULL,
                        filename   = NULL,
                        coordratio = NULL,
                        axisfs     = NULL,
                        xscale     = NULL,
                        plotmargin = NULL,
                        rmlegend   = F,
                        byest      = F) {
  
  p <- ggplot(df, aes(x = !!sym(var))) +
    geom_dotplot( aes( fill=Collection ),
                  binwidth   = binw,
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
      facet_grid(~ Estimand) +
      theme(strip.background = element_rect(fill = "lightgray", 
                                            linewidth = 0),
            panel.spacing = unit(2, "lines"))
  }
  
  if (rmlegend == T) {
    p <- p +
      theme(legend.position = 'none')
  }
  
  if (!is.null(axisfs)) {
    p <- p +
      theme(axis.text  = element_text(size = axisfs * 0.8),
            axis.title = element_text(size = axisfs),
            strip.text = element_text(size = axisfs))
  }
  
  if (!is.null(xscale)) {
    p <- p +
      scale_x_continuous(breaks = xscale)
  }
  
  if (!is.null(plotmargin)) {
    p <- p +
      theme(plot.margin = margin(plotmargin))
  }
  
  if (!is.null(filename)) {
    ggsave(paste0(PLOT_PATH, filename, '.png'), plot = p)
  }
  
  return(p)
}




### With covariates ----
ate_summary %>%
  filter(Covariates == "all") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              subtitle = 'With covariates',
              filename = 'ate_range_covar_240220',
              binw = 0.02) 

ate_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var        = 'Range',
              xlab       = 'Range of Point Estimates of β',
              ylab       = 'Number of Outcomes',
              filename   = 'ate_range_covar_paper_240513',
              coordratio = .15)



### Without covariates ----
ate_summary %>%
  filter(Covariates == "none") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'ate_range_240220')


### With "best" covariate
ate_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              subtitle = 'With single best covariate', 
              filename = 'ate_range_best')





## 1b. Effect estimates by estimand ----

# calculate proportion missing--some estimands have no defined
# estimators for some of the outcomes
combined %>%
  group_by(Project, Trt_arm, Outcome, Covariates, Estimand) %>%
  summarize(n = n(),
            N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            .groups = "drop") %>%
  mutate(PropMiss = Nmiss / (N + Nmiss)) %>%
  filter(PropMiss > 0) %>%
  arrange(desc(PropMiss))

# What estimands are in the dataset?
table( combined$Estimand )

ate_summary_estimand <- combined %>%
  # combining finite vs super population estimands 
  mutate(Estimand = str_replace_all(Estimand, 'SP\\*', 'SP')) %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            #WhichMin = Estimator[which.min(SE_hat_std)],
            #WhichMax = Estimator[which.max(SE_hat_std)],
            .groups = "drop") %>%
  mutate(Range = Max - Min)

write.xlsx(ate_summary_estimand, plot_path( 'ate_summary_estimand.xlsx' ) )

ate_summary_estimand %>%
  filter(Covariates == "best") %>%
  mutate(Estimand = str_remove_all(Estimand, 'SP(\\*|)_')) %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β',
              ylab     = 'Number of Outcomes',
              filename = 'ate_range_covar_byestimand_paper_240513',
              coordratio = .15,
              byest    = TRUE )



# Research question 2 ----
# To what extent can the choice of estimator of the standard error of the overall 
# average effect (SEˆ(β)) result in a different standard error? 
#   a. To what extent do the estimated standard errors vary across all estimators? 
#   b. For different types of estimands, to what extent do estimated standard errors 
#      vary across estimators that target that type of estimand? 

se_summary <- combined %>%
  # replace values of Inf with NA
  mutate(SE_hat_std = ifelse(is.infinite(SE_hat_std), NA, SE_hat_std)) %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = Max / Min)

write.xlsx(se_summary, plot_path( 'se_summary.xlsx') )


## 2a. Range of SE estimates across all estimators ----

### With covariates ----
se_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'With covariates',
              filename = 'se_ratio_covar_240220')

se_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_paper_240513') 

# remove outliers above 50
se_summary %>%
  filter(Covariates == "best") %>%
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
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T),
            #WhichMin = Estimator[which.min(SE_hat_std)],
            #WhichMax = Estimator[which.max(SE_hat_std)], 
            .groups ="drop" ) %>%
  mutate(Ratio = Max / Min) 

write.xlsx(se_summary_estimand, plot_path( 'se_summary_estimand.xlsx') )

se_summary_estimand %>%
  filter(Covariates == "best") %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP(\\*|)_')) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_byestimand_paper_240513',
              byest    = T,
              coordratio = 60) 

# remove extreme outliers
se_summary_estimand %>%
  filter(Covariates == "best") %>%
  filter(Ratio < 10) %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP_')) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_byestimand_rmoutlier_240513',
              byest    = T ) 




# Other checks ----

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
# ggplot(combined %>% filter(Covariates == "best"), aes(mresid_ate, fill = Outcome)) +
#   facet_wrap(~ Estimator) +
#   geom_dotplot()
# 
# ggplot(combined %>% filter(Covariates == "best"), aes(mresid_se, fill = Outcome)) +
#   facet_wrap(~ Estimator) +
#   geom_dotplot()

# median approach
combined2 <- combined %>%
  # filter(Covariates == "best") %>%
  group_by(UniqueID, Covariates) %>%
  mutate(median_ate   = median(ATE_hat_std, na.rm = T),
         medresid_ate = abs( ATE_hat_std - median_ate) )

plt <- ggplot( filter(combined2, Covariates == "best"), 
               aes(medresid_ate)) +
  facet_wrap(~ Estimator) +
  geom_dotplot() +
  theme(legend.position = 'none') +
  labs(title    = 'Variation in ATE by estimator',
       subtitle = 'With covariates',
       x        = 'Difference from median ATE') +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())
plt
ggsave(plot_path('ATE_median_byestimator_covar_240220.png'),
       width  = 9,
       height = 7)

c2arm <- combined2 %>%
  filter(Covariates == "best") %>% 
  mutate(Trt_arm = ifelse(is.na(Trt_arm), '', Trt_arm))
plt2 <- ggplot( c2arm, aes(medresid_ate)) +
  facet_wrap(vars(Project, Outcome, Trt_arm)) +
  geom_dotplot() +
  theme(legend.position = 'none') +
  labs(title    = 'Variation in ATE by UniqueID',
       subtitle = 'With covariates',
       x        = 'Difference from median ATE') +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text   = element_text(size = 6, margin = margin(.1, .1, .1, .1, 'mm')))
plt2
ggsave(paste0('/data/share/pse/IES CRT/tables/ATE_median_byuniqueID_covar_240220.png'),
       width  = 11,
       height = 8) 



## Assessing impact of covariates ----

ggplot( combined2, 
        aes(x = medresid_ate,
            fill = factor(Covariates), group = factor(Covariates),
            alpha = factor(Covariates))) +
  facet_wrap( ~ UniqueID ) +
  geom_dotplot() +
  labs(title    = 'Variation in ATE by UniqueID',
       subtitle = 'With and without covariates',
       x        = 'Difference from median ATE') +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text   = element_text(size = 6, margin = margin(.1, .1, .1, .1, 'mm')),
        legend.position = 'bottom')+
  scale_fill_manual(values = c('black', 'black', 'black'), name = 'Covariates') +
  scale_alpha_manual(values = c(1, 0.5, 0.2), name = 'Covariates')




ggplot( combined2, 
        aes( x = Covariates, y = medresid_ate ) ) +
  facet_wrap( ~ Estimator ) +
  geom_boxplot() +
  coord_flip() 



# (Old) SREE proposal dotplot ----

# ate_summary_sree <- combined %>%
#   filter(Estimator != 'LR_FI_CRVE_Cluster') %>%
#   group_by(Project, Trt_arm, Outcome, Covariates) %>%
#   summarize(N     = n(),
#             Nmiss = sum(is.na(ATE_hat_std)),
#             Min   = min(ATE_hat_std, na.rm = T),
#             Max   = max(ATE_hat_std, na.rm = T)) %>%
#   mutate(Range = Max - Min)

# ate_summary_sree %>%
#   filter(Covariates == "best") %>%
#   makeDotplot(var      = 'Range',
#               xlab     = paste0('Range across ', max(.$N), ' Point Estimates'), 
#               ylab     = 'Count', 
#               title    = 'Range of ATE Estimates across all Estimators',
#               subtitle = 'Each point represents one study-treatment arm-outcome combination',
#               filename = 'ate_range_covar_sree2024')

# SREE 2024 presentation ----

## ATE range, best covariate
ate_summary_sree <- combined %>%
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$', Estimator)) %>%
  filter(Collection == 'Curated') %>%
  group_by(Collection, Project, Trt_arm, Outcome, Covariates, J, K, ncov.1) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = 'drop') %>%
  mutate(Range = Max - Min,
         df    =  J - (2 * K) - ncov.1 - 1)

ate_summary_sree %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = paste0('Range across ', max(.$N), ' Point Estimates'),
              ylab     = 'Count',
              # title    = 'Range of ATE Estimates across all Estimators',
              # subtitle = 'Each point represents one study-treatment arm-outcome combination',
              filename = 'ate_range_bestcovar_sree2024',
              rmlegend = T,
              axisfs   = 20,
              coordratio = .1,
              binw     = .01)

## ATE range by estimand, best covariate
ate_summary_estimand_sree <- combined %>%
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$', Estimator)) %>%
  filter(Collection == 'Curated') %>%
  group_by(Collection, Project, Trt_arm, Outcome, Covariates, Estimand) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = 'drop') %>%
  mutate(Range = Max - Min)

ate_summary_estimand_sree %>%
  filter(Covariates == "best") %>%
  mutate(Estimand = case_when(Estimand == 'SP_Cluster_Cluster' ~ 'Cluster Weighted Estimand',
                              Estimand == 'SP_Person_Person'   ~ 'Person Weighted Estimand')) %>%
  makeDotplot(var        = 'Range',
              xlab       = 'Range across Point Estimates',
              ylab       = 'Count',
              # title      = 'Range of ATE Estimates across all Estimators, Grouped by Estimand',
              # subtitle   = 'Each point represents one study-treatment arm-outcome combination',
              filename   = 'ate_range_bestcovar_estimand_sree2024',
              coordratio = 0.10,
              byest      = T,  
              rmlegend   = T,
              axisfs     = 20,
              binw       = .01,
              xscale     = seq(0, .2, .02),
              plotmargin = c(5, 20, 5, 20))

## DF vs range, best vs. all covar
(ate_summary_sree %>%
    filter(Covariates %in% c('best', 'all')) %>%
    filter(df < 100) %>%
    mutate(Covariates = ifelse(Covariates == 'best', 'Best covariate', 'All covariates')) %>%
    ggplot(aes(x = df, y = Range)) +
    geom_point() + 
    facet_wrap(~ Covariates) +
    labs(x = "Degrees of Freedom", y = "Range") +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", 
                                          linewidth = 0),
          axis.title = element_text(size = 20),
          axis.text  = element_text(size = 16),
          strip.text = element_text(size = 20),
          panel.spacing = unit(2, "lines"))) %>%
  ggsave(plot     = ., 
         filename = plot_path('dfvsrange_sree2024.png'),
         width = 14)

