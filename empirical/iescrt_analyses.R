

source( here::here( "empirical/load_empirical_results.R" ) )

# Drop the block estimands and the HT/Raj estimators
table( combined$Estimand, useNA = 'always' )
combined <- combined %>%
  dplyr::filter( !(Estimand %in% c( 'SP_Person_Block', 'SP_Cluster_Block', "SP*_Cluster_Block" ) ) ) %>%
  dplyr::filter( !(Estimator %in% c( "DB_HT", "DB_Raj" ) ) )


# Research question 1 ----
# To what extent can the choice of estimator of the overall average treatment effect (β) result in a different impact estimate?
#    a. To what extent do effect estimates vary across all estimators? 
#    b. For different types of estimands (e.g., site average vs individual, or finite population vs not),
#       to what extent do effect estimates vary across estimators that target that type of estimand? 

## 1a. Range of effect estimates across all estimands ----

ate_summary <- combined %>%
  group_by(UniqueID, Collection,Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = "drop") %>%
  mutate(Range = Max - Min)
ate_summary


# How many outcomes do we have?
all_outcomes <- ate_summary %>%
  dplyr::select( UniqueID, Collection, Project ) %>%
  unique()
table( all_outcomes$Collection )

# How many studies do we have?
all_studies <- ate_summary %>%
  dplyr::select( Collection, Project ) %>%
  unique() %>%
  print( n = 100)
all_studies
table( all_studies$Collection )

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
                        tcbreaks   = NULL,
                        tcratio    = F,
                        dim        = NULL,
                        rmlegend   = F,
                        byest      = F) {
  
  p <- ggplot(df, aes(x = !!sym(var))) +
    geom_dotplot( aes( fill=Collection ),
                  binwidth   = binw,
                  dotsize    = 0.8,
                  stackratio = 1.2,
                  stackgroups = TRUE,
                  method     = 'histodot') +
    theme_minimal() +
    labs(x        = xlab,
         y        = ylab,
         title    = title,
         subtitle = subtitle) +
    theme(axis.text.y        = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  
  if (!is.null(coordratio)) {
    p <- p +
      coord_fixed(ratio = coordratio)
  }
  
  if (byest == T) {
    p <- p +
      facet_grid(~ Estimand) +
      theme(strip.background = element_rect(fill = "lightgray"),
            panel.spacing = unit(2, "lines"))
  }
  
  if (rmlegend == T) {
    p <- p +
      theme(legend.position = 'none')
  }
  
  
  if (!is.null(xscale)) {
    p <- p +
      scale_x_continuous(breaks = xscale) +
      xlim(c(min(xscale), max(xscale)))
  }
  
  if (!is.null(plotmargin)) {
    p <- p +
      theme(plot.margin = margin(plotmargin))
  }
  
  if (!is.null(tcbreaks)) {
    if (tcratio) {
      p <- p +
        scale_x_continuous(breaks = tcbreaks, 
                           labels = c(paste0(tcbreaks[1:(length(tcbreaks) - 1)], 'X'), paste0(tcbreaks[length(tcbreaks)], 'X+')))
      
    } else {
    p <- p +
      scale_x_continuous(breaks = tcbreaks, 
                         labels = c(tcbreaks[1:(length(tcbreaks) - 1)], paste0(tcbreaks[length(tcbreaks)], '+')))
    }
  }
  
  if (!is.null(axisfs)) {
    p <- p +
      theme(axis.text  = element_text(size = axisfs * 0.8),
            axis.title = element_text(size = axisfs),
            strip.text = element_text(size = axisfs))
  }
  
  if (!is.null(filename)) {
    if (!is.null(dim)) {
      ggsave( paste0(PLOT_PATH, filename, '.png'), plot = p,
              width = dim[1],
              height = dim[2])
    } else {
      ggsave( paste0(PLOT_PATH, filename, '.png'), plot = p)
    }
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
              subtitle = 'With all covariates',
              filename = 'ate_range_covar_240220',
              binw = 0.01) 


### With "best" covariate
ate_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              filename   = 'ate_range_covar_paper_240513',
              subtitle = 'With single best covariate',
              binw = 0.01 )



### Without covariates ----
ate_summary %>%
  filter(Covariates == "none") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Range of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'ate_range_240220',
              binw = 0.01 )





## 1b. Effect estimates by estimand ----

# calculate proportion missing--some estimands have no defined
# estimators for some of the outcomes
combined %>%
  group_by(UniqueID, Covariates, Estimand) %>%
  summarize(n = n(),
            N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            .groups = "drop") %>%
  mutate(PropMiss = Nmiss / (N + Nmiss)) %>%
  filter(PropMiss > 0) %>%
  arrange(desc(PropMiss))

# What estimands are in the dataset?
table( combined$Estimand, useNA = "always" )

ate_summary_estimand <- combined %>%
  # combining finite vs super population estimands 
  mutate(Estimand = str_replace_all(Estimand, 'SP\\*', 'SP')) %>%
  group_by(UniqueID, Collection, Covariates, Estimand) %>%
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




# Q: How much does the range go down by doing within estimand?

ate_summary
aa <- ate_summary_estimand %>%
  group_by( UniqueID, Covariates ) %>%
  summarise( rng = mean( Range ) )
aa
aaa <- ate_summary %>%
  left_join( aa, by=c( "UniqueID", "Covariates" ) ) %>%
  mutate( rngdiff =  rng / Range ) 

aaa %>%
  group_by(Covariates) %>%
  summarize(
    mean_range = mean(Range, na.rm = TRUE),
    min = min(rngdiff, na.rm = TRUE),
    max = max(rngdiff, na.rm = TRUE),
    iqr_value = IQR(rngdiff, na.rm = TRUE),
    mean = mean(rngdiff, na.rm = TRUE),
    sd = sd(rngdiff, na.rm = TRUE)
  )

ggplot( aaa, aes( x = Range, y = rng, col = Covariates ) ) +
  geom_point() +
  geom_abline( slope = 1, intercept = 0 ) +
  coord_fixed() +
  facet_wrap( ~ Covariates ) +
  theme_minimal() +
  labs( x = "Range of Point Estimates of β", y = "Mean Range within Estimand" )

# Commentary: The range of point estimates of β is reduced by around
# 40% on average when we look within estimand. But do both estimands
# have a lot of range within?


ate_summary
aa <- ate_summary_estimand %>%
  left_join( dplyr::select( ate_summary, UniqueID, Covariates, Range ), 
             by = c( "UniqueID", "Covariates" ) )
aa

ggplot( aa, aes( x = Range.y, y = Range.x, col = Estimand ) ) +
  facet_grid( Estimand ~ Covariates  ) +
  geom_point() +
  geom_abline( slope = 1, intercept = 0 ) +
  coord_fixed() +
  theme_minimal() +
  scale_x_log10() +
  scale_y_log10() +
  labs( x = "Range of Point Estimates of β", y = "Mean Range within Estimand" ) +
  theme(legend.position = "none")


aa %>%
  group_by(Estimand, Covariates) %>%
  mutate( rngdiff = Range.x / Range.y ) %>%
  summarize(
    min = min(rngdiff, na.rm = TRUE),
    max = max(rngdiff, na.rm = TRUE),
    iqr_value = IQR(rngdiff, na.rm = TRUE),
    mean = mean(rngdiff, na.rm = TRUE),
    sd = sd(rngdiff, na.rm = TRUE)
  )

# Commentary: The ranges seem similarly reduced within each estimand
# with the exception of no covariate adjustment, where the
# Person-Person has a much narrower range.
#
# It is very clear that covariate adjustment does increase,
# substantially, the range of estimates across the estimators.  This
# is unlike the Multisite case where covariate adjustment only had
# mild effects on differences between estimators.






# Research question 2 (the SE) ----
# To what extent can the choice of estimator of the standard error of the overall 
# average effect (SEˆ(β)) result in a different standard error? 
#   a. To what extent do the estimated standard errors vary across all estimators? 
#   b. For different types of estimands, to what extent do estimated standard errors 
#      vary across estimators that target that type of estimand? 


summary( combined$SE_hat_std )
summary( combined$df )
combined$df[ is.na( combined$df ) ] <- 0

# To make this workable we drop all estimates with a df of less than 5, and top-code all final ratios at 5

se_summary <- combined %>%
  dplyr::filter( df >= 5 ) %>%
  # replace values of Inf with NA
  mutate(SE_hat_std = ifelse(is.infinite(SE_hat_std), NA, SE_hat_std)) %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

write.xlsx(se_summary, plot_path( 'se_summary.xlsx') )


## 2a. Range of SE estimates across all estimators ----

### With covariates ----
se_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'With covariates (windsorized at 5)',
              filename = 'se_ratio_covar_240220')

se_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_paper_240513') 



### Without covariates ----
se_summary %>%
  filter(Covariates == "none") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'se_ratio_240220')

# making version without outlier
se_summary %>%
  filter(Covariates == "none") %>%
  filter(Ratio < 100) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of SE Estimates', 
              ylab     = 'Number of Outcomes', 
              title    = 'Ratio of Estimates across all Estimands', 
              subtitle = 'Without covariates', 
              filename = 'se_ratio_rmoutlier_240220')



## 2b. SE estimates by estimand ----

# Drop anything with small degrees of freedom
se_summary_estimand <- combined %>%
  dplyr::filter( df >= 5 ) %>%
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
  mutate(Ratio = pmin( 5, Max / Min)  )

write.xlsx(se_summary_estimand, plot_path( 'se_summary_estimand.xlsx') )

se_summary_estimand %>%
  filter(Covariates == "best") %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP(\\*|)_')) %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'se_ratio_covar_byestimand_paper_240513',
              byest    = T )

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


# How much variation in range and SE if we disallow the interacted estimators? ----

table( combined$Estimator )

cc_noFI <- combined %>%
  dplyr::filter( !grepl( '_FI_', Estimator ) )
table( cc_noFI$Estimator, cc_noFI$Estimand)


cc_noFI <- combined %>%
  dplyr::filter( !grepl( '_FI_', Estimator ) ) %>%
  group_by( UniqueID, Collection, Covariates, Estimand ) %>%
  summarise( n = n(),
             N     = sum(!is.na(SE_hat_std)),
             Nmiss = sum(is.na(SE_hat_std)),
             Min   = min(SE_hat_std, na.rm = T),
             Max   = max(SE_hat_std, na.rm = T),
             Ratio = Max / Min,
             Range = max(ATE_hat_std,na.rm=TRUE) - min(ATE_hat_std,na.rm=TRUE) )
filter( cc_noFI, n == Nmiss )

cc_noFI %>%
  filter(Covariates == "best") %>%
  # mutate(Estimand = str_remove_all(Estimand, 'SP(\\*|)_')) %>%
  makeDotplot(var      = 'Ratio',
              title    = "Ratio of SE for estimators without interactions",
              xlab     = 'Ratio of Largest to Smallest Estimates of SE(β)', 
              ylab     = 'Number of Outcomes', 
              filename = 'ratio_no_interacted_included',
              byest    = T )

cc_noFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Largest to Smallest Estimates of β',
              title    = "Range of point estimates for estimators without interactions",
              subtitle = "best covariate, no interacted estimators",
              ylab     = 'Number of Outcomes', 
              filename = 'range_no_interacted_included',
              byest    = T )


# Commentary:
# (1) Even without the interacted estimators, the standard errors can vary
# by 20% or more quite easily.  But this is much less of a big deal.
#
# (2) The range of point estimates still varies a lot for both estimands

cc_noFI
ate_summary_estimand
a_comp <- left_join( ate_summary_estimand, 
                     dplyr::select( cc_noFI, UniqueID, Covariates, Estimand, Ratio, Range ), 
                     by=c( "UniqueID", "Covariates", "Estimand" ) ) %>%
  mutate( range.red = Range.y / Range.x,
          SE.red = Ratio.y / Ratio.x )
a_comp






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

combined2 %>%
  group_by( Estimator ) %>%
  summarise( avg = mean( medresid_ate, na.rm = T ) ) %>%
  arrange( avg )


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
ggsave(plot_path( "ATE_median_byuniqueID_covar_240220.png" ),
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
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$|^DBi', Estimator)) %>%
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
              binw     = .01,
              dim      = c(9, 6))

## ATE range by estimand, best covariate
ate_summary_estimand_sree <- combined %>%
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$|^DBi', Estimator)) %>%
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
              xscale     = seq(0, .15, .025),
              dim        = c(9, 5),
              plotmargin = c(5, 20, 5, 20))

## SE min/max ratio, best covariate
se_summary_sree <- combined %>%
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$|^DBi', Estimator)) %>%
  filter(Collection == 'Curated') %>%
  group_by(Collection, Project, Trt_arm, Outcome, Covariates, J, K, ncov.1) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T),
            .groups = 'drop') %>%
  mutate(Ratio    = Max / Min,
         #top code ratio
         Ratio_tc = ifelse(Ratio >= 4, 4, Ratio),
         df       =  J - (2 * K) - ncov.1 - 1)

# no top-coding
se_summary_sree %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              # title    = 'Max/Min Ratio of SE Estimates across all Estimators',
              # subtitle = 'Each point represents one study-treatment arm-outcome combination',
              filename = 'se_ratio_bestcovar_sree2024',
              rmlegend = T,
              axisfs   = 20)

# top-coded ratio
se_summary_sree %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio_tc',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              # title    = 'Max/Min Ratio of SE Estimates across all Estimators',
              # subtitle = 'Each point represents one study-treatment arm-outcome combination',
              filename = 'se_ratio_tc_bestcovar_sree2024',
              rmlegend = T,
              axisfs   = 20,
              binw     = 0.1,
              coordratio = 1,
              tcratio  = T,
              tcbreaks = c(1:4),
              dim      = c(9, 5))

## SE max/min ratio by estimand, best covariate
se_summary_estimand_sree <- combined %>%
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$|^DBi', Estimator)) %>%
  filter(Collection == 'Curated') %>%
  group_by(Collection, Project, Trt_arm, Outcome, Covariates, Estimand) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T),
            .groups = 'drop') %>%
  mutate(Ratio = Max / Min,
         Ratio_tc = ifelse(Ratio >= 4, 4, Ratio))

# no top-coding
se_summary_estimand_sree %>%
  filter(Covariates == "best") %>%
  mutate(Estimand = case_when(Estimand == 'SP_Cluster_Cluster' ~ 'Cluster Weighted Estimand',
                              Estimand == 'SP_Person_Person'   ~ 'Person Weighted Estimand')) %>%
  makeDotplot(var        = 'Ratio',
              xlab       = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab       = 'Count',
              # title      = 'Ratio of SE Estimates across all Estimators, Grouped by Estimand',
              # subtitle   = 'Each point represents one study-treatment arm-outcome combination',
              filename   = 'se_ratio_bestcovar_estimand_sree2024',
              byest      = T,  
              rmlegend   = T,
              axisfs     = 20)

# with top-coding
se_summary_estimand_sree %>%
  filter(Covariates == "best") %>%
  mutate(Estimand = case_when(Estimand == 'SP_Cluster_Cluster' ~ 'Cluster Weighted Estimand',
                              Estimand == 'SP_Person_Person'   ~ 'Person Weighted Estimand')) %>%
  makeDotplot(var        = 'Ratio_tc',
              xlab       = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab       = 'Count',
              # title      = 'Ratio of SE Estimates across all Estimators, Grouped by Estimand',
              # subtitle   = 'Each point represents one study-treatment arm-outcome combination',
              filename   = 'se_ratio_tc_bestcovar_estimand_sree2024',
              byest      = T,  
              rmlegend   = T,
              axisfs     = 20,
              binw       = 0.1,
              coordratio = 1.5,
              tcratio    = T, 
              tcbreaks   = c(1:4),
              dim        = c(9, 5))

## DF vs range, best vs. all covar
(ate_summary_sree %>%
    filter(Covariates %in% c('best', 'all')) %>%
    # filter(df < 100) %>%
    # topcode high DF 
    mutate(Covariates = ifelse(Covariates == 'best', 'Best covariate', 'All covariates'),
           df_tc = ifelse(df >= 50, 50, df)) %>%
    ggplot(aes(x = df_tc, y = Range)) +
    geom_point() + 
    facet_wrap(~ Covariates) +
    labs(x = "Degrees of Freedom", y = "Range") +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgray", 
                                          linewidth = 0),
          axis.title = element_text(size = 20),
          axis.text  = element_text(size = 16),
          strip.text = element_text(size = 20),
          panel.spacing = unit(2, "lines")) +
    scale_x_continuous(breaks = seq(-10, 50, by = 10), 
                       labels = c(seq(-10, 40, by = 10), "50+"))) %>%
  ggsave(plot     = ., 
         filename = plot_path('dfvsrange_sree2024.png'),
         width = 14)

## DF vs SE ratio, best covar ----
se_summary_sree %>%
  # filter(df < 100) %>%
  # topcode high DF 
  mutate(df_tc = ifelse(df >= 50, 50, df)) %>%
  ggplot(aes(x = df_tc, y = Ratio_tc)) +
  geom_smooth( se=FALSE ) +
  # geom_point( size = 4 ) + 
  geom_jitter(height = 0,
              width = .5,
              size = 5) +
  labs(x = "Degrees of Freedom", y = 'Ratio of Largest to Smallest SE across Point Estimates') +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgray", 
                                        linewidth = 0),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        strip.text = element_text(size = 20),
        panel.spacing = unit(2, "lines")) +
  scale_x_continuous(breaks = seq(-10, 50, by = 10), 
                     labels = c(seq(-10, 40, by = 10), "50+")) +
  scale_y_continuous(breaks = seq(1, 4, by = 1), 
                     labels = c(paste0(seq(1, 3, by = 1), 'X'), "4X+")) 

ggsave(  filename = plot_path('dfvsratio_sree2024.png'),
         width  = 14,
         height = 7)


## ATEs and CIs for a "typical" outcome ----

combined_wt <- combined %>%
  filter(Covariates == 'best') %>%
  filter(Collection == 'Curated') %>%
  filter(!grepl('_(Block|FIRC|RIRC|HT|Raj)$|^DBi', Estimator)) %>%
  mutate(df_app =  J - (2 * K) - ncov.1 - 1,
         Weight = factor(case_when(Estimator %in% c("DB_FI_Person_Person",
                                                    "Agg_FI_Person_Person",
                                                    "LR_FI_CRVE_Person") ~ 'Person-weight',
                                   Estimator %in% c("LR_FE_CRVE",
                                                    "Agg_FE_Person",
                                                    "DB_FE_Person") ~ 'Precision-weight (person)',
                                   # group together MLMs?
                                   Estimator %in% c("MLM_FI_Cluster",
                                                    "MLM_FE",
                                                    "MLM_RE") ~ 'MLM',
                                   Estimator %in% c("Agg_FI_Cluster_Cluster",
                                                    "DB_FI_Cluster_Cluster") ~ 'Cluster-weight',
                                   Estimator %in% c("DB_FE_Cluster",
                                                    "Agg_FE_Cluster") ~ 'Precision-weight (cluster)'),
                         levels = c('Cluster-weight', 'Precision-weight (cluster)',
                                    'MLM',
                                    'Precision-weight (person)', 'Person-weight'))) %>%
  
  # order estimators
  mutate(Estimator = factor(Estimator, levels = c("Agg_FI_Cluster_Cluster", "DB_FI_Cluster_Cluster", # cluster-weight
                                                  "Agg_FE_Cluster", "DB_FE_Cluster",  # precision-weight (cluster)
                                                  "MLM_FE", "MLM_RE", "MLM_FI_Cluster", # MLM    
                                                  "Agg_FE_Person", "DB_FE_Person", "LR_FE_CRVE", # Precision-weight (person)
                                                  "Agg_FI_Person_Person", "DB_FI_Person_Person", "LR_FI_CRVE_Person")), # person-weight
         Estimand = factor(case_when(Estimand == 'SP_Cluster_Cluster' ~ 'Cluster Weighted Estimand',
                                     Estimand == 'SP_Person_Person'   ~ 'Person Weighted Estimand')))

library(grid)

for (uid in unique(combined_wt$UniqueID)) {
  
  temp <- combined_wt %>%
    filter(UniqueID == uid) %>%
    mutate(CImin = ATE_hat_std - 1.645*SE_hat_std,
           CImax = ATE_hat_std + 1.645*SE_hat_std)
  
  yseq <- round(seq(floor(min(temp$CImin, na.rm = T)*10)/10, 
                    ceiling(max(temp$CImax, na.rm = T)*10)/10,
                    .05), 2)
  
  minbeta <- min(temp$ATE_hat_std, na.rm = T)
  maxbeta <- max(temp$ATE_hat_std, na.rm = T)
  
  minse <- min(temp$SE_hat_std, na.rm = T)
  maxse <- max(temp$SE_hat_std, na.rm = T)
  
  tagtext <- bquote(atop(
    atop('Approx. DF = ' * .(max(temp$df_app)),
         hat(beta) * ' range = ' * .(round(maxbeta - minbeta, 2))),
    atop('Max-Min SE Ratio = ' * .(round(maxse/minse, 2)), 
         phantom())))
  
  
  
  plottitle <- paste0(unique(temp$Project),
                      ifelse(is.na(unique(temp$Trt_arm)), '', paste0(', Arm: ', unique(temp$Trt_arm))))
  plotsubtitle <- unique(temp$Outcome_desc)
  
  p <- temp %>%
    ggplot(aes(x=Estimator, y=ATE_hat_std )) + 
    geom_hline(yintercept = yseq, col = "gray", lty="dotted") +
    geom_errorbar(aes(ymin=(ATE_hat_std - 1.645*SE_hat_std), ymax=(ATE_hat_std + 1.645*SE_hat_std), col = Estimand), width=.2, lwd=1.1) +
    geom_point(size=2) +
    geom_hline(yintercept = min(temp$ATE_hat_std, na.rm = T), lty = 'dashed') +
    geom_hline(yintercept = max(temp$ATE_hat_std, na.rm = T), lty = 'dashed') +
    facet_grid(. ~ Weight , scales = "free") +
    scale_y_continuous(breaks= yseq) +
    scale_x_discrete(breaks = levels(combined_wt$Estimator),
                     labels = parse(text = paste0('bold(hat(beta)[', levels(combined_wt$Estimator), '])')))  +
    theme_bw() +
    theme(legend.position = "right",
          legend.title = element_text(size=20, face = "bold"),
          legend.text = element_text(size=16, face="bold"),
          plot.title = element_text(face="bold", size = 20),
          # legend.key.width = unit(1, "cm"),
          strip.text.x = element_text(angle = 0, size=10),
          axis.text.x = element_text(face= "bold",
                                     angle=90,
                                     size = 16,
                                     vjust = 0,
                                     hjust=0, 
                                     margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 35, r = 0, b = 0, l = 0), size = 20, face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size = 20, face="bold"),
          axis.text.y = element_text(face="bold",
                                     margin = margin(t=0,r=0,b=0,l=0), 
                                     size =16, 
                                     #  family = "serif",
                                     vjust=0)) +
    ylab(expression(paste(hat(beta)," (Effect Size Units)"))) + xlab("Estimator") +
    #  coord_flip() +
    guides(col=guide_legend(
      label.position = "right", nrow =4,
      #  keywidth=1,
      # keyheight=.25,
      #  default.unit="inch",
      shape = 50, 
      legend.key.width = unit(10, "cm") , size = 2) ) +
    ggtitle(plottitle,
            subtitle = plotsubtitle) + labs(color='Estimands') +
    # annotate('text', x = Inf, y = -0.1, label = 'TEST')
    labs(tag = tagtext) +
    theme(plot.tag.position = c(.8, .5),
          plot.tag = element_text(hjust = 0, size = 16))
  
  ggsave(plot_path(paste0('ATEwCI_', uid, '.pdf')),
         plot = p,
         width = 14,
         height = 8)
  
}

# combine into a single pdf
pdf_files <- list.files(path = here::here( "empirical/tables" ),
                        pattern = '^ATEwCI',
                        full.names = T)

library(qpdf)
pdf_combine(input = pdf_files,
            output = plot_path('ATEwCI_COMBINED.pdf'))
