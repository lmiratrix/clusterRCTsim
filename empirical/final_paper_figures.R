# load data ----

# load all data
source( here::here( "empirical/load_empirical_results.R" ) )

# filter to curated studies
combined <- combined %>%
  filter(Collection == 'Curated')

# helper functions and summary datasets ----
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
      theme(legend.position = 'none') +
      scale_fill_manual(values = 'black')
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

## ATE ----
# summarize ATE by study-outcome
ate_summary <- combined %>%
  group_by(UniqueID, Collection,Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(ATE_hat_std)),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = "drop") %>%
  mutate(Range = Max - Min)

# summarize ATE by study-outcome-estimand
ate_summary_estimand <- combined %>%
  group_by(Collection, Project, Trt_arm, Outcome, Covariates, Estimand) %>%
  summarize(N     = n(),
            Nmiss = sum(is.na(ATE_hat_std)),
            Min   = min(ATE_hat_std, na.rm = T),
            Max   = max(ATE_hat_std, na.rm = T),
            .groups = 'drop') %>%
  mutate(Range = Max - Min)

## SE ----

se_summary <- combined %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

se_summary_estimand <- combined %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

## SE (excluding fully interacted models) ----

se_summary_exclFI <- combined %>%
  filter(!grepl('-FI', Estimator)) %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

se_summary_estimand_exclFI <- combined %>%
  filter(!grepl('-FI', Estimator)) %>%
  group_by(Collection, Project, Trt_arm, Outcome, Outcome_desc, Covariates, Estimand) %>%
  summarize(N     = sum(!is.na(SE_hat_std)),
            Nmiss = sum(is.na(SE_hat_std)),
            Min   = min(SE_hat_std, na.rm = T),
            Max   = max(SE_hat_std, na.rm = T)) %>%
  mutate(Ratio = pmin( 5, Max / Min) )

# Figures ----
## Figure 1: The Range (across 11 Estimators) of Estimates of β for 28 Outcomes ----
ate_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β', 
              ylab     = 'Number of Outcomes', 
              title    = 'Figure 1', 
              filename   = 'final paper figures/fig1_ATErange',
              subtitle = 'The Range (across 11 Estimators) of Estimates of β for 28 Outcomes',
              binw = 0.01,
              coordratio = 0.1,
              rmlegend = T)

## Figure 2: The Range of Estimates of β for all Estimators Targeting a Given Estimand for 28 Outcomes ----
ate_summary_estimand %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Range',
              xlab     = 'Range of Point Estimates of β',
              ylab     = 'Number of Outcomes',
              filename = 'final paper figures/fig2_ATErange_byestimand',
              title = 'Figure 2',
              subtitle = 'The Range of Estimates of β for all Estimators\nTargeting a Given Estimand for 28 Outcomes',
              coordratio = .05,
              rmlegend = T,
              byest    = TRUE )

## Figure 3: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) for 28 outcomes. ----
se_summary %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 3',
              subtitle = 'The Ratio (across 19 estimators) of the Largest to Smallest\nEstimates of SE(β) for 28 outcomes',
              filename = 'final paper figures/fig3_SEratio',
              rmlegend = T,
              binw     = 0.1,
              coordratio = 1,
              tcratio  = T,
              tcbreaks = c(1:5))

## Figure 4: The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β) among all estimators targeting a given estimand for 28 outcomes ----
se_summary_estimand %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 4',
              subtitle = 'The Ratio (across 19 estimators) of the Largest to Smallest Estimates of SE(β)\namong all estimators targeting a given estimand for 28 outcomes',
              filename = 'final paper figures/fig4_SEratio_byestimand',
              byest = T,
              rmlegend = T,
              binw     = 0.1,
              coordratio = 1.5,
              tcratio  = T,
              tcbreaks = c(1:5))

## SE plots excluding fully interacted models ----

se_summary_exclFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 3a',
              subtitle = 'The Ratio (across 10 estimators) of the Largest to Smallest\nEstimates of SE(β) for 28 outcomes\n(excluding fully interacted models)',
              filename = 'fig3a_SEratio_exclFI',
              rmlegend = T,
              binw     = 0.1,
              coordratio = 1)

se_summary_estimand_exclFI %>%
  filter(Covariates == "best") %>%
  makeDotplot(var      = 'Ratio',
              xlab     = 'Ratio of Largest to Smallest SE across Point Estimates',
              ylab     = 'Count',
              title    = 'Figure 4a',
              subtitle = 'The Ratio (across 10 estimators) of the Largest to Smallest Estimates of SE(β)\namong all estimators targeting a given estimand for 28 outcomes\n(excluding fully interacted models)',
              filename = 'fig4a_SEratio_byestimand_exclFI',
              byest = T,
              rmlegend = T,
              binw     = 0.1,
              coordratio = 1.5)

