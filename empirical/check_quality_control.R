

# Check our estimates against collection of reported effects from
# literature


source( here::here( "empirical/load_empirical_results.R" ) )

# Secure server data root; override via env var for a different mount point.
IES_CRT_ROOT <- Sys.getenv( "IES_CRT_DATA_ROOT", "/data/share/pse/IES CRT" )

## QC effect size estimates ----

# read in Excel file with reported effect sizes
report_es <- read_excel(file.path( IES_CRT_ROOT, "tables/Report vs Our Outcomes 240429.xlsx" )) %>%
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
  filter(Covariates == "best") %>%
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

ggsave(filename = file.path( IES_CRT_ROOT, "tables/QC/QC_EffectSize_240429.png" ),
       plot = plot_grid,
       width = 15,
       height = 15)
