# packages
#install.packages("devtools","kfigr","ggplot2","ggthemes","grid","gridExtra","lattice",
 #               "dplyr","caret","knitr","kableExtra","tidyverse","pBrackets","extrafont","ggrepel")
library(devtools)
library(kfigr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(lattice)
library(dplyr)
library(caret)
library(knitr)
library(kableExtra)
library(tidyverse)
library(pBrackets)
library(extrafont)
library(ggrepel)
library(staplr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


######################## Read in restricted and nonrestricted adjusted results ###################################
if ( FALSE ) {
    nonres = read.csv("T:/SIESM/Cross-Study/Estimands_Estimators/Output/CovResultsNonrestricted191219.csv")
    res = read.csv("T:/SIESM/Cross-Study/Estimands_Estimators/Output/CovResultsrestricted191219.csv")
    covres1 = rbind(nonres,res)
    
    ### Remove RIRC  and Luke's extra estimators ###
    table(covres1$method, useNA = "if")
    covres2 = filter(covres1, !(method %in% c("RIRC-adj", "hybrid_m","hybrid_p","plug_in_big")))
    table(covres2$method, useNA = "if")
    write.csv(covres2, "T:/SIESM/Cross-Study/Estimands_Estimators/Output/CovResultsAll191220.csv")
}


covres2 = read_csv("CovResultsAll191220.csv" )

nrow( covres2 )
covres2 = filter( covres2, study != "TennStar" )
nrow( covres2 )

covres2$X1 = NULL

######################## Create new/ modify variables ###########################




# Reorder estimators #
covres2$model = factor(covres2$method, 
                      levels = c("DB-FP-Persons-adj" ,
                                 "FE-IPTW-adj",
                                 "FE-Int-Persons-adj",
                                 "DB-SP-Persons-adj",
                                 
                                 'FE-adj',
                                 'FE-Het-adj',
                                 'FE-CR-adj',
                                 'FE-Club-adj',
                                 
                                 'RICC-adj',
                                 'FIRC-adj',
                                 # 'RIRC',
                                 
                                 'DB-FP-Sites-adj',
                                 'FE-IPTW-Sites-adj',
                                 'FE-Int-Sites-adj',
                                 'DB-SP-Sites-adj'))



# Create Weight Category Variable ####
covres2$weight[covres2$method %in% c("DB-FP-Persons-adj" ,
                                       "FE-IPTW-adj",
                                       "FE-Int-Persons-adj",
                                       "DB-SP-Persons-adj" )] <- "N-weight"


covres2$weight[covres2$method %in% c('FE-adj',
                                       'FE-Het-adj',
                                       'FE-CR-adj',
                                       'FE-Club-adj')] <- "Precision-weight"


covres2$weight[covres2$method %in% c('RICC-adj',
                                       'FIRC-adj' #,
                                      # 'RIRC-adj' 
                                     )] <- "MLM"


covres2$weight[covres2$method %in% c('DB-FP-Sites-adj',
                                       'FE-IPTW-Sites-adj',
                                       'FE-Int-Sites-adj',
                                       'DB-SP-Sites-adj' )] <- "Equal-weight"


covres2$weight = factor(covres2$weight,
                       levels = c("N-weight","Precision-weight","MLM","Equal-weight"))


table(covres2$weight, useNA = "always")
table(covres2$method, useNA = "always")
table(covres2$method, covres2$weight, useNA = "always")
table(covres2$study, useNA = "always")
table(covres2$outcome, useNA = "always")


# Create Estimands Category Variable ####
covres2$estimand[covres2$method %in% c("DB-FP-Persons-adj" ,
                                       "FE-IPTW-adj",
                                       "FE-Int-Persons-adj",
                                       "FE-adj",
                                       "FE-Het-adj",
                                       "RICC-adj" )] <- "Person-Finite"

                                     
covres2$estimand[covres2$method %in% c('DB-SP-Persons-adj',
                                       'FE-CR-adj',
                                       'FE-Club-adj' )] <- "Person-Super"
                                      
                                     
covres2$estimand[covres2$method %in% c('DB-FP-Sites-adj',
                                       'FE-IPTW-Sites-adj',
                                       'FE-Int-Sites-adj' )] <- "Site-Finite"


covres2$estimand[covres2$method %in% c('FIRC-adj',
                                     'DB-SP-Sites-adj' )] <- "Site-Super"


covres2$estimand = factor(covres2$estimand,
                        levels = c("Person-Finite","Person-Super","Site-Finite","Site-Super"))


table(covres2$estimand, useNA = "always")
table(covres2$method, useNA = "always")
table(covres2$method, covres2$estimand, useNA = "always")
table(covres2$study, useNA = "always")
table(covres2$outcome, useNA = "always")


# For facetting#

covres2$finite2[covres2$estimand %in% c("Person-Finite", "Site-Finite" )] <- 1
covres2$finite2[! (covres2$estimand %in% c("Person-Finite", "Site-Finite" ))] <- 0
covres2$finite2 = factor(covres2$finite2, levels = c(0,1), labels = c('Superpop', 'Finite'))
table(covres2$finite2, useNA = "always")


covres2$finite3[covres2$estimand %in% c("Person-Finite" , "Site-Finite")] <- 1
covres2$finite3[! (covres2$estimand %in% c("Person-Finite" , "Site-Finite"))] <- 0
covres2$finite3 = factor(covres2$finite3, levels = c(1,0), labels = c('Finite', 'Superpop'))
table(covres2$finite3, useNA = "always")


covres2$person2[covres2$estimand %in% c("Person-Finite", "Person-Super" )] <- 1
covres2$person2[! (covres2$estimand %in% c("Person-Finite" , "Person-Super" ))] <- 0
covres2$person2 = factor(covres2$person2, levels = c(0,1), labels = c('Site', 'Person'))
table(covres2$person2, useNA = "always")


covres2$person3[covres2$estimand %in% c("Person-Finite", "Person-Super" )] <- 1
covres2$person3[! (covres2$estimand %in% c("Person-Finite", "Person-Super" ))] <- 0
covres2$person3 = factor(covres2$person3, levels = c(1,0), labels = c('Person', 'Site'))
table(covres2$person3, useNA = "always")

table(covres2$estimand,covres2$finite2,  useNA = "always")



# Label Outcomes #
outcome_labels = c(
  `R1TTARCR_wiz`        = 'Learning Communities - Targeted Credits, 1 Sem',
  `all_t_9_100_wiz`     = 'Early College High Schools - 9th Grade on Track',
  `c02mcret_wiz`        = 'Performance Based Scholarships - Total Credits, 2 Sem',
  `c06mcret_wiz`        = 'Performance Based Scholarships - Total Credits, 6 Sem',
  `c06mvdeg_wiz`        = 'Performance Based Scholarships - Degrees Earned, 6 Sem',
  `C3TOTCR_wiz`         = 'Learning Communities - Total Credits, 3 Sem',
  `C3TTARCR_wiz`        = 'Learning Communities - Targeted Credits, 3 Sem',
  `CHRONAB_100_wiz`     = 'Communities in Schools - Chronic Absenteeism',
  `earnt8q_2015_an_wiz` = 'Welfare to Work - Avg. Annual Earnings, qtrs 1-8',
  `eqjberna_2015_wiz`   = 'Career Academies - Avg. Annual Earnings, yrs 5-8',
  `EQJBMEMP_an_wiz`     = 'Career Academies - Avg. Months Worked, yrs 5-8',
  `FAIL_100_wiz`        = 'Communities in Schools - Fail', 
  `grad5yr_100_wiz`     = 'Early College High Schools - Earned Degree, 5 yrs',
  `lqjberna_2015_wiz`   = 'Career Academies - Avg. Annual Earnings, yrs 1-4', 
  `LQJBMEMP_an_wiz`     = 'Career Academies - Avg. months worked annually, yrs 1-4',
  `PQHSGRDGED100_wiz`   = 'Career Academies - Earned HS Diploma/GED, yr 5', 
  `PQPSEDDEG100_wiz`    = 'Career Academies - Enrolled in Post-secondary',
  `R1TOTCR_wiz`         = 'Learning Communities - Total Credits, 1 Sem',
  `tmathss_sepSD`       = 'Tennessee STAR - Math Scores',
  `treadss_sepSD`       = 'Tennessee STAR - Reading Scores',
  `R01BTVCRAT_wiz_e1`   = 'EASE Cohort 1 - Enrolled in the first summer after intervention',
  `R01BTMCRET_wiz_e1`   = 'EASE Cohort 1 - Credits earned in first summer after intervention',
  `R01BTVCRAT_wiz_e2`   = 'EASE Cohort 2 - Enrolled in the first summer after intervention',
  `R01BTMCRET_wiz_e2`   = 'EASE Cohort 2 - Credits earned in first summer after intervention' ,
  `ZGRCCSS_POST_wiz`    = 'Enhanced Reading Opportunities - Grade Reading Comprehension',
  `ZGRRVSS_POST_wiz`    = 'Enhanced Reading Opportunities - Grade Reading Vocabulary',
  `PCCRERCOREYR1_wiz`   = 'Enhanced Reading Opportunities - % of Required Credits Earned Y1',
  `PCCRERCOREYR2_wiz`   = 'Enhanced Reading Opportunities - % of Required Credits Earned Y2 ',
  `as_y1read_brz_wiz`   = 'After School Reading - SAT Total Reading',
  `as_y1math_brz_wiz`   = 'After School Math - SAT Total Math',
  `s03chkl_ext_wiz`     = 'Head Start Impact Study - Externalizing behavior problems',
  `s03irt_ppvtml_wiz`   = 'Head Start Impact Study - PPVT-III receptive vocabulary',
  `s03lei_selfr_wiz`    = 'Head Start Impact Study - Self-regulation skills',
  `s03wj3ap_w_wiz`      = 'Head Start Impact Study - WJ-III AP early numeracy',
  `s03wj3lw_w_wiz`      = 'Head Start Impact Study - WJ-III LW early reading',
  `s03wj3oc_w_wiz`      = 'Head Start Impact Study - WJ-III OC oral comprehension') 


covres2$B = covres2$tau

#### Make individaul graphs #####
if ( FALSE ) {
    
    #### Create a new folder with today's date for storing graphs ####
    orig.dir <- getwd()  ## I'm going to demo in a tempdir
    orig.dir
    #setwd("T:/SIESM/Cross-Study/Estimands_Estimators/Output")
    list.dirs()
    newdir <- paste("Graphs_Original", Sys.Date(), sep = "")
    dir.create(newdir)
    setwd(newdir)
    list.dirs()
    
#### loop for graphing each outcome ####
count = 0
for( i in outcome_labels){
  label = i
  count = count + 1
  varx = names(outcome_labels)[count]
  print(count)
  print(varx)
  print(label)
  
  temp.data = filter(covres2, outcome == varx)
  
  
  Plot <- ggplot(temp.data, aes(x=model, y=B )) + 
    geom_hline(yintercept = seq(-.20,.45,.05), col = "gray", lty="dotted") +
    geom_errorbar(aes(ymin=(B - 1.645*SE), ymax=(B + 1.645*SE), col = estimand), width=.2, lwd=1.1) +
    geom_point(size=2) +
    facet_grid(. ~ weight , scales = "free") +
    scale_y_continuous(breaks= c(-.20,-.15,-.10,-.05,0,.05,.10,.15,.20,.25,.30,.35,.40,.45),limits = c(-.20,.45)) +
    scale_x_discrete(breaks = c("DB-FP-Persons-adj" ,
                                "FE-IPTW-adj",
                                "FE-Int-Persons-adj",
                                "DB-SP-Persons-adj",
                                
                                'FE-adj',
                                'FE-Het-adj',
                                'FE-CR-adj',
                                'FE-Club-adj',
                                
                                'RICC-adj',
                                'FIRC-adj',
                                # 'RIRC',
                                
                                'DB-FP-Sites-adj',
                                'FE-IPTW-Sites-adj',
                                'FE-Int-Sites-adj',
                                'DB-SP-Sites-adj'), 
                     
                     labels = expression(bold(hat(beta)[DB-FP-Persons-adj]),
                                         bold(hat(beta)[FE-IPTW-adj]),
                                         bold(hat(beta)[FE-Int-Persons-adj]),
                                         bold(hat(beta)[DB-SP-Persons-adj]),
                                         
                                         bold(hat(beta)[FE-adj]),
                                         bold(hat(beta)[FE-Het-adj]),
                                         bold(hat(beta)[FE-CR-adj]),
                                         bold(hat(beta)[FE-Club-adj]),
                                         
                                         bold(hat(beta)[RICC-adj]),
                                         bold(hat(beta)[FIRC-adj]),
                               #          bold(hat(beta)[ML-RIRC]),
                                         
                                         bold(hat(beta)[DB-FP-Sites-adj]),
                                         bold(hat(beta)[FE-IPTW-Sites-adj]),
                                         bold(hat(beta)[FE-Int-Sites-adj]),
                                         bold(hat(beta)[DB-SP-Sites-adj])))  +
    theme_bw() +
    theme(legend.position = "right",
          legend.title = element_text(size=16, face = "bold"),
          legend.text = element_text(size=13, face="bold"),
          plot.title = element_text(face="bold", size = 18),
          # legend.key.width = unit(1, "cm"),
          strip.text.x = element_text(angle = 0, size=14),
          axis.text.x = element_text(face= "bold",
                                     angle=90,
                                     size = 16,
                                     vjust = 0,
                                     hjust=0, 
                                     margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 35, r = 0, b = 0, l = 0), size = 16, face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),size=16, face="bold"),
          axis.text.y = element_text(face="bold",
                                     margin = margin(t=0,r=0,b=0,l=0), 
                                     size =14, 
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
    ggtitle(label) + labs(color='Estimands')
  
  
  fname = paste0("plot_for_",varx,".pdf")
  ggsave(fname,plot = Plot, width = 10, height = 6)
  
}

# staple_pdf(input_directory = newdir )

}



#### Histograms for two main research questions ####



#### Covariate Adjusted Figure 3 in main paper ####
nrow( covres2 )
question1B = covres2 %>% 
    group_by(study, outcome, person3 ) %>%  #, N, S) %>% 
    summarise(range.B = max(B) - min(B),
              # range.SE = max(SE) - min(SE),
              ratio.SE = max(SE)/min(SE)) 

question1B
nrow( question1B )


p1b = ggplot(question1B, aes(x=range.B)) +
    geom_dotplot(fill="black", binwidth = .005, dotsize = .5, stackratio = .625)  +
    facet_grid( . ~ person3 ) + theme_gray() +
    scale_x_continuous(limits = c(0,.10), breaks = seq(0,.10,.025)) +
    scale_y_continuous(limits = c(0,2)) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 14, face = "plain"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),size = 14),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),size = 14),
          plot.title = element_text(size=16, face = "bold"),
          text = element_text(size=14, face = "bold"),
          panel.spacing.x = unit(1.5, "cm"),
          strip.text = element_text(size = 14, face = "bold")) +
    #    labs(title = "Range of Estimates within a Given Estimand\n") + ylab("Count") +
    labs( x="Range of estimates", y = "Number of outcomes") 
p1b
ggsave( "range_within_estimand_B_cov.pdf", height = 5.5, width = 14.8 )




#### Max to Min SE ratio with estimand plot (Figure 6) ####

question1 = covres2 %>% 
    group_by(study, outcome, estimand, finite3, person3 ) %>%  #, N, S) %>% 
    summarise(range.B = max(B) - min(B),
              # range.SE = max(SE) - min(SE),
              ratio.SE = max(SE)/min(SE)) 

question1
nrow( question1 )

p2 = ggplot(question1,aes(x=ratio.SE)) +
    geom_dotplot(fill="black", dotsize = .75, stackratio = .75, binwidth=0.05) +
    facet_grid(person3 ~ finite3) + theme_gray() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 14, face = "plain"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),size = 16),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),size = 16),
          plot.title = element_text(size=18, face = "bold"),
          text = element_text(size=14, face = "bold"),
          panel.spacing.x = unit(1.5, "cm"),
          strip.text = element_text(size = 16, face = "bold")) +
    scale_x_continuous( breaks = c(1.0, 1.5, 2.0 ) ) + #labs(title = "Ratio of SE within a Given Estimand\n")  + 
    ylab("Number Outcomes") + xlab("Ratio of SE Estimates")
p2
ggsave( "ratio_SE_within_estimand_cov.pdf", width=6, height=4  )


#### Other plots (not used) #####

## Range of scores across estimands
question2 = covres2 %>% 
  group_by(study, outcome) %>% 
  summarise(range.B = max(B) - min(B),
            range.SE = max(SE) - min(SE),
            ratio.SE = max(SE)/min(SE))


p3 = ggplot(question2,aes(x=range.B)) +
  geom_dotplot(fill="black", binwidth = .005, dotsize = .48)  +
  scale_x_continuous(limits = c(0,.10), breaks = seq(0,.10,.025)) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        plot.title = element_text(size=14, face = "bold"),
        text = element_text(size=12, face = "bold"),
        panel.spacing.x = unit(1.5, "cm"),
        strip.text = element_text(size = 13, face = "bold")) +
  labs(title = "Range of Estimates across all Estimands") + ylab("Count") + xlab("Range")

p3



p4 = ggplot(question2,aes(x=ratio.SE)) +
  geom_dotplot(fill="black", binwidth = .005, dotsize = 10)  +
  scale_x_continuous(limits = c(1,4.5), breaks = seq(1,4.5,.25)) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        plot.title = element_text(size=14, face = "bold"),
        text = element_text(size=12, face = "bold"),
        panel.spacing.x = unit(1.5, "cm"),
        strip.text = element_text(size = 13, face = "bold")) +
  labs(title = "Ratio of SE across all Estimands") + ylab("Count") + xlab("Ratio")
p4



# Save plots #
if ( FALSE ) {
    
file = newdir <- paste("T:/SIESM/Cross-Study/Estimands_Estimators/Output/Main Histograms Covariates", Sys.Date(),".pdf", sep = "")
pdf(file, height = 10, width=10)
p1
p2
p3
p4
dev.off()

}

setwd(orig.dir)
