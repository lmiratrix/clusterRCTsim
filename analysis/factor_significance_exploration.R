


knitr::opts_chunk$set(
  echo = FALSE,          # Don't show code in output by default
  message = FALSE,       # Suppress messages
  warning = FALSE,       # Suppress warnings
  fig.width = 7,         # Default figure width
  fig.height = 5,        # Default figure height
  fig.align = "center"   # Center figures
)


require(tidyverse)
require(rlang)
require(stargazer)
require( ggthemes )
require(broom)
require(knitr)
require(kableExtra)
require(patchwork)
require(gridExtra)
require(pander)
require(rpart)
require( rpart.plot )
require(RColorBrewer)
require(cowplot)

select <- dplyr::select


my_t = theme_tufte() + 
  theme( legend.position="bottom", 
         legend.direction="horizontal",
         legend.key.width=unit(1,"cm"),
         panel.border = element_blank(), 
         axis.title=element_text(size=12))
theme_set( my_t )

DGP_vars <- c("K", "J.bar", "p.bar", "sd.ATE", "ICC", "variable.p", "n.bar", "sd.cluster", "ICC.cluster", 
              "size.impact.cor", "num.cov.block","cv.nbar", "R2.cluster", 
              "num.cov.cluster", "cv.cluster", "num.cov.indiv", 
              "empirical", "cv.block", "nbar.J.cor", "size.intercept.cor")

source( here::here( "analysis/R/load_and_prep_data.R" ) )

source( here::here( "analysis/R/create_custom_boxplot.R" ) )

source( here::here( "analysis/R/create_analysis_tree.R" ) )


names( sruns.clean )
table( sruns.clean$group )

sruns = filter( sruns.clean, group=="full_factor")
nrow( sruns ) / n_distinct(sruns$ID) 
table( sruns$method )


# Do covariates impact bias? ----
table( sruns$num.cov.indiv )
table( sruns$num.cov.block )
table( sruns$num.cov.cluster )


sest <- sruns %>%
  filter(method %in% core.methods) %>% 
  mutate(method = str_remove(method, "-db$"),    
         method_ordered = str_remove(method_ordered, "-db$"))
table( sest$method )



# Prep the data ----

names( sest )
factors <- c( "block_type", "J.bar", 
              "ICC.cluster",
              "num.cov.cluster", "p.bar" )

factors <- c( "block_type", "J.bar", 
              "ICC.cluster",
              "num.cov.cluster", "p.bar" )

sest <- filter( sest,
                cv.cluster != 0,
                cv.nbar != 0 )

sest %>% 
  group_by(across (all_of( factors ) )) %>%
  summarise( n = n(), .groups = "drop" ) %>%
  pull( n )

# drop to a full factor experiment
names(sest)
table(sest$cv.block)
names( sest )
length( factors )
form = paste0( "bias.s.estimand.mag ~ (", paste( factors, collapse = " + " ), ")^5" )
form
M <- lm( form, data=sest )
coef(M)
length(coef(M))
head(sest)
form
anova_table <- aov( formula = as.formula(form),
                    data = sest )
knitr::kable( summary(anova_table)[[1]],
              digits = c(0,4,4,1,5) )

library(lsr)
eta <- etaSquared(anova_table, anova=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("source") %>%
  mutate( order = 1 + str_count(source, ":" ) ) %>%
  group_by( order ) %>%
  arrange( -eta.sq, .by_group = TRUE ) %>%
  relocate( order ) 

sum( eta$eta.sq.part, na.rm=TRUE )
sum( eta$eta.sq )

eta %>%
  knitr::kable( digits = 3 )

# NOTE: Unbalanced design causing this to not work?


