

# Demo of checking specific prediction with the simulation results




# Load the workspace ----


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


select <- dplyr::select


my_t = theme_tufte() + 
  theme( legend.position="bottom", 
         legend.direction="horizontal",
         legend.key.width=unit(1,"cm"),
         panel.border = element_blank(), 
         axis.title=element_text(size=12))
theme_set( my_t )

DGP_vars <- c("K", "J.bar", "p.bar", "sd.ATE", "ICC", "variable.p", "n.bar", "sd.school", "ICC.school", 
              "size.impact.cor", "num.cov.district","cv.nbar", "R2.school", 
              "num.cov.school", "cv.school", "num.cov.indiv", 
              "empirical", "cv.district", "nbar.J.cor", "size.intercept.cor")

source( here::here( "analysis/Simulation_Write_Up/load_and_prep_data.R" ) )


simdat = sruns.clean

# Check combos of simulations
table( cv.school = simdat$cv.school, 
       sd.school = simdat$sd.school,
       size.impact.cor.school = simdat$size.impact.cor.school )


table( simdat$variable.p, simdat$proptx.impact.cor )

table( simdat$size.impact.cor, simdat$sd.school )


simdat <- simdat %>%
  filter( !(size.impact.cor.school != 0 & sd.school == 0) ) %>%
  filter( !(size.impact.cor.school != 0 & cv.school == 0) ) 



# Filter to possible ones
dat 




# Check 1-- From Mike ----



# From Mike Memo on June 3, 2025


# For example, I expect LRi-FE and LRapw-FE to be biased when there’s a
# correlation between B_k and N_k * p_k * (1 – p_k).

# -	Can we confirm this is the case under DGPs that hold all else constant?  

#    -	I think that we should expect zero or minimal bias under all other DGPs.  Is that correct?

DGP_vars

dat = simdat %>%
  filter( method %in% core.methods ) 
#%>%
#  mutate( bias_pot = size.impact.cor != 0 | proptx.impact.cor != 0 )

dat %>%
  dplyr::filter( block_type=="Same", 
                 sd.school != 0 ) %>%
  ggplot( aes( x = method, y = bias.s.estimand ) ) +
  facet_grid( size.impact.cor.school ~ cv.school, 
              labeller = label_both ) +
  geom_boxplot( outlier.size = 0.25 ) +
  coord_flip()

# The above plot confirms, at the school level, the bias.


# Now for proportion treated

dat %>%
  dplyr::filter( sd.school != 0, block_type != "Same" ) %>%
  ggplot( aes( x = method, y = bias.s.estimand ) ) +
  facet_grid( block_type + size.impact.cor.school ~ cv.school, 
              labeller = label_both ) +
  geom_boxplot( outlier.size = 0.25 ) +
  geom_hline( yintercept = 0, col="red" ) +
  coord_flip()


# Or this
DGP_vars
table( dat$ICC.district )
table( dat$sd.ATE )
table( dat$empirical )
table( dat$cv.district, dat$size.impact.cor )
table( dat$variable.p, dat$proptx.impact.cor )
table( dat$cv.nbar, dat$size.impact.cor )
table( dat$cv.nbar, dat$size.impact.cor.school )
table( dat$ICC.school )
table( dat$nbar.J.cor, dat$block_type )
table( dat$size.intercept.cor, dat$block_type )
table( dat$cv.nbar, dat$block_type )
table( dat$p.bar, dat$block_type )
table( dat$sd.ATE, dat$block_type )

dat %>%
  dplyr::filter( cv.school == 0,
                 
                 size.impact.cor.school == 0 ) %>%
  ggplot( aes( x = method, y = bias.s.estimand ) ) +
  facet_grid(  block_type ~ p.bar + sd.school, 
               labeller = label_both ) +
  geom_boxplot( outlier.size = 0.25 ) +
  geom_hline( yintercept = 0, col="red" ) +
  coord_flip()

dd <- dat %>%
  dplyr::filter( method == "LRapw-FE-db",
                 block_type == "Mess",
                 cv.school == 0,
                 size.impact.cor.school == 0 )
summary( dd$bias.s.estimand )

source( here::here( "analysis/Simulation_Write_Up/create_analysis_tree.R" ) )

DGP_vars

table( simdat$n.bar )
DGP_vars = setdiff(DGP_vars, c( "empirical", "R", "K", "n.bar",
                                "ICC.district", "num.cov.district" ) )


aa <- arrange( dd, -bias.s.estimand ) %>%
  dplyr::select( all_of( c( "bias.s.estimand", DGP_vars ) ) ) %>%
  select(where(~ n_distinct(.) > 1)) %>%
  mutate( bias.s.estimand = round( bias.s.estimand, 3 ) )

aa[-1] %>% summarise(across(everything(), n_distinct))
3 * 2 * 2 * 2 * 3

create_analysis_tree( aa,
                      "bias.s.estimand",
                      predictor_vars = colnames(aa)[-1],
                      num_leaves = 15, palette_name = "Greens", 
                      plot_tree = TRUE, initial_cp = 0.001,
                      tree_title = NULL)

knitr::kable( aa )
# view( aa )

# COMMENTS: The above table makes no sense.  Why are some things
# biased up, some down, and some about not?
#
# Could it be the random chance of the randomly generated districts
# for each of the scenarios?
#
#
# Next step: Rerun sim with a 0 bias and 0.20 bias to see if it
# replicates?  is this somehow monte carlo error???
#
# Or rerun for all scenarios for just the estimator given?



