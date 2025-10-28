

# Demo of checking specific prediction with the simulation results




# Load the workspace ----

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

source( here::here( "analysis/Simulation_Write_Up/load_and_prep_data.R" ) )


simdat = sruns.clean

# Check combos of simulations
table( cv.cluster = simdat$cv.cluster, 
       sd.cluster = simdat$sd.cluster,
       size.impact.cor.cluster = simdat$size.impact.cor.cluster )


table( simdat$variable.p, simdat$proptx.impact.cor )

table( simdat$size.impact.cor, simdat$sd.cluster )


simdat <- simdat %>%
  filter( !(size.impact.cor.cluster != 0 & sd.cluster == 0) ) %>%
  filter( !(size.impact.cor.cluster != 0 & cv.cluster == 0) ) 



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
                 sd.cluster != 0 ) %>%
  ggplot( aes( x = method, y = bias.s.estimand ) ) +
  facet_grid( size.impact.cor.cluster ~ cv.cluster, 
              labeller = label_both ) +
  geom_boxplot( outlier.size = 0.25 ) +
  coord_flip()

# The above plot confirms, at the cluster level, the bias.


# Now for proportion treated

dat %>%
  dplyr::filter( sd.cluster != 0, block_type != "Same" ) %>%
  ggplot( aes( x = method, y = bias.s.estimand ) ) +
  facet_grid( block_type + size.impact.cor.cluster ~ cv.cluster, 
              labeller = label_both ) +
  geom_boxplot( outlier.size = 0.25 ) +
  geom_hline( yintercept = 0, col="red" ) +
  coord_flip()


# Or this
DGP_vars
table( dat$ICC.block )
table( dat$sd.ATE )
table( dat$empirical )
table( dat$cv.block, dat$size.impact.cor )
table( dat$variable.p, dat$proptx.impact.cor )
table( dat$cv.nbar, dat$size.impact.cor )
table( dat$cv.nbar, dat$size.impact.cor.cluster )
table( dat$ICC.cluster )
table( dat$nbar.J.cor, dat$block_type )
table( dat$size.intercept.cor, dat$block_type )
table( dat$cv.nbar, dat$block_type )
table( dat$p.bar, dat$block_type )
table( dat$sd.ATE, dat$block_type )

dat %>%
  dplyr::filter( cv.cluster == 0,
                 
                 size.impact.cor.cluster == 0 ) %>%
  ggplot( aes( x = method, y = bias.s.estimand ) ) +
  facet_grid(  block_type ~ p.bar + sd.cluster, 
               labeller = label_both ) +
  geom_boxplot( outlier.size = 0.25 ) +
  geom_hline( yintercept = 0, col="red" ) +
  coord_flip()

dd <- dat %>%
  dplyr::filter( method == "LRapw-FE-db",
                 block_type == "Mess",
                 cv.cluster == 0,
                 size.impact.cor.cluster == 0 )
summary( dd$bias.s.estimand )

source( here::here( "analysis/Simulation_Write_Up/create_analysis_tree.R" ) )

DGP_vars

table( simdat$n.bar )
DGP_vars = setdiff(DGP_vars, c( "empirical", "R", "K", "n.bar",
                                "ICC.block", "num.cov.block" ) )


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
# Could it be the random chance of the randomly generated blocks
# for each of the scenarios?
#
#
# Next step: Rerun sim with a 0 bias and 0.20 bias to see if it
# replicates?  is this somehow monte carlo error???
#
# Or rerun for all scenarios for just the estimator given?



