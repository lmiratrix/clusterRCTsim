
# Comparing the fingerprint method and the rerandomization method for
# calibrated simulation
#
# This runs both on fakeCRT to see if the results agree


source( here::here( "R/calibrated_simulation_code.R" ) )

data( fakeCRT )
head( fakeCRT )

set.seed( 505054 )
R = 300


# Fingerprint ----

data( fakeCRT )
dd = describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data = fakeCRT,
                          control_formula = ~ X.jk + C.ijk )
dd

struc <- get_structure( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
struc

sim_res <- run_fingerprint_simulation( structure = struc, 
                                       parameters = dd,
                                       R = R,
                                       include_MLM = FALSE )
sim_res

ggplot( sim_res, aes( method, SE ) ) +
  geom_point( aes( y = ESEhat, color = "SE_hat") ) +
  geom_errorbar( aes( ymin = Q10, ymax = Q90, color = "SE_hat" ), width = 0 ) +
  geom_point( aes( color = "true SE" ), size = 2 ) +
  coord_flip()

sL <- sim_res %>%
  pivot_longer( cols = c(EATE, ESEhat, sdSEhat, SE), 
                names_to="metric", values_to="value" )
sL

ggplot( sL, aes( method, value, pch=weight )) +
  facet_wrap( ~ metric, scales="free_x" ) +
  geom_hline( yintercept = 0 ) +
  geom_point() +
  coord_flip()

sim_res


# Rerandomization ----


sim_res_rr <- run_rerandomize_simulation( Yobs ~ T.x | S.id | D.id, 
                                       data=fakeCRT,
                                       R = R,
                                       include_MLM = FALSE )
sim_res_rr

ggplot( sim_res_rr, aes( method, SE ) ) +
  geom_point( aes( y = ESEhat, color = "SE_hat") ) +
  geom_errorbar( aes( ymin = Q10, ymax = Q90, color = "SE_hat" ), width = 0 ) +
  geom_point( aes( color = "true SE" ), size = 2 ) +
  coord_flip()

sL <- sim_res_rr %>%
  pivot_longer( cols = c(EATE, ESEhat, sdSEhat, SE), names_to="metric", values_to="value" )
sL

ggplot( sL, aes( method, value, pch=weight )) +
  facet_wrap( ~ metric, scales="free_x" ) +
  geom_hline( yintercept = 0 ) +
  geom_point() +
  coord_flip()

sim_res_rr



# Compare the two simulation result sets ----

sim_res
sim_res_rr

all_sims = left_join( select( sim_res, -Q10, -Q90 ), 
                      select( sim_res_rr, method, SE, EATE, ESEhat, sdSEhat ), 
                      by="method",
                      suffix=c( "_f", "_r" ) )
all_sims <- all_sims %>%
  relocate( EATE_f, EATE_r,
            SE_f, SE_r,
            ESEhat_f, ESEhat_r,
            sdSEhat_f, sdSEhat_r,
            .after = method )
all_sims

ratios <- all_sims %>%
  mutate( ESE_rat = ESEhat_f / ESEhat_r,
          SE_rat = SE_f / SE_r,
          delta_ATE = EATE_f - EATE_r,
          cov_SEhat = sdSEhat_f / SE_f,
          sdSE_rat = sdSEhat_f / sdSEhat_r ) %>%
  dplyr::select( method, delta_ATE, SE_rat, ESE_rat, cov_SEhat, sdSE_rat )
ratios

# TAKE AWAY: It seems like fingerprint has more randomness in it, due
# to the regenerating the block effects each iteration.
# Rerandomization keeps those fixed as we permute within block.

# Interestingly, this makes the sd of the SE MUCH higher for the
# fingerprint vs. the rerandomization.  (Only 6% of the true SE).  Why
# is that?  I have no idea? (That said, the sd of SE-hat is small?)



