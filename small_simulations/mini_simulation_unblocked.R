


# Simple simulation as sanity check of SEs etc. This simulation
# generates code for a single district CRT (i.e., we just have
# clusters).
#
# We can thus use the blkvar package and then re-assign tx with
# randomizr, so we can have clusters of unequal size, etc.

#
# NOTE: No districts means FEWER ESTIMATORS in the list from
# compare_methods




set.seed( 1039 )

library( tidyverse )
library( PUMP )
# devtools::install_github( "https://github.com/lmiratrix/blkvar" )

library( clusterRCT )



one_run <- function() {

    sim.data <- blkvar::generate_multilevel_data( n.bar = 50, J = 30, p = 0.3,
                                          size.impact.correlate = 1,
                                          variable.n = TRUE )
    sim.data <- mutate( sim.data,
                        Z = randomizr::cluster_ra(clusters = sid),
                        Yobs = ifelse( Z == 1, Y1, Y0 ) )

    c1 <- clusterRCT::compare_methods( Yobs ~ Z | sid,
                                       data=sim.data,
                                       include_method_characteristics = FALSE )

    c1$tau_indiv = mean( sim.data$Y1 - sim.data$Y0 )
    c1$tau_clust = sim.data %>%
        group_by(sid) %>%
        summarize(tau = mean(Y1-Y0)) %>%
        summarize(tau = mean(tau)) %>%
        pull(tau)

    c1
}

one_run()


rps = map_df( 1:100, ~ one_run(), .id = "runID", .progress=TRUE )

head( rps )


res <- rps %>% group_by( method ) %>%
    summarise( sdtau_indiv = sd( tau_indiv ),
               tau_indiv = mean( tau_indiv ),
               sdtau_clust = sd( tau_clust ),
               tau_clust = mean( tau_clust ),
               EATE = mean( ATE_hat ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ),
               bias_indiv = EATE - tau_indiv,
               bias_clust = EATE - tau_clust,
               SEb = sd(ATE_hat - tau_indiv) / sqrt(n()),
               bI_t = bias_indiv / SEb,
               bI_c = bias_clust / ( sd(ATE_hat - tau_clust) / sqrt(n()) ),
               biased = round( pmin( abs(bI_t), abs(bI_c) ) ) ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) ) %>%
    relocate( method, tau_indiv, sdtau_indiv, tau_clust )


res %>%
    dplyr::select( -sdtau_indiv )

res %>%
    dplyr::select( method, EATE, SE, bias_indiv, bias_clust, calib, biased )

# CONCLUSIONS: Estimators are, for the most part, either targeting
# individual or cluster average impacts, and the standard errors all
# seem reasonably calibrated under a superpopulation model.

# RESULTS:
# A tibble: 6 × 9
# method         sdtau_indiv tau_indiv sdtau_clust tau_clust  EATE    SE ESE_hat calib
# <glue>               <dbl>     <dbl>       <dbl>     <dbl> <dbl> <dbl>   <dbl> <dbl>
# 1 Agg-FE-cluster       0.109     0.338       0.102     0.184 0.190 0.339   0.314 0.927
# 2 Agg-FE-person        0.109     0.338       0.102     0.184 0.355 0.378   0.347 0.917
# 3 DB_clust             0.109     0.338       0.102     0.184 0.190 0.339   0.314 0.927
# 4 DB_indiv             0.109     0.338       0.102     0.184 0.355 0.378   0.341 0.902
# 5 LR w/ CRSE           0.109     0.338       0.102     0.184 0.355 0.378   0.347 0.917
# 6 MLM-NoFE             0.109     0.338       0.102     0.184 0.193 0.339   0.314 0.927

