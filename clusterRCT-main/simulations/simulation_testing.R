

# Proof-of-concept simulation:
#  1. Generate district params
#  2. For each district:
#      - blkvar -> reassign tx with randomizr

# district params
#  - number of clusters
#  - proportion of clusters treated
#  - district-level average intercept
#  - district-level average treatment effect
#  - (and correlations of the above)


set.seed( 1039 )
library( tidyverse )
library( blkvar )
library( clusterRCT )


# generate correlated district-level ATEs & num. clusters / prop. tx clusters
#  - in Online Supplement B 1.1: sort noised random effects, order site sizes
gen_corr_districts <- function(K, eps=0.05) {
    # idea:
    #  1. sort (noised) random district effects
    #  2. order num. clusters & prop. tx clusters

    # values for J
    num_clusters <- rnorm(K, mean=15, sd=2) %>% round() %>% sort()
    # values for proportion treated
    prop_tx_clusters <- runif(K, min=0.25, max=0.75) %>% sort()

    # district-level parameters
    grid <- tibble(
        gamma.10 = rnorm(K, mean=0.2, sd=0.1),
    ) %>%
        arrange(gamma.10 + rnorm(K, mean=0, sd=eps)) %>%
        mutate(did = 1:n(), .before=everything(),
               num_clusters = num_clusters,
               prop_tx_clusters = prop_tx_clusters)

    # make student-level dataset
    sim.data <- grid %>%
        rowwise() %>%
        mutate(sim_data = list(
            generate_multilevel_data( n.bar = 50, J = num_clusters, p = 0.3,
                                      gamma.10 = gamma.10,
                                      size.impact.correlate = 1,
                                      variable.n = TRUE )))
    # rerandomize treatment
    sim.data <- sim.data %>%
        rowwise() %>%
        mutate(sim_data = list(
            sim_data %>%
                mutate(
                    Z = randomizr::cluster_ra(clusters = sid, prob = prop_tx_clusters),
                    Yobs = ifelse( Z == 1, Y1, Y0 ) )
        ))

    sim.data <- sim.data %>%
        unnest(cols=sim_data)


    # ensure that sids are unique
    sid_key <- sim.data %>%
        select(did, sid) %>%
        unique() %>%
        mutate(sid_new = 1:n())
    sim.data %>%
        left_join(sid_key, by=c("did", "sid")) %>%
        mutate(sid = sid_new) %>%
        select(-sid_new)
}

# check correlations
#  - note: they vary a lot, and can't be cranked to 1
#    - J will have repeats
#    - prop_tx is randomized by randomizr
if (F) {
    gen_corr_districts(K=20, eps=0.0) %>%
        group_by(did, sid) %>%
        summarize(ate = mean(Y1-Y0),
                  tx = mean(Z)) %>%
        summarize(ate = mean(ate),
                  J = length(unique(sid)),
                  prop_tx = mean(tx)) %>%
        summarize(cor_ate_J = cor(ate, J, method="spearman"),
                  cor_ate_tx = cor(ate, prop_tx, method="spearman"))
        # pivot_longer(c(J, prop_tx)) %>%
        # ggplot(aes(x=ate, y=value, color=name)) +
        # geom_point() +
        # facet_wrap(~name)
}

# check that sd(Y0) = 1
if (F) {
    # note: from blkvar, sd(Y0) = 0.97ish, 0.95ish with size.impact.correlate = 1
    benchmark <- map_dbl(1:100,
                         function(x) {
                             generate_multilevel_data( n.bar = 50, J = 25, p = 0.3,
                                                       gamma.10 = 0.2,
                                                       size.impact.correlate = 1,
                                                       variable.n = TRUE ) %>%
                                 # mutate(
                                 #     Z = randomizr::cluster_ra(clusters = sid, prob = 0.5),
                                 #     Yobs = ifelse( Z == 1, Y1, Y0 ) ) %>%
                                 summarize(sd = sd(Y0)) %>%
                                 pull(sd)
                         },
                         .progress = T)
    hist(benchmark)
    mean(benchmark)

    # so gen_corr_districts seems to preserve this property!
    bar <- map_dbl(1:100,
                   ~gen_corr_districts(K=10, eps=0.01) %>%
                       summarize(sd = sd(Y0)) %>%
                       pull(sd),
                   .progress=T)
    hist(bar)
    mean(bar)
}


one_run <- function() {

    # super hacky, just want to ensure that there are at least
    #  2 treated and 2 control clusters so that design-based
    #  variance estimators don't break!
    repeat {
        sim.data <- gen_corr_districts(K=5)
        min_txco <- sim.data %>%
            group_by(did, sid) %>%
            summarize(Z = mean(Z), .groups="drop_last") %>%
            summarize(ntx = sum(Z),
                      nco = sum(1-Z)) %>%
            select(-did) %>%
            min()
        if (min_txco > 1) {
            break
        }
    }

    c1 <- clusterRCT::compare_methods( Yobs ~ Z | sid | did, data=sim.data,
                                       include_method_characteristics = FALSE )

    c1 %>%
        mutate(
            tau_indiv = mean(sim.data$Y1-sim.data$Y0),
            tau_clust = sim.data %>%
                group_by(sid) %>%
                summarize(tau = mean(Y1-Y0), .groups="drop_last") %>%
                summarize(tau = mean(tau)) %>%
                pull(tau)
        )
}

if ( FALSE ) {
    one_run()
}



rps = map_df( 1:100, ~ one_run(), .id = "runID", .progress=TRUE )
head( rps )


rps %>% group_by( method ) %>%
    summarise( tau_indiv = mean( tau_indiv ),
               tau_clust = mean( tau_clust ),
               EATE = mean( ATE_hat ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ) ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) ) %>%
    knitr::kable(digits=2)

# RESULTS:
# # A tibble: 8 x 7
#   method         tau_indiv tau_clust  EATE    SE ESE_hat calib
#   <chr>              <dbl>     <dbl> <dbl> <dbl>   <dbl> <dbl>
# 1 DB_clust (FE)      0.351     0.213 0.203 0.206   0.205 0.995
# 2 DB_clust (int)     0.351     0.213 0.205 0.205   0.205 0.997
# 3 DB_indiv (FE)      0.351     0.213 0.329 0.221   0.219 0.989
# 4 DB_indiv (int)     0.351     0.213 0.329 0.219   0.209 0.956
# 5 LR (agg)           0.351     0.213 0.203 0.206   0.205 0.998
# 6 LR (agg, wt)       0.351     0.213 0.329 0.221   0.222 1.00
# 7 LR w/ CRSE         0.351     0.213 0.329 0.221   0.222 1.00
# 8 MLM                0.351     0.213 0.206 0.205   0.204 0.993




