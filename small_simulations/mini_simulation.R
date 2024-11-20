


# Simple simulation as sanity check of SEs etc
#
# This sim has equal sized blocks and clusters and all that.



set.seed( 1039 )
library( tidyverse )
# library( PUMP )
library( clusterRCT )


ATE = 0.125
tx_het = TRUE

model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 30                          # number of schools
    , K = 10                          # number of districts
    , nbar = 10                       # number of individuals per school
    , S.id = NULL                     # N-length vector of school assignments
    , D.id = NULL                     # N-length vector of district assignments
    , Xi0 = 0                         # scalar grand mean outcome under no treatment
    , MDES = ATE            # minimum detectable effect size
    , R2.3 = 0.1              # percent of district variation
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.2*tx_het           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.0*tx_het          # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)

if ( FALSE ) {
    sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )
    class( sim.data )
    sd( sim.data$Yobs )
    mean( sim.data$Y1 - sim.data$Y0 ) / sd( sim.data$Yobs )
}



one_run <- function() {

    sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5, include_POs = TRUE )

    c1 <- clusterRCT::compare_methods( Yobs ~ T.x | S.id | D.id, data=sim.data,
                                       include_MLM = FALSE,
                                       include_LM = TRUE,
                                       include_agg = FALSE,
                                       include_method_characteristics = FALSE )

    c1$tau = mean( sim.data$Y1 - sim.data$Y0 )
    c1$sdY0 = sd( sim.data$Y0 )

    c1
}

if ( FALSE ) {
    # testing
    one_run()
}

R = 100
rps = map_df( 1:R, ~ one_run(), .id = "runID", .progress=TRUE )

head( rps )


res <- rps %>%
    group_by( method ) %>%
    summarise( Etau = mean( tau ),
               sdtau = sd( tau ),
               EATE = mean( ATE_hat ),
               SE = sd( ATE_hat ),
               ESE_hat = sqrt( mean( SE_hat^2 ) ),
               bias = mean( EATE - tau ),
               SEb = sd(ATE_hat-tau) / sqrt(n()),
               bias_t = bias / SEb ) %>%
    mutate( calib = sqrt( ESE_hat^2 / SE^2 ) )


res %>%
    arrange( calib ) %>%
    knitr::kable(digits=2)



