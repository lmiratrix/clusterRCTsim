#
# Design-based methods follow those found in Schochet JASA paper
# "Design-Based Ratio Estimators and Central Limit Theorems for
# Clustered, Blocked RCTs" Peter Z. Schochet, Nicole E. Pashley, Luke
# W. Miratrix, and Tim Kautz
#
#
# Also implements the Middleton and Aronow design based Raj estimator.



get_overall_ATE <- function( mod, weights ) {
    cc = coef(mod)

    J = length( weights )

    ids <- grep( "Z:", names(cc) )
    stopifnot(length(ids) == J)
    ATE_hats <- cc[ids]

    weighted.mean( ATE_hats, w = weights )
}



# If we ran an individual level regression, need to aggregate
# residuals to calculate Standard Errors.
agg_residuals <- function( data ) {

    data_agg <- data %>%
        group_by( blockID, clusterID, Z ) %>%
        summarise( n = n(),
                   .weight = sum( .weight ),
                   resid = mean( resid ) )
    data_agg
}


# Notes on Schochet paper notation
#
# m_b^1 = Total number of clusters assigned to tx group in block b





# Implement formula (9) and following equations from Schochet et al
# paper.
#
# See page 2139-40 (section 4.2)
#
# This return a dataframe, one row per block, of standard errors for
# impact estimates.
#
# @param adt Aggregated (cluster-level) data.
# @param v Degrees of freedom.  Use v=# covariates.
schochet_variance_formula_block <- function( adt, v, aggregated=TRUE ) {

    require( tidyr )
    require( dplyr )

    stopifnot( !is.null( adt$resid ) )
    stopifnot( !is.null( adt$Z ) )

    if ( !aggregated ) {
        adt = agg_residuals( adt )
    }

    # This calculates the sum part of the s2 calculations.
    calc_s2_partial <- function( .weight, resid ) {

        w_bar = mean( .weight )
        tot = sum( ( (.weight^2) / (w_bar^2) ) * resid^2 )

        tot
    }

    W = sum( adt$.weight )
    adtw <- adt %>%
        group_by( blockID, Z ) %>%
        summarize( m = n(),
                   n = sum( n ),
                   totwt = sum( .weight ),
                   s2 = calc_s2_partial( .weight, resid ),
                   .groups = "drop" ) %>%
        dplyr::group_by( blockID ) %>%
        dplyr::mutate( p = m / sum(m),
                       weight = sum( totwt ) ) %>%
        dplyr::select( -totwt ) %>%
        ungroup() %>%
        dplyr::mutate( q = weight / W ) %>%
        dplyr::group_by( Z ) %>%
        dplyr::mutate( sdf = (m - v * p * q - 1),
                       s2 = s2 / sdf ) %>%
        ungroup()

    # calc_s2_partial can give negative numbers on order of 10^-32
    adtw$s2[ adtw$s2 < 0 ] = 0

    adtw$s2[ (adtw$sdf <= 0) | is.infinite(adtw$s2) ] = NA

    pt2 <- adtw %>%
        tidyr::pivot_wider( names_from = "Z",
                            values_from = c( s2, sdf, m, n, p ) ) %>%
        dplyr::mutate( varD = s2_1 / m_1 + s2_0 / m_0,
                       SE_hat = sqrt( varD ),
                       n = n_0 + n_1,
                       m = m_0 + m_1 )


    pt2 <- pt2 %>%
        dplyr::select( -p_0 ) %>%
        rename( p = p_1 )

    pt2$weight = pt2$weight / sum(pt2$weight)   # normalize weights to avoid overflow

    pt2
}




# From Schochet et al equation (15) for the restricted model with no
# interaction term for tx and block.
schochet_FE_variance_formula <- function( adt, v, aggregated=TRUE ) {

    stopifnot( !is.null( adt$resid ) )
    stopifnot( !is.null( adt$Z ) )

    if ( !aggregated ) {
        adt = agg_residuals( adt )
    }

    m = nrow(adt)
    h = length( unique(adt$blockID) )

    SE_ATE <- adt %>%
        group_by( blockID ) %>%
        summarize( p_b = mean(Z),
                   m_b = n(),
                   wbar_b = mean(.weight),
                   num = sum(.weight^2 * (Z-p_b)^2 * resid^2)) %>%
        summarize( SE_ATE = sqrt( m / (m-h-v-1) * sum(num) /
                                      sum(m_b*p_b*(1-p_b)*wbar_b)^2 )) %>%
        pull(SE_ATE)

    df = m - v - h - 1

    list( SE_hat = SE_ATE, df = df )
}



#' Estimate ATE via design-based approach
#'
#' This follows "Design-Based Ratio Estimators and Central Limit
#' Theorems for Clustered, Blocked RCTs" by Peter Z. Schochet, Nicole
#' E. Pashley, Luke W. Miratrix, and Tim Kautz.
#'
#' It runs a weighted linear regression, as discussed in that paper.
#' The paper discusses two possabilities: running on the aggregate, or
#' running on the individual.
#'
#' In this implementation, we run the regression on the
#' cluster-aggregated data.
#'
#' @inheritParams linear_model_estimators
#' @param aggregated TRUE means data is already aggregated (and in
#'   canonical form).  FALSE means it is not.
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
design_based_estimators <- function( formula,
                                     data = NULL,
                                     control_formula = NULL,
                                     weight = c( "Person", "Cluster" ),
                                     aggregated = FALSE,
                                     include_block_estimates = FALSE ) {

    require( estimatr )

    # Determine which version of the estimator we are doing.
    weight <- match.arg(weight)

    if ( !aggregated ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula = control_formula )
        # Collapse to clusters
        data_agg = aggregate_data( data, control_formula )
        control_formula = attr( data_agg, "control_formula" )
    } else {
        data_agg = data
    }
    has_block = "blockID" %in% names( data )


    # Make our weights variable
    suff = ""
    if (weight == "Person") {
        data_agg$.weight <- data_agg$n
        suff = "Person"
    } else {
        # data_agg$.weight <- 1 / data_agg$n
        data_agg$.weight <- 1
        suff = "Cluster"
    }

    v = number_controls(control_formula)
    J = length( unique( data$clusterID ) )

    # block level stuff: only if we have a blockID do we fit an
    # interacted model or allow for fixed effects for each block.
    ATE_int = NA
    if ( has_block ) {
        K = length( unique( data$blockID ) )

        # Fully interacted model
        form = make_regression_formula( Yobs = "Ybar", interacted = TRUE,
                                        control_formula = control_formula )
        mod = lm( form, data=data_agg, weights = .weight )
        no_estimate = sum( is.na( coef(mod) ) ) > 0

        data_agg$resid = residuals(mod)
        block_tab = schochet_variance_formula_block( data_agg, v = v )

        # Translate Schochet notation to our notation.
        block_tab <- rename( block_tab,
                             J = m )

        ests = generate_all_interacted_estimates( mod, data = data_agg,
                                                  use_full_vcov = FALSE,
                                                  SE_table = block_tab,
                                                  method = glue::glue("DB_FI_{suff}"),
                                                  aggregated = TRUE,
                                                  include_block_estimates = include_block_estimates )
        # df = m - 2h - v* = #clusters - 2 #blocks - #covariates
        df = J - 2*K - v
        ests$df = df
        if ( no_estimate ) {
            ests$ATE_hat = NA
            ests$SE_hat = NA
            ests$p_value = NA
        }

        # Random debugging code -- here for future checking.
        if ( FALSE ) {
            aa <- aggregation_estimators( formula = NULL, data_agg, aggregated = TRUE )
            aa

            mm <- data_agg %>% group_by( blockID, Z ) %>%
                summarise( n = n(),
                           s2 = var( resid ),
                           s2Y = var( Ybar ), .groups = "drop" ) %>%
                pivot_wider( names_from=Z, values_from=c(n, s2, s2Y ) ) %>%
                mutate( varD = s2_1 / n_1 + s2_0 / n_0,
                        SE_hat = sqrt( varD ),
                        w = (n_0+n_1) / sum( n_0+n_1) )
            mm

            mm$SE_hat / block_tab$SE_hat
            sqrt( with( mm, sum( w^2 * varD ) ) )
            mm$J = block_tab$J
            mm$n = block_tab$n
            generate_all_interacted_estimates( mod, data = data_agg,
                                               use_full_vcov = FALSE,
                                               SE_table = mm,
                                               method = "tmp",
                                               aggregated = TRUE,
                                               include_block_estimates = TRUE )
        }


        # fixed effects model
        form = make_regression_formula( Yobs = "Ybar", FE = TRUE,
                                        control_formula = control_formula )
        mod_FE = lm( form, data=data_agg, weights = .weight )
        ATE_FE = coef(mod_FE)[ "Z" ]
        data_agg$resid = residuals(mod_FE)
        SE_FE = schochet_FE_variance_formula( data_agg, v = v )

        ests$weight = paste( suff, ests$weight, sep="_" )

        FE <- tibble(
            method = c( glue::glue( "DB_FE_{suff}" ) ),
            weight = suff,
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE_FE$SE_hat ),
            df = c( SE_FE$df ),
            p_value = two_sided_p( ATE_hat, SE_hat, df )
        )

        # Pack up results
        ests <- bind_rows( ests, FE )

        if ( include_block_estimates ) {
            blocks <- attr( ests, "blocks" )
            blocks = left_join( blocks, block_tab, by="blockID" )
            attr( ests, "blocks" ) <- blocks
        }

        # clean up NaNs
        ests$ATE_hat[ is.nan(ests$ATE_hat) ] = NA
        ests$SE_hat[ is.nan(ests$SE_hat) ] = NA
        ests$SE_hat[ ests$df < 1 ] = NA
        ests$p_value[ is.na(ests$SE_hat) ] = NA
        ests

    } else {
        # No Blocks, so we use simple Cluster randomized estimators.
        form = make_regression_formula( Yobs = "Ybar", FE = FALSE,
                                        control_formula = control_formula )
        mod = lm( form, data=data_agg, weights = .weight )
        ATE_FE = coef(mod)[ "Z" ]

        data_agg$blockID = "Ind"
        data_agg$resid = residuals( mod )
        SE = schochet_variance_formula_block( data_agg, v = v )

        df = nrow(data_agg) - 2 - v

        tibble(
            method = glue::glue( "DB_{suff}" ),
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE$SE_hat ),
            df = df,
            p_value = two_sided_p(ATE_hat, SE_hat, df )
        )
    }
}





#' Estimate ATE via design-based approach
#'
#' This function follows `design_based_estimators`, but runs the
#' regression on the individual level data.
#'
#' @inheritParams design_based_estimators
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
design_based_estimators_individual <- function( formula,
                                                data = NULL,
                                                control_formula = NULL,
                                                weight = c( "Person", "Cluster" ),
                                                include_block_estimates = FALSE ) {

    require( estimatr )

    # Determine which version of the estimator we are doing.
    weight <- match.arg(weight)

    has_block = "blockID" %in% names( data )

    # Make our weights variable
    suff = ""
    if (weight == "Person") {
        data$.weight <- 1
        suff = "Person"
    } else {
        data <- data %>%
            group_by( clusterID ) %>%
            mutate( .weight = 1 / n() ) %>%
            ungroup()
        suff = "Cluster"
    }

    v = number_controls(control_formula)
    J = length( unique( data$clusterID ) )

    # block level stuff: only if we have a blockID do we fit an
    # interacted model or allow for fixed effects for each block.
    ATE_int = NA
    if ( has_block ) {
        K = length( unique( data$blockID ) )

        # Fully interacted model
        form = make_regression_formula( interacted = TRUE,
                                        control_formula = control_formula )
        mod = lm( form, data=data, weights = .weight )
        no_estimate = sum( is.na( coef(mod) ) ) > 0

        data$resid = residuals(mod)
        block_tab = schochet_variance_formula_block( data, v = v, aggregated=FALSE )

        # Translate Schochet notation to our notation.
        block_tab <- rename( block_tab,
                             J = m )

        ests = generate_all_interacted_estimates( mod, data = data,
                                                  use_full_vcov = FALSE,
                                                  SE_table = block_tab,
                                                  method = glue::glue("DBi_FI_{suff}"),
                                                  aggregated = FALSE,
                                                  include_block_estimates = include_block_estimates )
        df = J - 2*K - v
        ests$df = df
        if ( no_estimate || df < 0 ) {
            ests$ATE_hat = NA
            ests$SE_hat = NA
            ests$p_value = NA
        }

        # fixed effects model
        form = make_regression_formula( FE = TRUE,
                                        control_formula = control_formula )
        mod_FE = lm( form, data=data, weights = .weight )
        ATE_FE = coef(mod_FE)[ "Z" ]
        data$resid = residuals(mod_FE)
        SE_FE = schochet_FE_variance_formula( data, v = v, aggregated=FALSE )

        ests$weight = paste( suff, ests$weight, sep="_" )

        FE <- tibble(
            method = c( glue::glue( "DBi_FE_{suff}" ) ),
            weight = suff,
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE_FE$SE_hat ),
            df = c( SE_FE$df ),
            p_value = two_sided_p( ATE_hat, SE_hat, df )
        )

        # Pack up results
        ests <- bind_rows( ests, FE )

        if ( include_block_estimates ) {
            blocks <- attr( ests, "blocks" )
            blocks = left_join( blocks, block_tab, by="blockID" )
            attr( ests, "blocks" ) <- blocks
        }

        # clean up NaNs
        ests$ATE_hat[ is.nan(ests$ATE_hat) ] = NA
        ests$SE_hat[ is.nan(ests$SE_hat) ] = NA
        ests$SE_hat[ ests$df < 1 ] = NA
        ests$p_value[ is.na(ests$SE_hat) ] = NA

        ests

    } else {
        # No Blocks, so we use simple Cluster randomized estimators.
        form = make_regression_formula( FE = FALSE,
                                        control_formula = control_formula )
        mod = lm( form, data=data, weights = .weight )
        ATE_FE = coef(mod)[ "Z" ]

        data$blockID = "Ind"
        data$resid = residuals( mod )
        SE = schochet_variance_formula_block( data, v = v, aggregated = FALSE )

        df = J - 2 - v

        tibble(
            method = glue::glue( "DBi_{suff}" ),
            ATE_hat = c( ATE_FE ),
            SE_hat = c( SE$SE_hat ),
            df = df,
            p_value = two_sided_p(ATE_hat, SE_hat, df )
        )
    }
}







#' Estimate ATE via Middleton & Aronow unbiased approach
#'
#' This follows Middleton, J. A., & Aronow, P. M. (2015). Unbiased
#' Estimation of the Average Treatment Effect in Cluster-Randomized
#' Experiments. 6(1–2), 312. doi: 10.1515/spp-2013-0002
#'
#' They introduce a Raj Difference estimator that is an extension of a
#' Horvitz-Thompson estimator where, in effect, we divide the sum of
#' the outcomes by a fixed constant rather than the realized sample
#' size to avoid biasing our estimate.
#'
#' Note that the control_formula is ignored--these estimators do not
#' adjust for additioanl controls at this time.
#'
#' @inheritParams linear_model_estimators
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
middleton_aronow_estimator <- function( formula,
                                        data = NULL,
                                        control_formula = NULL,
                                        aggregated = FALSE,
                                        include_block_estimates = FALSE ) {

    #if ( !is.null( control_formula ) ) {
    #    message( "Note: Adjusting for control variables not yet implemented in the HT and Raj estimators" )
    #}

    if ( !aggregated ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula = control_formula )
        # Collapse to clusters
        data_agg = aggregate_data( data, control_formula )
        control_formula = attr( data_agg, "control_formula" )
    } else {
        data_agg = data
    }
    has_block = "blockID" %in% names( data )

    if ( !has_block ) {
        data_agg$blockID = "Uni"
    }
    n_block = length( unique( data_agg$blockID ) )

    M = nrow( data_agg )
    m_t = sum( data_agg$Z )

    data_agg <- data_agg %>%
        mutate( Ytot = Ybar * n )

    # Calculate relationship of cluster total and cluster size
    alpha = 0
    if ( sd( data_agg$n ) > 0 ) {
        Malpha = lm( Ytot ~ n, data=data_agg )
        alpha = coef(Malpha)[[2]]
    }

    # Calculate totals and U
    data_agg <- data_agg %>%
        group_by( blockID ) %>%
        mutate( n_bar = mean( n ) ) %>%
        ungroup() %>%
        mutate( U = Ytot - alpha*(n - n_bar) )

    # Calculate point estimates
    tots <- data_agg %>%
        group_by( blockID, Z ) %>%
        summarise( Ybar = mean(Ytot),
                   Ubar = mean(U),
                   S2Y = var(Ytot),
                   S2U = var(U),
                   J = n(),
                   N = sum(n) )

    tots <- pivot_wider( tots,
                         names_from = Z,
                         values_from = c( Ybar, Ubar, S2Y, S2U, J, N ) ) %>%
        mutate( J = J_0 + J_1,
                N = N_0 + N_1,
                ATE_HT = (J/N) * (Ybar_1 - Ybar_0),
                ATE_Raj = (J/N) * (Ubar_1 - Ubar_0),
                SE_HT = ((J^2)/(N^2)) * (S2Y_1/J_1 + S2Y_0/J_0),
                SE_Raj = ((J^2)/(N^2)) * (S2U_1/J_1 + S2U_0/J_0) )
    stopifnot( nrow( tots ) == n_block )

    ATE_HT = calc_agg_estimate( tots$N, tots$ATE_HT, tots$SE_HT, df=rep( NA, length(tots$N) ) )
    ATE_Raj = calc_agg_estimate( tots$N, tots$ATE_Raj, tots$SE_Raj, df=rep( NA, length(tots$N) ) )

    df = M - 2*n_block

    res <- tibble(
        method = c( "DB_HT", "DB_Raj" ),
        ATE_hat = c( ATE_HT$ATE_hat, ATE_Raj$ATE_hat ),
        SE_hat = c( ATE_HT$SE_hat, ATE_Raj$SE_hat ),
        df = c( df, df ),
        p_value = two_sided_p(ATE_hat, SE_hat, df )
    )

    if ( include_block_estimates ) {
        attr( res, "blocks" ) <- tots
    }

    res
}






#' Estimate ATE via canonical weights on the cluster means
#'
#' These estimators use the cluster means and the table of weights
#' shown in the main paper to estimate the ATE across all combinations
#' of weights at the block and the cluster level.
#'
#' There are no standard errors for this method.  It is purely to
#' validate how the other estimators can be represented as given with
#' the weight formula in the paper.
#'
#' @inheritParams linear_model_estimators
#'
#' @return tibble of estimates using different varieties of the
#'   methods described in the paper.
#'
#' @export
canonical_weight_estimators <- function( formula,
                                         data = NULL,
                                         data_agg = NULL ) {

    stopifnot( !is.null( data ) )
    data = make_canonical_data( formula=formula, data=data )

    if ( is.null( data_agg ) ) {
        # Collapse to clusters
        data_agg = aggregate_data( data )
    }

    has_block = "blockID" %in% names( data )
    if ( !has_block ) {
        data_agg$blockID = "Uni"
    }

    calc_est <- function( data_agg, cw, bw ) {
        if ( length( bw ) == 1 ) {
            bw = rep( bw, nrow(data_agg) )
        }
        if ( length( cw ) == 1 ) {
            cw = rep( cw, nrow(data_agg) )
        }
        if ( length( bw ) != nrow(data_agg) ) {
            stop( "Block weight vector must be same length as data" )
        }
        if ( length( cw ) != nrow(data_agg) ) {
            stop( "Cluster weight vector must be same length as data" )
        }

        data_agg$cw = cw
        data_agg$bw = bw
        s_dat <- data_agg %>%
            group_by( blockID, Z, bw ) %>%
            summarise( Ybar = weighted.mean( Ybar, w = cw ), .groups = "drop" )
        stopifnot( nrow( s_dat ) == length( unique( data_agg$blockID ) ) * 2 )
        s_dat <- s_dat %>%
            group_by( Z ) %>%
            summarise( ATE = weighted.mean( Ybar, w = bw ) )
        if ( nrow( s_dat ) != 2 ) {
            return( NA )
        }
        s_dat$ATE[[2]] - s_dat$ATE[[1]]
    }

    data_agg <- data_agg %>%
        group_by( blockID ) %>%
        mutate( N = sum( n ),
                p = weighted.mean( Z, w=n ),
                pc = mean( Z ),
                J = n() )

    person_person = calc_est( data_agg,
                              cw = data_agg$n,
                              bw = data_agg$N )

    person_FE = calc_est( data_agg,
                          cw = data_agg$n,
                          bw = data_agg$N * data_agg$p * (1-data_agg$p) )

    person_block = calc_est( data_agg,
                             cw = data_agg$n,
                             bw = 1 )

    cluster_cluster = calc_est( data_agg,
                                cw = 1,
                                bw = data_agg$J )
    cluster_FE = calc_est( data_agg,
                           cw = 1,
                           bw = data_agg$J * data_agg$pc * (1-data_agg$pc) )

    cluster_block = calc_est( data_agg,
                              cw = 1,
                              bw = 1 )

    form = make_regression_formula( FE = has_block, cluster_RE = TRUE )
    M0 = lmer( form, data=data )
    tau2 = VarCorr(M0)$clusterID[1,1]
    sigma2 = sigma(M0)^2
    formRE = make_regression_formula( FE = FALSE, cluster_RE = TRUE )
    formRE = update( formRE, . ~ . + (1+Z|blockID) )
    M1 = lmer( formRE, data=data )
    eta2 = VarCorr(M1)$blockID[ 2, 2 ]

    cw_MLM = 1 / (tau2 + sigma2 / data_agg$J)
    bw_RIRC = 1 / (eta2 + tau2 / (data_agg$J * data_agg$pc * (1-data_agg$pc) ))

    MLM_FI_cluster = calc_est( data_agg,
                               cw = cw_MLM,
                               bw = data_agg$J )
    MLM_FE = calc_est( data_agg,
                       cw = cw_MLM,
                       bw = data_agg$N * data_agg$p * (1-data_agg$p) )
    MLM_RIRC = calc_est( data_agg,
                         cw = cw_MLM,
                         bw = bw_RIRC )
    MLM_FI_block = calc_est( data_agg,
                             cw = cw_MLM,
                             bw = 1 )

    rs <- tibble( method = c( "Person_Person", "Person_FE", "Person_Block",
                              "MLM_FI_Cluster", "MLM_FE", "MLM_RIRC", "MLM_Block",
                              "Cluster_Cluster", "Cluster_FE", "Cluster_Block" ),
                  ATE_hat = c( person_person, person_FE, person_block,
                               MLM_FI_cluster, MLM_FE, MLM_RIRC, MLM_FI_block,
                               cluster_cluster, cluster_FE, cluster_block ) )
    rs$method = paste0( "WT_", rs$method )

    return( rs )
}






#### Testing #####

if ( FALSE ) {


    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of schools
        , K = 10                          # number of districts
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of district assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of district variation
        , ICC.3 = 0.2 # district intraclass correlation
        , omega.3 = 0.1           # ratio of district effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.1           # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )


    sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
    rm( model.params.list )

    head( sim.data )
    table( sim.data$D.id)
    aa = design_based_estimators( Yobs ~ T.x | S.id | D.id, data=sim.data )
    aa

    # No block-level.
    sim.data$D.id = NULL
    bb = design_based_estimators( Yobs ~ T.x | S.id, data=sim.data )
    bb


    # Unbalanced data
    data( "fakeCRT" )
    head( fakeCRT )
    bb = design_based_estimators( Yobs ~ T.x | S.id | D.id, data = fakeCRT )


    compare_methods( Yobs ~ T.x | S.id | D.id, data = fakeCRT , include_MLM = FALSE, include_agg = FALSE ) %>%
        arrange( SE_hat )
}











