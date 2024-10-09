


#' Describe a cluster RCT dataset
#'
#' Calculate various statistics describing total clusters, etc.  If
#' formula is NULL, data assumed to be in canonical form.
#'
#' Will patch data via patch_data() to handle missing values.
#'
#' @inheritParams compare_methods
#' @seealso [patch_data_set()]
#'
#' @export
describe_clusterRCT <- function( formula = NULL,
                                 data = NULL,
                                 control_formula = NULL,
                                 warn_missing = TRUE ) {

    if ( !is.null( formula ) ) {
        # Deal with multiple outcomes
        lhs_vars <- formula.tools::lhs.vars(formula)
        if (length(lhs_vars) > 1) {
            res <- purrr::map(
                lhs_vars,
                function(v) {
                    # see https://stackoverflow.com/questions/37824625/how-to-modify-the-left-side-of-a-formula
                    outcome <- as.name(v)
                    formula1 <- as.formula(bquote(.(outcome) ~ .(formula.tools::rhs(formula))))
                    res <- describe_clusterRCT(
                        formula1,
                        data = data,
                        control_formula = control_formula ) %>%
                        dplyr::mutate(outcome = v, .before=n)
                }) %>%
                purrr::list_rbind()
            return(res)
        } else {
            data = make_canonical_data( formula=formula, data=data,
                                        control_formula = control_formula,
                                        drop_missing = FALSE,
                                        warn_missing = FALSE )
        }
    }

    cnames = NULL
    if ( !is.null( control_formula ) ) {
        cnames = setdiff( colnames(data), c("Yobs", "Z", "clusterID", "blockID" ) )
    }

    # count up missing data and report, then patch dataset.
    miss = colSums(is.na(data))

    data = patch_data_set( NULL, data=data, control_formula = control_formula,
                           warn_missing = warn_missing )

    is_blocked = "blockID" %in% colnames(data)

    if ( is_blocked ) {
        data = patch_singleton_blocks(NULL, data=data,
                                      warn_missing = warn_missing )
    }

    control_formula = attr( data, "control_formula" )

    if ( is_blocked ) {
        K = length( unique( data$blockID ) )
    } else {
        K = 1
        data$blockID = ".single"
    }

    # Get cluster sizes
    sizes = data %>%
        group_by( blockID, clusterID, Z ) %>%
        summarise( n = n(), .groups = "drop" )

    cstat = sizes %>%
        summarise( J = n(),
                   nbar = mean(n),
                   ncv = sd(n) / nbar,
                   n.25 = quantile( n, 0.25 ),
                   n.75 = quantile( n, 0.75 ),
                   n.IQR = n.75 - n.25,
                   p.tx = mean(Z) )

    sstat = sizes %>%
        group_by( blockID ) %>%
        summarise( J = n(),
                   n = sum(n),
                   p.tx = mean(Z) ) %>%
        summarise( K = n(),
                   Jbar = mean(J),
                   Jcv = sd(J) / Jbar,
                   J.25 = quantile( J, 0.25 ),
                   J.75 = quantile( J, 0.75 ),
                   J.IQR = J.75 - J.25,
                   n_block = mean(n),
                   n_block_cv = sd( n ) / n_block,
                   tx.avg = mean( p.tx ),
                   tx.cv = sd( p.tx ) / tx.avg,
                   tx.25 = quantile( p.tx, 0.25 ),
                   tx.75 = quantile( p.tx, 0.75 ),
                   tx.IQR = tx.75 - tx.25 )


    ICCs = calc_ICCs(data, is_blocked = is_blocked )

    stats <- tibble( n = nrow(data) )
    stats = bind_cols( stats, cstat, sstat )
    stats$cluster_ICC = ICCs$C.ICC
    stats$block_ICC = ICCs$S.ICC
    stats$sdY0 = sd( data$Yobs[ data$Z == 0] )

    # Calculate some notes about block distribution.
    notes = c()
    tb = table( sizes$blockID, sizes$Z )

    stats$matched_pairs = FALSE
    stats$num_singletons = sum( tb == 1 )

    if ( any( tb == 1 ) ) {
        if ( all( tb <= 1 ) ) {
            notes = c( notes, "matched pairs design" )
            stats$matched_pairs = TRUE
        } else {
            nblk = apply( tb == 1, 1, max )
            notes = c( notes, glue::glue( "{sum(nblk)} blocks have at least one singleton treated or control clusters" ) )
            nblk = apply( tb == 1, 1, min )
            notes = c( notes, glue::glue( "{sum(nblk)} blocks have exactly 2 clusters, one treated and one control" ) )
        }
    }

    stats$num_doubletons = sum( tb == 2 )
    if ( any( tb == 2 ) ) {
        notes = c( notes, glue::glue( "{sum(tb==2)} within-block treatment arms have exactly 2 clusters." ) )
    }
    attr( stats, "notes" ) <- notes

    class( stats ) <- c( "clusterRCTstats", class( stats ) )

    # R2 values
    if ( !is.null( cnames ) ) {
        stats <- bind_cols( stats, calc_covariate_R2s(data, is_blocked = is_blocked ) )
    }

    if ( sum( miss ) > 0 ) {
        stats$.missing = list( miss )
    }

    stats
}



#' Calculate ICCs
#'
#' Use a multilevel model to calculate ICCs.  Treated units included,
#' using a constant treatment effect model.
#'
calc_ICCs <- function( data, is_blocked = TRUE ) {

    require( lme4 )

    if ( is_blocked ) {
        form = Yobs ~ 1 + Z + (1 | blockID ) + (1 | clusterID )
    } else {
        form = Yobs ~ 1 + Z + (1 | clusterID )
    }
    #if ( !is.null( control_formula ) ) {
    #    form = update( control_formula, Yobs ~ . + 1 + Z + (1 | blockID ) + (1 | clusterID ) )
    #}
    M = lmer( form, data=data )
    #arm::display(M)

    a = VarCorr( M )

    sigma2 = sigma( M )^2

    if ( is_blocked ) {
        C.ICC = as.numeric( a$clusterID )
        S.ICC = as.numeric( a$blockID )
        tvar = (C.ICC + S.ICC + sigma2)
        S.ICC = S.ICC / tvar
        C.ICC = C.ICC / tvar
    } else {
        C.ICC = as.numeric( a$clusterID )
        S.ICC = NA
        tvar = (C.ICC + sigma2)
        C.ISS = C.ICC / tvar
    }

    list( C.ICC = C.ICC, S.ICC = S.ICC )
}


#' Calculate R2 of level 1, level 2, and block.
#'
#' All variables are group mean centered to make level 1 only and
#' level 2 only covariates.
#'
#' For R2.1: Calculate by comparing model with fixed effects for
#' cluster vs fixed effects for cluster AND level 1 variables.
#'
#' For R2.2: Calculate by comparing multilevel model with fixed
#' effects for block and random effects for cluster to model with
#' level 2 covariates as well.
#'
#' block is always the fixed effects for block, so no covariates
#' allowed at that level.  In all models block fixed effects are
#' included, so the level 2 covariates R2 are additional explanatory
#' power beyond their representing block differences (e.g., they are
#' considered block centered).
#'
#' @importFrom purrr map_dbl
#' @param data Data in canonical form.
#' @param pooled Calculate pooled Tx and Co groups, attempting to
#'   remove systematic shift of average treatment effect.  If FALSE,
#'   subset to control observations only.
#'
calc_covariate_R2s <- function( data, pooled = TRUE, is_blocked = TRUE ) {

    # Calculating R2.1
    if ( !pooled ) {
        data = filter( data, Z == 0 )
        data$Z = NULL
    }

    # Covert covariates to numerical, if they are not already.
    if ( is_blocked ) {
        dm = model.matrix( Yobs ~ . - clusterID - blockID, data=data ) %>%
            as.data.frame() %>%
            dplyr::select( -`(Intercept)` )
    } else {
        # NOTE: This annoying thing is neccessary as model.matrix
        # bails before we can subtract out blockID, so prior block
        # doesn't work.
        data$blockID = NULL
        dm = model.matrix( Yobs ~ . - clusterID, data=data ) %>%
            as.data.frame() %>%
            dplyr::select( -`(Intercept)` )

    }
    cnames = colnames(dm)

    dm$clusterID = data$clusterID
    if ( !is_blocked ) {
        dm$blockID = ".single"
    } else {
        dm$blockID = data$blockID
    }

    # Add group means and group mean-centered versions of all
    # covariates
    result <- dm %>%
        group_by(blockID, clusterID) %>%
        mutate(
            across( everything(),
                    list(mn = ~ mean(.),
                         cent = ~. - mean(.)),
                    .names = "{.col}_{.fn}")
        ) %>%
        ungroup() %>%
        dplyr::select( -all_of( cnames ) )

    # Drop centered covariates with no variation (i.e., level 2
    # covariates get dropped from level 1 consideration)
    cents = which( endsWith( colnames(result), "_cent" ) )
    sds = purrr::map_dbl( cents, ~ sd( result[[ . ]] ) )
    ncov.1 = sum( sds != 0 )
    cov.1 = paste0( colnames(result)[ c( cents[sds != 0] ) ], collapse = ", " )

    result = result[ -c( cents[sds == 0] ) ]
    ncov.2 = sum( endsWith( colnames(result), "_mn" ) ) - 1
    cov.2 = paste0( setdiff( colnames(result)[ endsWith( colnames(result), "_mn" ) ], "Z_mn" ),
                    collapse = ", " )

    result$Yobs = data$Yobs

    noLvl2 = result %>%
        dplyr::select( !ends_with( "_mn" ) )
    noLvl1 = result %>%
        dplyr::select( !ends_with("_cent" ) )

    # Calculate R2.1
    if ( is_blocked ) {
        # browser()
        Mfull = lmer( Yobs ~ 1 + . - blockID - clusterID + (1|clusterID), data=noLvl2 )
        Mnull = lmer( Yobs ~ 1 + (1|clusterID), data=noLvl2 )
        R2.1 = (sigma( Mnull )^2 - sigma(Mfull)^2) / sigma(Mnull)^2


        # Calculate R2.2
        M2null = lmer( Yobs ~ 1 + . - clusterID + blockID + (1|clusterID), data=noLvl2 )
        M2full = lmer( Yobs ~ 1 + . - clusterID + blockID + (1|clusterID), data=result )

        vFull = as.numeric( VarCorr( M2full )$clusterID )
        vNull = as.numeric( VarCorr( M2null )$clusterID )
        R2.2 = pmax( 0, (vNull - vFull) / vNull )

    } else {
        # browser()
        noLvl2$blockID = NULL
        Mfull = lmer( Yobs ~ 1 + . - clusterID + (1|clusterID), data=noLvl2 )
        Mnull = lmer( Yobs ~ 1 + (1|clusterID), data=noLvl2 )
        R2.1 = (sigma( Mnull )^2 - sigma(Mfull)^2) / sigma(Mnull)^2


        # Calculate R2.2
        result$blockID = NULL
        M2null = lmer( Yobs ~ 1 + . - clusterID  + (1|clusterID), data=noLvl2 )
        M2full = lmer( Yobs ~ 1 + . - clusterID  + (1|clusterID), data=result )

        vFull = as.numeric( VarCorr( M2full )$clusterID )
        vNull = as.numeric( VarCorr( M2null )$clusterID )
        R2.2 = pmax( 0, (vNull - vFull) / vNull )
    }

    # Pack and ship!
    R2s <- tibble( R2.1 = R2.1,
                   ncov.1 = ncov.1,
                   R2.2 = R2.2,
                   ncov.2 = ncov.2
    )

    # Add Hedges SEs
    calc_hedge_SE <- function( R2, n ) {
        R2 = pmax( 0.01, pmin( R2, 0.99 ) )
        sqrt( 4 * R2 * (1-R2)^2 / n )
    }
    R2s$SER2.1 = calc_hedge_SE( R2.1, nrow(result) )
    R2s$SER2.2 = calc_hedge_SE( R2.2, length(unique(result$clusterID)) )

    # Add in covariate lists
    R2s$cov.1 = cov.1
    R2s$cov.2 = cov.2

    R2s
}



#' Make table of statistics by block
#'
#' Given individual/school/block data of a blocked, cluster RCT,
#' calculate statistics for each block (block).
#'
#' @param check_data_integrity TRUE means runs some checks and give
#'   errors if data fails them (e.g., incorrectly processed treatment
#'   vector.). FALSE means calculate statistics without these checks.
#' @export
make_block_table <- function(  formula = NULL,
                               data = NULL,
                               control_formula = NULL,
                               check_data_integrity = FALSE ) {

    data = make_canonical_data( formula=formula, data=data,
                                control_formula = control_formula )

    if ( check_data_integrity ) {
        check_data_integrity( data )
    }

    n <- nrow( data )

    K = length( unique( data$blockID ) )

    sizes = data %>%
        group_by( blockID, clusterID, Z ) %>%
        summarise( n = n(), .groups = "drop" )

    sstat = sizes %>% group_by( blockID ) %>%
        summarise( J = n(),
                   nbar = mean(n),
                   ncv = sd(n) / nbar,
                   n.25 = quantile( n, 0.25 ),
                   n.75 = quantile( n, 0.75 ),
                   n.IQR = n.75 - n.25,
                   p.tx = mean(Z) )
    sstat
}





#' Pretty print result from describe cluster RCT call
#'
#' @export
#' @param x A clusterRCTstats object.
#' @param ... No extra options passed.
#' @family clusterRCTstats
#'
print.clusterRCTstats <- function( x, ... ) {
    stopifnot( is.clusterRCTstats(x) )

    if ( nrow( x ) > 1 ) {
        rr <- print.data.frame( x, ... )
        return( rr )
    }

    if (  is.null( x$K ) || x$K == 1 ) {
        scat( "Cluster RCT: %d units in %d clusters\n", x$n, x$J )
    } else {
        scat( "Cluster RCT: %d units in %d clusters across %d blocks\n",
              x$n, x$J, x$K )
    }


    if ( !is.null( x$K ) && x$K != 1 ) {
        scat( "Block Statistics:\n\tAvg units/block: %.2f (coef var %.2f)\n\tAvg clusters/block: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f-%.1f w/ IQR = %.1f\n\tICC: %.2f\n\tAvg tx: %.2f (coef var %.2f)\n\t25-75 Quantile tx: %.2f-%.2f w/ IQR %.2f\n",
              x$n_block, x$n_block_cv,
              x$Jbar, x$Jcv,
              x$J.25, x$J.75, x$J.IQR,
              x$block_ICC,
              x$tx.avg, x$tx.cv,
              x$tx.25, x$tx.75, x$tx.IQR)
        if ( hasName( x, "R2.3" ) ) {
            scat( "\tR2.3: %.2f (%d fixed effects)\n", x$R2.3, x$K)
        }
    }

    scat( "Cluster Statistics:\n\tAvg size: %.2f (coef var %.2f)\n\t25-75 Quantiles: %.1f-%.1f w/ IQR = %.1f\n\tICC: %.2f\n\tProp treated: %.2f\n",
          x$nbar, x$ncv,
          x$n.25, x$n.75, x$n.IQR,
          x$cluster_ICC,
          x$p.tx )
    if ( hasName( x, "R2.2" ) ) {
        scat( "\tR2.2: %.2f (%d covariates)\n", x$R2.2, x$ncov.2)
    }


    scat( "Unit Statistics:\n\tstddev( Y0 ): %.2f\n",
          x$sdY0 )
    if ( hasName( x, "R2.1" ) ) {
        scat( "\tR2.1: %.2f (%d covariates)\n", x$R2.1, x$ncov.1)
    }

    if ( ".missing" %in% names(x) ) {
        scat( "missing data counts:\n" )
        print( x$.missing[[1]] )
    }

    notes = attr( x, "notes" )
    if ( !is.null(notes) && ( length( notes ) > 0 ) ) {
        cat( "Notes:\n" )
        for ( n in notes ) {
            cat( "\t", n, "\n" )
        }
    }

    invisible( x )

}



#' Is object a clusterRCTstats object?
#'
#' @export
#' @aliases clusterRCTstats
#' @param x the object to check.
#' @family clusterRCTstats
is.clusterRCTstats = function( x ) {
    inherits(x, "clusterRCTstats")
}




#' Cast cluster RCT info result to data.frame
#'
#' @export
#' @aliases clusterRCTstats
#' @param x the clusterRCTstats object to covert
#' @family clusterRCTstats
as.data.frame.clusterRCTstats = function( x ) {
    x$.missing = NULL
    class(x) = "list"
    x <- as.data.frame( x )
    rownames(x) = 1:nrow(x)

    x
}





count_block_sizes <- function( clusterID, blockID ) {
    tt = tibble( clusterID = clusterID, blockID = blockID )
    t2 <- tt %>%
        group_by( blockID ) %>%
        summarize( n = n(),
                   J = length( unique( clusterID ) ) )

    t2
}




#### Testing code ####

if ( FALSE ) {

    library( tidyverse )
    library( clusterRCT )

    set.seed( 1039 )

    # Make some data
    model.params.list <- list(
        M = 1                            # number of outcomes
        , J = 30                          # number of clusters
        , K = 10                          # number of blocks
        , nbar = 10                       # number of individuals per school
        , S.id = NULL                     # N-length vector of school assignments
        , D.id = NULL                     # N-length vector of block assignments
        , Xi0 = 0                         # scalar grand mean outcome under no treatment
        , MDES = 0.125            # minimum detectable effect size
        , R2.3 = 0.1              # percent of block variation
        , ICC.3 = 0.2             # block intraclass correlation
        , omega.3 = 0.2           # ratio of block effect size variability
        , R2.2 = 0.1              # percent of school variation
        , ICC.2 = 0.2             # school intraclass correlation
        , omega.2 = 0.0          # ratio of school effect size variability
        , R2.1 = 0.1    # percent of indiv variation explained
    )

    data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )

    data = slice_sample( data, n = nrow(data) / 2 )
    dd = sample( unique(data$S.id), 10 )
    data = filter( data, !( S.id %in% dd ) )
    data <- clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=data )
    head( data )


    # Describe the data
    dd <- describe_clusterRCT( data )
    class( dd )
    dd2 = dd

    b = bind_rows( dd, dd2 )
    class(b)
    as.data.frame(b)

    data( "fakeCRT" )
    head( fakeCRT )
    dsc = describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    class( dsc )
    print( dsc )

    # Check R2 calculations
    head( fakeCRT )
    dsc = describe_clusterRCT( formula = Yobs ~ T.x | S.id | D.id, data=fakeCRT,
                               control_formula = ~ V.k + X.jk + C.ijk )
    as.data.frame(dsc)
    print( dsc )



    # Looking at R2 calculation

    dd = data.frame( cid = rep( 1:10, 1:10 ) )
    dd = mutate( dd,
                 Ybar = cid,
                 Z = as.numeric( cid %% 3 == 1 ) )
    dd2 = bind_rows( dd, dd )
    dd$cid = dd$cid + 100
    dd3 = bind_rows( dd, dd2 )
    dd$cid = dd$cid - 100
    dd2$Ybar = dd2$Ybar + 1
    dd3$Ybar = dd3$Ybar + 2
    dd = bind_rows( A = dd, B = dd2, C = dd3, .id = "sid" )
    dd$Y = dd$Ybar + rnorm( nrow(dd) )
    dd <- mutate( dd,
                  Z = ifelse( cid <= as.numeric(as.factor(sid))*2, 1, 0 ) )
    head(dd)
    skimr::skim(dd)

    make_block_table(  Y ~ Z | cid | sid, data=dd )
    describe_clusterRCT( Y ~ Z | cid | sid, data=dd)

    dd <- dd %>%
        group_by( sid ) %>%
        mutate(   X1 = Y + rnorm(n(), sd=3),
                  X2 = cid,
                  X3 = mean(Y),
                  X4 = sample( c( "A", "B", "C"), n(), replace=TRUE ),
                  Y = Y + 1 * Z ) %>%
        ungroup()


    data = clusterRCT:::make_canonical_data( Y ~ Z | cid | sid, data=dd,
                                             control_formula = ~ X1 + X2 + X3 + X4 )
    head( data )


    calc_covariate_R2s( data )

}
