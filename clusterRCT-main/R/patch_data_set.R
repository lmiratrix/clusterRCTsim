

# Patch a dataset with mean imputation of the covariates


#' Fill in missing values and convert data to canonical form
#'
#' This method will fill in means for all missing covariate values and
#' return a dataset where rows with missing outcome, cluster or block
#' membership, and treatment assignment are all dropped and everything
#' else is filled in.  To do this, it will convert categorical
#' covariates to a set of dummy variables (dropping the reference
#' group).
#'
#' It will drop missing data indicators that are co-linear with prior
#' missing data indicators.
#'
#' If there are categorical covariates, then it will convert those to
#' dummy variables regardless of missingness.  The `control_formula`
#' attribute handed back will give the updated control_formula in this
#' case.
#'
#' This attribute will also have all missing value indicators added in
#' as covariates.
#'
#' @export
patch_data_set <- function( formula = NULL,
                            data = NULL,
                            control_formula = NULL,
                            warn_missing = FALSE ) {

    if ( !is.null(formula) ) {
        data = make_canonical_data( formula=formula, data=data,
                                    control_formula = control_formula,
                                    drop_missing = FALSE,
                                    warn_missing = FALSE )
    }

    mod_data = NULL
    if ( !is.null( control_formula ) ) {
        mod_data <- clusterRCT:::expand_control_variables( data, control_formula )
        data = mod_data$data
        control_formula = mod_data$control_formula
    }

    any_missing = FALSE
    n = nrow(data)

    data = filter( data,
                   !is.na( Z ), !is.na( Yobs ),
                   !is.na( clusterID ) )
    if ( is.factor(data$clusterID) ) {
        data$clusterID = droplevels(data$clusterID)
    }

    if ( "blockID" %in% colnames(data) ) {
        data <- filter( data, !is.na( blockID ) )
        if ( is.factor(data$blockID) ) {
            data$blockID = droplevels(data$blockID)
        }
    }

    if ( nrow(data) < n ) {
        any_missing = "dropped"
    }

    if ( is.null( control_formula ) ) {
        if ( warn_missing && (any_missing != FALSE )) {
            warning( glue::glue( "{n - nrow(data)} rows with missing data dropped (of {n} total rows) in patch_data_set" ),
                     call. = FALSE )
        }
        return( data )
    }

    # Locate missing data and do imputation on the columns with
    # missingness.  Also generate a dummy var.  Make sure they are not
    # colinear.
    miss = colSums(is.na(data))

    which_miss = which( miss > 0 )
    n_imputed = 0
    if ( length( which_miss ) > 0 ) {
        any_missing = ifelse( any_missing == FALSE, "imputed", "dropped and imputed" )

        missInd = 0 + is.na( data[ , which_miss, drop=FALSE ] )
        colnames(missInd) = paste0( colnames(missInd) , "_mi" )
        qr = qr( missInd )

        n_imputed = sum( missInd )

        missInd <- missInd[ , qr$pivot[ 1:qr$rank ], drop=FALSE ] %>%
            as.data.frame()

        data <- data %>%
            mutate( across(all_of(which_miss), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

        data <- bind_cols(data, missInd)

        mod_data$control_formula <- update( mod_data$control_formula,
                                            as.formula( paste0( " ~ . + ", paste( colnames(missInd), collapse = " + " ) ) ) )
    }

    if ( !is.null( mod_data ) ) {
        attr( data, "control_formula" ) <- mod_data$control_formula
    } else {
        attr( data, "control_formula" ) <- NULL
    }

    if ( warn_missing && (any_missing != FALSE) ) {
        if ( n_imputed > 0 ) {
            warning( glue::glue( "Missing data {any_missing} in patch_data_set. {n-nrow(data)} rows of {n} dropped.  {n_imputed} values imputed." ),
                 call. = FALSE )
        } else {
            warning( glue::glue( "Missing data {any_missing} in patch_data_set. {n-nrow(data)} rows of {n} dropped." ),
                     call. = FALSE )
        }
    }


    data
}





#' Drop or combine all-tx or all-co blocks
#'
#' Given dataset with some blocks that are all treated or all control,
#' drop those blocks or replace the block ID with canonical new block
#' ID shared by all such blocks.
#'
#' Missing values in the block ID are all considered singletons.
#'
#' Also, depending on pool_clusters, pool the clusters in each of
#' these identified blocks into single clusters.
#'
#' @param drop_data Drop the troublesome blocks if TRUE, pool them if
#'   FALSE.
#' @param pool_clusters If pooling blocks rather than dropping them,
#'   also pool clusters in 100% tx or 100% co blocks into single
#'   cluster.
#' @param warn_missing Say something if anything happens.
#' @return data with the blockID column modified or troublesome rows
#'   of data missing.
#' @export
#'
patch_singleton_blocks <- function( formula = NULL, data,
                                    drop_data = TRUE,
                                    pool_clusters = TRUE,
                                    warn_missing = TRUE ) {


    if ( is.null( formula ) ) {
        formula = Yobs ~ Z | clusterID | blockID
    }

    missing = identify_singleton_blocks( formula, data )
    if ( is.null( missing ) ) {
        return( data )
    }

    parts = deconstruct_var_formula(formula, data)

    blkID = data[[ as.character(parts$blockID) ]]
    drp = is.na( blkID ) | ( blkID %in% missing$blockID )

    if ( drop_data ) {
        if ( warn_missing ) {
            warning( glue::glue( "Dropping all-tx and/or all-co blocks. {sum(drp)} of {length(drp)} rows of data dropped." ),
                 call. = FALSE )
        }
        data = data[ !drp, ]
        if ( is.factor( data[[ parts$blockID ]] ) ) {
            data[parts$blockID] = droplevels( data[[parts$blockID]] )
        }
        return( data )
    }

    if ( warn_missing ) {
        warning( glue::glue( "Pooling all-tx and/or all-co blocks into singleton block. {sum(drp)} rows of data affected." ),
             call. = FALSE )
    }

    new_blk_id = data[[parts$blockID]]
    is_fac = is.factor(new_blk_id)
    new_blk_id = as.character(new_blk_id)

    if ( pool_clusters ) {
        cid = data[[ parts$clusterID ]]
        is_fac_cid = is.factor(cid)
        cid = as.character(cid)
        for ( n in missing$blockID ) {
            cid[ new_blk_id == n ] = paste0( ".", n )
        }
        if ( is_fac_cid ) {
            cid = as.factor(cid)
        }
        data[[ parts$clusterID ]] = cid
    }

    if ( ".pooled" %in% unique( new_blk_id ) ) {
        warning( "Block ID already has a '.pooled' block.  This is likely to cause problems." )
    }

    new_blk_id[ drp ] = ".pooled"

    if ( is_fac ) {
        new_blk_id = as.factor(new_blk_id)
    }
    data[[parts$blockID]] = new_blk_id

    # Now check if still have singleton blocks
    missing = identify_singleton_blocks( formula, data )
    if ( !is.null( missing ) ) {
        warning( "Still have singleton blocks after pooling." )
    }

    return( data )
}






#### For debugging and demonstration ####

if ( FALSE ) {

    data <- tribble( ~ X, ~ R, ~ Q, ~ W,
                     1,  2, 3, 4,
                     NA, 2, 4, 5,
                     NA, NA,4, 6,
                     2,  3, NA, 7)

    data$T.x = c( 1, 1, 0, 0 )
    data$S.id = c( 10, 5, 10, 5 )
    data = bind_rows( data, data, data, data, data)
    data$D.id = rep( c("A","B"), c(8, 12) )
    data$Y = rnorm( nrow(data), sd=10 )
    data$X[1:5] = 1:5
    data$W[5:10] = 4
    data$Y[12] = NA
    data$T.x[9] = NA
    data$D.id[11] = NA
    data$S.id[2] = NA
    data$RR = sample( c("A","B","C"), nrow(data), replace=TRUE )
    data$RR[17] = NA
    data

    ptch <- patch_data_set( Y ~ T.x | S.id | D.id, data=data,
                            control_formula = ~ X + R + Q + W + RR )
    ptch




    #### Code for debugging ####
    data <- clusterRCT:::make_canonical_data( Y ~ T.x | S.id | D.id, data=data,
                                              control_formula = ~ X + R + Q + W + RR,
                                              drop_missing = FALSE )
    data

    formula = Y ~ T.x | S.id | D.id
    control_formula = ~ X + R + Q + W + RR




    library(tidyverse)

    # Replace missing values with the mean of each column
    data <- data %>%
        mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
    data

    # Now, 'data' contains mean-imputed values for missing data


    #data[is.na(data)] <- sapply(data, function(x) mean(x, na.rm = TRUE))
    data

    data( "fakeCRT" )
    data <- clusterRCT:::make_canonical_data( Yobs ~ T.x | S.id | D.id, data=fakeCRT )
    head( data )



    data( "fakeBrokeCRT" )
    fakeBrokeCRT
    patch_singleton_blocks( Yobs ~ T.x | S.id | D.id, data= fakeBrokeCRT )

}



