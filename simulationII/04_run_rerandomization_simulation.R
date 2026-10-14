


#
# Rerandomization simulation on all datasets
#
# Output lands directly in simulationII/results_rerand_sim_public/,
# named jg_crt_{id}.rds - this is where
# simulationII/load_master_rerandomization_results.R expects to find it
# (matched against the UniqueID column of
# empirical/results/study_descriptions.csv).


# load libraries
library(tidyverse)
library(glue)
library(clusterRCT)
library(readxl)
library(tictoc)

# clear memory
rm(list = ls())


source( here::here( "prep_gilbert_data/load_all_datasets.R" ) )

targets = metadata %>%
    filter( cluster_rct == 1, blocked == 1 )

dfs <- dfs[ names(dfs) %in% targets$dataset_id ]

ids = names(dfs)

OUTPUT_DIR = here::here( "simulationII/results_rerand_sim_public" )
dir.create( OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE )

source( here::here( "prep_gilbert_data/processing_functions.R" ) )

SIZE_THRESHOLD = 30000


do_one_simulation <- function( id, R = 2, ... ) {

    clean_dat = grab_data( i = id, dataset_list = dfs, size_threshold = SIZE_THRESHOLD )

    if ( is.null( clean_dat ) ) {
        return( NULL )
    }

    res <- run_rerandomize_simulation( y_std ~ treat | cluster_id | block_id,
                                       data = clean_dat$data,
                                       R = R,
                                       summarize_results = FALSE,
                                       parallel = FALSE, ... )

    nm = glue::glue( "jg_crt_{id}" )
    write_rds( res,
               file.path( OUTPUT_DIR, glue::glue( "{nm}.rds" ) ) )
    res
}



# Run the simulation on a single dataset ----

if ( FALSE ) {
    ids
    ids[[5]]
    do_one_simulation( ids[[5]], control_formula = ~ std_baseline )


    do_one_simulation( "7", R = 20 )

    do_one_simulation( "2", R = 2 )

    sim <- do_one_simulation( "16", R = 20 )
    skimr::skim( sim )
    table( sim$runID )

    clean_dat = grab_data( i = "16", dataset_list = dfs )
    filter( clean_dat$data, is.na( cluster_id ) )
    filter( clean_dat$data, is.na( block_id ) )

    id = ids[[1]]
    clean_dat = grab_data( i = id, dataset_list = dfs )
    clean_dat

    R = 2
    clean_dat$data$y2 = clean_dat$data$y + rnorm( nrow( clean_dat$data ) )

    res <- run_rerandomize_simulation( y2 + y_std ~ treat | cluster_id | block_id,
                                       data = clean_dat$data,
                                       R = R,
                                       summarize_results = FALSE,
                                       include_MLM = FALSE )

    nm = glue::glue( "jg_crt_{id}" )
    write_rds( res,
               file.path( OUTPUT_DIR, glue::glue( "{nm}.rds" ) ) )
    res
}



# Run the simulation on all datasets ----

if ( TRUE ) {
    library( furrr )
    future::nbrOfWorkers()
    if ( future::nbrOfWorkers() == 1 ) {
        plan(multisession, workers = parallel::detectCores() - 1 )
    }
    rps = furrr::future_map( ids, do_one_simulation,
                             R = 20,
                             .options = furrr_options(seed = TRUE),
                             .progress = TRUE )
    plan( future::plan( "sequential" ) )
    gc()
} else {
    rps = map( ids, do_one_simulation,
               R = 20,
               .progress = TRUE )
}
