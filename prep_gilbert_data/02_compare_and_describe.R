

# Read all the datasets in the "raw" folder and run the
# compare_methods and describe_clusterRCT functions on them.  Save the
# results to a csv file along with merged in metadata.




# Load the data files into global variables ----

source( here::here( "prep_gilbert_data/load_all_datasets.R" ) )

dir.create( here::here( "prep_gilbert_data/data/clean" ), showWarnings = FALSE, recursive = TRUE )


# Compare all the datasets below a given size threshold ----
if ( FALSE ) {

    # threshold for skipping datasets to run quickly (for debugging, etc.)
    # Set to NULL for no threshold
    size_threshold = NULL # 5000


    #size_threshold = 1000


    # Run compare on all the datasets ----


    # map across the data sets
    tic()
    compare_tab <- map_dfr( ids, compare_all_cov,
                            size_threshold = size_threshold,
                            dataset_list = dfs,
                            .progress = TRUE)
    cli::cli_alert("Finished compare of {length(ids)} datasets" )
    toc()


    write_csv(compare_tab, here::here( "prep_gilbert_data/data/clean/jg_compare.csv") )

}


# Describe all the datasets that were analyzed ----

if ( !file.exists( here::here( "prep_gilbert_data/data/clean/jg_compare.csv") ) ) {
    cli::cli_abort( c(
        "Can't find {.file jg_compare.csv}.",
        "i" = "This file is produced by the {.code if (FALSE) {{ ... }} } block above (lines ~18-43), which is disabled by default because it is slow (it runs {.code compare_methods()} on every dataset).",
        "i" = "On a fresh checkout, open this script, temporarily flip that block to {.code if (TRUE)}, and run it once to generate {.file jg_compare.csv}. You can set {.code size_threshold} (inside that block) to a number like 5000 to skip large/slow datasets for a quicker test run.",
        "i" = "Once {.file jg_compare.csv} exists, flip the block back to {.code if (FALSE)} so you don't accidentally re-run the expensive step on every source()."
    ) )
}

compare_tab <- read_csv( here::here( "prep_gilbert_data/data/clean/jg_compare.csv") )
dim( compare_tab )
head( compare_tab )

rem_id = unique( compare_tab$id )
cli::cli_alert( "Starting to describe remaining {length(rem_id)} datasets" )

# Describe all datasets we have previously and successfully compared
describe_tab <- map_dfr( unique( compare_tab$id ),
                         describe_all_cov,
                         dataset_list = dfs,
                         .progress = TRUE)
describe_tab <- describe_tab  |>
    select( -c(Trt_arm:Outcome_desc, .missing),
            -Pretest)

write_csv( describe_tab, here::here( "prep_gilbert_data/data/clean/jg_describe.csv") )



# combine the compare and describe results ----

vars_c = colnames( compare_tab )
vars_d = colnames( describe_tab )
vars_m = colnames( metadata )
over_d = intersect( vars_c, vars_d )
# NOTE: original script referenced an undefined variable "over" here
# (almost certainly meant "over_d" itself) - this would have thrown an
# error ("object 'over' not found") if run. Fixed below.
over_d <- setdiff( over_d,
                   c( "Project", "Covariates") )
over_m = intersect( vars_c, vars_m )
# compare_tab's "Outcome" is just the internal standardized outcome
# variable name (always "y_std"), not the human-readable outcome label.
# The case-sensitive intersect() above doesn't catch it against
# metadata's "outcome" column, so it survives into combined and later
# collides with the rename() below (which is trying to bring in
# metadata's real, human-readable outcome name as "Outcome"). Drop it
# here so metadata's "outcome" wins.
over_m <- union( over_m, "Outcome" )
over_m <- setdiff( over_m,
                   c( "dataset_id") )

combined <- compare_tab |>
    dplyr::select( -all_of(c( over_d, over_m) ) ) |>
    left_join(describe_tab,
              by = c( "Project", "Covariates" ) ) %>%
    relocate( Project, method, Covariates ) |>
    mutate(dataset_id = parse_number(Project)) |>
    left_join(metadata, by = "dataset_id")

head( combined )

skimr::skim( compare_tab )


# clean up variable names
combined <- combined |>
    rename( Outcome = outcome )
names( combined )

write_csv(combined, here::here( "prep_gilbert_data/data/clean/jg_results.csv") )


# Sanity check on the results ----


# Quick sanity checks on the final combined results

comb = read_csv( here::here( "empirical/results/jg_results.csv" ) )
head( comb )

table( comb$Project, comb$blocked )

dd <- comb %>%
    dplyr::select( Project, blocked ) %>%
    unique()
table( dd$blocked )




if ( FALSE ) {

    # copy jg_results.csv and jg_describe.csv to the final location used
    # by the rest of clusterRCTsim (empirical/results/)
    file.copy( here::here( "prep_gilbert_data/data/clean/jg_results.csv"),
               here::here( "empirical/results/jg_results.csv" ),
               overwrite = TRUE )

    file.copy( here::here( "prep_gilbert_data/data/clean/jg_describe.csv"),
               here::here( "empirical/results/jg_describe.csv" ),
               overwrite = TRUE )

}
