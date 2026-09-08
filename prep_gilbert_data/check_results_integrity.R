
# Quick sanity checks on the final combined results (renamed from the
# original "integrety_checks.R").

library( tidyverse )

comb = read_csv( here::here( "empirical/results/jg_results.csv" ) )
head( comb )

table( comb$Project, comb$blocked )


dd <- comb %>%
    dplyr::select( Project, blocked ) %>%
    unique()
table( dd$blocked )
