

# This file looks at the simulation scenarios and evaluates how
# closely the generated grids correspond to the target parameters.

library( tidyverse )

scenarios = read_rds( here::here( "results/scenarios.rds") )
scenarios

# Each scenario has a stored grid of blocks
knitr::kable( scenarios[3,] %>% dplyr::select(-grid) %>% pivot_longer( cols=everything() ) )
scenarios$grid[[3]]

# Comparing the actual number of clusters to the target number of clusters


nrow(scenarios )
ugh <- scenarios %>%
  dplyr::select( ID, J.bar, grid ) %>%
  unnest( cols = grid ) %>%
  group_by( ID, J.bar ) %>%
  summarise( J.bar.true = mean( J ),
             n.bar.true = mean( n.bar ), .groups = "drop" )


mean( ugh$J.bar.true == ugh$J.bar )

summary( ugh$J.bar.true - ugh$J.bar )

ggplot( ugh,
        aes( J.bar, J.bar.true, group=J.bar ) ) +
  geom_boxplot() +
  geom_abline( slope=1, intercept=0, color="red", linetype="dashed" )


dels <- abs( ugh$J.bar.true - ugh$J.bar ) / ugh$J.bar
summary( dels )

mean( dels > 0.5 )


dels <- abs( ugh$J.bar.true - ugh$J.bar )
mean( dels > 5 )
