

require(tidyverse)
require(rlang)
require(stargazer)
require( ggthemes )
require(broom)
require(knitr)
require(kableExtra)
require(patchwork)
require(gridExtra)
require(pander)
require(rpart)
require( rpart.plot )
require(RColorBrewer)
require(cowplot)

select <- dplyr::select


DGP_vars <- c("K", "J.bar", "p.bar", "sd.block", "ICC.block", "variable.p", "n.bar",
              "sd.cluster", "ICC.cluster",
              "size.impact.cor", "num.cov.block","cv.nbar", "R2.cluster",
              "num.cov.cluster", "cv.cluster", "num.cov.indiv",
              "empirical", "cv.block", "nbar.J.cor", "size.intercept.cor")


# Utility plotting functions
source( here::here( "analysis/R/create_custom_boxplot.R" ) )
source( here::here( "analysis/R/create_analysis_tree.R" ) )

# And load the data
source( here::here( "analysis/load_and_prep_data.R" ) )

sruns.clean <- sruns.clean %>%
  filter( group == "fullfactor" ) %>%
  select( -group ) 

SErats <- sruns.clean %>%
  group_by( across( c( ID, block_type, all_of( DGP_vars ) ) ) ) %>%
  summarize( SE_ratio = max( SE.s ) / min( SE.s ) ) %>%
  filter( sd.cluster == "0", sd.block == "0" ) %>%
  ungroup()

mean( is.na( SErats$SE_ratio ) )
nrow( SErats )
mean( SErats$SE_ratio > 1.3 )
quantile( SErats$SE_ratio, 0.95 )
quantile( SErats$SE_ratio, 0.99 )
max( SErats$SE_ratio )

nrow(SErats)
filter( SErats, SE_ratio > 1.3 ) %>%
  select( ID, block_type, J.bar, SE_ratio ) %>%
  group_by( block_type ) %>%
  summarise( n = n(),
             mn_rat = mean(SE_ratio),
             max_rat = max(SE_ratio) ) 

aa <- filter( SErats, SE_ratio > 1.3, block_type=="Same" ) %>%
  dplyr::select( ID, all_of( DGP_vars ), SE_ratio ) 
table( aa$cv.cluster )
#View(aa)
aa$ID
1 - 1 / aa$SE_ratio


library( partykit )
cat <- create_analysis_tree_v2(
  SErats,
  outcome_var = "SE_ratio",
  predictor_vars = c( "block_type", DGP_vars ),
  tree_title = "",
  num_leaves = 5
)



