

require(tidyverse)
require(rlang)
require( ggthemes )
require(knitr)

library(rpart)
library( rpart.plot )

my_t = theme_tufte() + 
  theme( legend.position="bottom", 
         legend.direction="horizontal",
         legend.key.width=unit(1,"cm"),
         panel.border = element_blank(), 
         axis.title=element_text(size=12))
theme_set( my_t )


# Load simulation results
source( here::here( "analysis/Simulation_Write_Up/load_and_prep_data.R" ) )

# Restruct to the core methods
sruns.sub <- plot_data %>% 
  filter(method %in% core.methods) 

# Focus on MLM only (for example)
#sruns.sub <- sruns.sub %>%
#  filter( method %in% MLM.methods )


# Make relative SE
sruns.sub <- sruns.sub  %>% 
  mutate( SE.s.norm=(SE.s/SE.s.median) )

sruns.sub <- sruns.sub %>%
  mutate(
    is_FI = str_detect(method, "FI"),
    is_MLM = str_detect(method, "MLM"),
    is_cw = str_detect(method, "cw"),
    is_pw = str_detect(method, "pw")
  )

# Make the tree ----

# Rename to make variables readable on printed trees
sruns.sub <- sruns.sub %>%
  rename( blk_type = block_type,
          cor_J_tau = size.impact.cor.school,
          num.cov = num.cov.school,
          cv.sch = cv.school,
          sd.sch = sd.school,
          ICC.sch = ICC.school )


# QUESTION: Include method or not?
form = as.formula( "abs(bias.s.estimand) ~ is_FI + is_MLM + is_cw + is_pw + blk_type + J.bar + cor_J_tau + 
                                                 sd.sch + ICC.sch + cv.sch + num.cov+p.bar" )
form


# Initial fit with a low cp
tree <- rpart(form, data = sruns.sub, cp = 0.001)


cp_table <- tree$cptable %>%
  as_tibble()
cp_table

min_xerror <- which.min(cp_table$xerror)
xerror_se = cp_table$xstd[[ min_xerror ]]
xerror_target = cp_table$xerror[[min_xerror]] + xerror_se
cp_best = min( cp_table$CP[  cp_table$xerror <= xerror_target ] )

num_leaves = 10
row <- pmin( nrow(cp_table),
             1 + sum( !( cp_table[ , "nsplit" ] >= ( num_leaves - 1 ) ) ) )
#cp_table[(row-2):(row+2),]
cp_table[row,]
cp = cp_table[ row, "CP" ]
cp
cp_best

tree_p = prune( tree, cp = max( cp, cp_best ) )


# Make an index for coloring
node_means <-  tree_p$frame$yval
norm_means <- (node_means - min(node_means)) / 
  (max(node_means) - min(node_means)) * 100
norm_means = round( norm_means ) + 1


# Create color mapping
library(RColorBrewer)  # Make sure this is loaded
palette_name <- "Greens" 
cols <- colorRampPalette(brewer.pal(9, palette_name))(131)
node_colors <- cols[round(norm_means) + 1]

# Plot with colored nodes
if ( FALSE ) {
  prp(tree_p, 
      node.fun = function(x, labs, digits, varlen) {
        paste0(labs, "\n", round(x$frame$yval[x$nodenum], 2))
      },
      box.col = node_colors,
      border.col = "black",
      fallen.leaves = TRUE,
      cex = 1.2)
}

leaf_fun <- function(x, labs, digits, varlen) {
  ifelse(x$frame$var == "<leaf>",
         sprintf("%.2f", x$frame$yval),
         labs)
}

prp( tree_p,
     box.col = node_colors,
     border.col = "black",
     node.fun = leaf_fun,
     #fallen.leaves = TRUE,
     varlen = -10)
