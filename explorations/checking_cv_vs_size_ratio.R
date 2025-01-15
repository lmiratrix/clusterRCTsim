
# This code calculates the relationship between the coefficient of
# variation (that we get from describe_clusterRCT) and the size.ratio
# we use for our DGP (for how much cluster size varies)


library( tidyverse )
library( clusterRCT )

gen_cv <- function( n.bar, size.ratio ) {
  ns <- blkvar::block_distn( 100,
                                      n.bar = n.bar,
                                      size.ratio = size.ratio,
                                      min.size = 4 )

  sd( ns ) / mean( ns )  
}

sr <- expand_grid( size.ratio = seq( 0, 1, length.out = 200 ),
                   n.bar = c( 6, 10, 20, 40, 100 ) )

sr$cv = pmap_dbl( sr, gen_cv )

ggplot( sr, aes( cv, size.ratio, color = factor( n.bar ) ) ) +
  geom_smooth( se = FALSE ) +
  theme_minimal() +
  labs( y = "Size Ratio", x = "CV" ) +
  coord_fixed() +
  geom_abline( intercept = 0, slope = 1, linetype = "dashed" )


head( sr )
# fit a loess curve to predict cv from size ratio
loess_mod <- loess( size.ratio ~ n.bar + cv, data = sr ) #filter( sr, n.bar > 10 ) ),
                    
sr$size.ratio.loess = predict( loess_mod,
                       newdata = sr )
                       
ggplot( sr, aes( size.ratio.loess, size.ratio, color = factor( n.bar ) ) ) +
  geom_smooth( se = FALSE ) +
  theme_minimal() +
  labs( x = "Predicted size.ratio", y = "size.ratio" ) +
  coord_fixed() +
  geom_abline( intercept = 0, slope = 1, linetype = "dashed" )



# Look at actual data

dt <- readRDS( "empirical/results/master240528_combined.RDS" )
# View( dt )
ggplot( dt, aes(Jcv ) ) +
  geom_histogram()

tb = tibble( cv = c( 0, 0.5, 1.0 ),
             n.bar = 20 )
tb$size.ratio = predict( loess_mod, newdata=tb )
tb

# The above are our levels of size.ratio we will use for our
# simulation.


