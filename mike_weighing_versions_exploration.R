

# Mike's test of two different weighting schemes

# This script generates two districts, one with big clusters and one
# with small clusters.
#
# We then compare two different weighting scheemes, one within block
# and one overall.  They give the same point estimates as long as
# there is no covariate adjustment.




dat <- generate_districts(
  # district-level params
  K = 2, J.bar = 20, p.bar=0.5,
  ATE.district = 0.3, sd.district = 0.3,
  ICC.district = 0.3,
  variable.J = TRUE,
  variable.p = FALSE,
  R2.district = 0.2,
  
  # correlation params at district level
  size.ratio.district = 0.6, 
  size.impact.cor.district = 0,
  proptx.impact.cor.district = 0,
  size.intercept.cor.district = 0,
  
  # school-level params
  n.bar = 50, sd.school = 1, ICC.school = 0.20,
  variable.n = TRUE,
  R2.school = 0.2,
  
  # correlation params at school level
  size.ratio.school = 0.6,
  size.impact.cor.school = 0,
  cor.strength.school = 0.75,
  
  # Individual-level params
  R2.indiv = 0.2 )


head( dat )

table( table( dat$did ) )
table( dat$did, dat$sid )

d2 = filter( dat, did == 1 )
dat = bind_rows( dat, d2, d2, d2 )

table( dat$did )


require( estimatr )

weight = "Cluster"


head( dat )
control_formula = ~ X + W
control_formula = NULL
data = clusterRCT:::make_canonical_data( formula= Yobs ~ Z | sid | did, data=dat, 
                                         control_formula=control_formula )

# Make our weights variable
if (weight == "Person") {
  data$.weight <- 1
} else {
  data <- data %>%
    group_by( clusterID ) %>%
    mutate( .weight = 1 / n() ) %>%
    ungroup()
}

head( data )
table( data$.weight )

form = clusterRCT:::make_regression_formula( FE = TRUE, interacted = TRUE,
                                control_formula = control_formula )
form

M0.int <- estimatr::lm_robust( form, data=data, clusters=clusterID, weights = .weight )

tidy( M0.int )


ests <- generate_all_interacted_estimates( M0.int, data,
                                           method = "test",
                                           weight = weight,
                                           se_method = "crve",
                                           use_full_vcov = TRUE )

ests


# Now the "mike weights"

data
data <- data %>% 
  group_by( blockID ) %>%
  mutate( .weight_mike = .weight * n() / length(unique(clusterID) ) )
head( data )

M0.mike <- estimatr::lm_robust( form, data=data, clusters=clusterID, weights = .weight_mike )

tidy( M0.mike )

tidy( M0.int )




