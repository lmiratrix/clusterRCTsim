## code to prepare `fakeCRT` dataset goes here


set.seed( 1039 )
library( tidyverse )
library( PUMP )

library( clusterRCT )

model.params.list <- list(
    M = 1                            # number of outcomes
    , J = 30                          # number of schools
    , K = 10                          # number of districts
    , nbar = 10                       # number of individuals per school
    , S.id = NULL                     # N-length vector of school assignments
    , D.id = NULL                     # N-length vector of district assignments
    , Xi0 = 0                         # scalar grand mean outcome under no treatment
    , MDES = 0.125            # minimum detectable effect size
    , R2.3 = 0.1              # percent of district variation
    , ICC.3 = 0.2             # district intraclass correlation
    , omega.3 = 0.1           # ratio of district effect size variability
    , R2.2 = 0.1              # percent of school variation
    , ICC.2 = 0.2             # school intraclass correlation
    , omega.2 = 0.1           # ratio of school effect size variability
    , R2.1 = 0.1    # percent of indiv variation explained
)


sim.data <- PUMP::gen_sim_data( d_m = "d3.2_m3ff2rc", model.params.list, Tbar = 0.5 )
fakeCRT = slice_sample(sim.data, n = nrow(sim.data)*0.5 )
head( fakeCRT )
k = length( unique( fakeCRT$S.id ) )
k
uu = unique( fakeCRT$S.id )
levels(uu)
cids = sample( uu, 200, prob = 50 + as.numeric( uu ) )
fakeCRT = filter( fakeCRT, S.id %in% cids )

library( lme4 )
M <- lmer( Yobs ~ 1 + T.x + (1|S.id) + (1|D.id), data=fakeCRT )
arm::display(M)

fakeCRT$S.id = droplevels( fakeCRT$S.id )

describe_clusterRCT( Yobs ~ T.x | S.id | D.id, data=fakeCRT )

v <- (1.06^2 + 0.56^2+ 0.66^2)
v
0.19 / v

0.66^2 / (1.06^2 + 0.56^2+ 0.66^2)
0.56^2 / (1.06^2 + 0.56^2+ 0.66^2)

head( fakeCRT )
fakeCRT$X = sample( LETTERS[1:4], nrow(fakeCRT), replace=TRUE )
fakeCRT = as_tibble(fakeCRT )

sdY0 = sd( fakeCRT$Yobs[ fakeCRT$T.x==1 ] )
fakeCRT <- mutate( fakeCRT,
                    Yobs = Yobs / sdY0 )
head(fakeCRT)


usethis::use_data(fakeCRT, overwrite = TRUE)


#### Make a fake data with missing outcomes and covariate values ####
table( fakeCRT$D.id )


fakeCRT$T.x[ fakeCRT$D.id == 2 ] = 1
fakeCRT$T.x[ fakeCRT$D.id == 9 ] = 0
nrow(fakeCRT)
mss = sample( nrow(fakeCRT), 100 )

table( fakeCRT$S.id, fakeCRT$D.id )

fakeCRT$Yobs[mss[1:50]] = NA
fakeCRT$D.id[ fakeCRT$D.id == 3 ] = NA
fakeCRT$S.id[ fakeCRT$S.id %in% c(164, 80, 192 ) ] = NA
fakeCRT$X[ mss[40:70]] = NA
fakeCRT$C.ijk[ mss[60:100] ] = NA
fakeCRT$T.x[ fakeCRT$D.id == 5 & fakeCRT$T.x == 0 ] = NA
fakeCRT$Yobs[ fakeCRT$D.id == 6 & fakeCRT$T.x == 1 ] = NA

fakeBrokeCRT = fakeCRT
usethis::use_data(fakeBrokeCRT, overwrite = TRUE)




#### A second fake dataset ####

source( here::here( "simulations/mini_simulation_varied.R" ) )
fakeCRT2 = make_simple_block_example( tx_scale = 0.4)

dd <- fakeCRT2 %>% group_by( T.x ) %>%
    summarise( n = n(),
               sd = sd( Yobs ) )
dd

usethis::use_data(fakeCRT2, overwrite = TRUE)

