

# Looking at bias of person weighted....

source( here::here( "R/data_generation_code.R" ) )

grid <- build_districts( K = 10, J.bar = 20, n.bar = 20, p.bar = 0.5,
                         ATE = 0.2, sd.ATE = 0.3,
                         ICC = 0.2,
                         variable.J = TRUE, variable.p = FALSE, variable.nbar = TRUE,
                         num.cov.district = 3, size.intercept.cor = -0.8,
                         size.impact.cor = 0, nbar.J.cor = -0.9, 
                         empirical = TRUE ) %>%
  arrange( J )

grid
grid$J[1] = 4
grid$n.bar[1] = 100
grid$district.tau[1] = -2
#ggplot( grid, aes( n.bar, J ) ) +
#  geom_point()
cor( grid$J, grid$n.bar )  
describe_district_grid( grid )

ac = single.CRT.fixed.district( grid = grid, variable.n = TRUE, size.impact.cor.school = 0.75, 
                                sd.school = 0.5, n.runs = 50, 
                                just.describe.data = TRUE )
ac

a = single.CRT.fixed.district( grid = grid, variable.n = FALSE, size.impact.cor.school = 0, 
                               sd.school = 0.5, n.runs = 100, 
                               .progress = TRUE )
table( a$subrun )
a  

a %>% 
  filter( disfavored != 1 ) %>%
  group_by( method, weight ) %>%
  summarise( SE.b = sd( ATE.hat ) / sqrt( n() ),
             bias.stu = mean( ATE.hat - ATE.finite.student ),
             bias.sch = mean( ATE.hat - ATE.finite.school ) ) %>%
  arrange( weight ) %>%
  knitr::kable( digits = 2 )


