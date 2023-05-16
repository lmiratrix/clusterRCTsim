# ---
# title: "Examine coverage of the different estimators"
# author: "Miratrix, Weiss, Henderson"
# date: "`r Sys.Date()`"
# ---


# Using the results from clean_simulation_results, look at coverage

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library( tidyverse )
library( ggthemes )

my_theme = theme_tufte() + theme( legend.position="bottom", 
                                  legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                  legend.box = "vertical", 
                                  panel.border=element_rect(size=0.3, fill=NA))
theme_set( my_theme )


if ( FALSE )  {
  # Make temp file of coverage results.
  cat( "Loading the big data file\n" )
  load( "../results/processed_simulation_results.RData" )
  
  head( sruns )
  
  
  
  cov = sruns %>% select( method, ID, cover.fp, cover.fs, cover.sp, cover.ss )
  cov = gather( cov, cover.fp, cover.fs, cover.sp, cover.ss, key="estimand", value="coverage" )
  head( cov )
  
  cov = merge( cov, scenarios, by="ID" )
  nrow( cov )
  
  head( cov )
  
  write_rds(cov, path="../results/coverage_results.rds" )
}

# Load coverage results (from aggregated simulation)
cov = read_rds("../results/coverage_results.rds" )



ggplot( cov, aes( method, coverage ) ) +
  facet_wrap( ~ estimand ) +
  geom_boxplot() +
  geom_hline(yintercept=0.95, col="red" ) 


head( cov )
cov = merge( cov, blkvar::method_characteristics(), by="method" )

head( cov )
cov$coverage = 100 * cov$coverage

# Keep correct coverage rates aligned with target estimand.
covsub = filter( cov, (estimand=="cover.fp" & population == "finite" & weight=="person") |
                  (estimand=="cover.fs" & population == "finite" & weight=="site") |
                (estimand=="cover.sp" & population == "superpop" & weight=="person") |
                (estimand=="cover.ss" & population == "superpop" & weight=="site") )
nrow( covsub )
table( covsub$method )


# Coverage of estimators for target estimand.
ggplot( covsub, aes( method, coverage, col=population ) ) +
  facet_grid(  weight + population ~ J, scales="free_y" ) +
  geom_boxplot() +
  geom_hline(yintercept=95, lty=1, col="red" ) +
  geom_hline(yintercept=90, lty=2, col="red" ) +
  coord_flip(ylim=c(80,100) ) + 
  labs( x="", y="coverage" )
ggsave( filename="../figures/coverage.pdf", width=8, height=5)

head( covsub )
summary( covsub$coverage )
sum( is.na( covsub$coverage ) )
summary( covsub$R )


# What is the estimation precision of the coverages?
SE.cov = mutate( covsub,
                 SE.cov = sqrt( coverage * (100-coverage) / R ) )
summary( SE.cov$SE.cov )


head( covsub )
nrow( covsub )
c2 = filter( covsub, dependence == "indep", proptx.dependence=="indep" )
nrow( c2 )
ggplot( c2, aes( method, coverage, col=population ) ) +
  facet_grid(  weight + population ~ J, scales="free_y" ) +
  geom_boxplot() +
  geom_hline(yintercept=95, lty=1, col="red" ) +
  geom_hline(yintercept=90, lty=2, col="red" ) +
  coord_flip(ylim=c(80,100) ) + 
  labs( x="", y="coverage" ) 
ggsave( filename="../figures/coverage_nobias.pdf", width=8, height=5)




##### Comparing FIRC to DB-SP-Sites #####

head( covsub )
m2_cov = covsub %>% 
  filter( method %in% c( "FIRC", "DB-SP-Sites" ) ) %>%
  dplyr::select( method, ID, J, n.bar, dependence, proptx.dependence, tau, coverage ) %>%
  pivot_wider( names_from=method, values_from=c( coverage ) )
head( m2_cov )


m2_cov = mutate( m2_cov, rel_cov = pmin( 0, 0.95 - `DB-SP-Sites` ) - pmin( 0, 0.95 - FIRC ) )
qplot( m2_cov$rel_cov, binwidth = 0.2 ) +
  labs( x = "Relative coverage of DB-SP-Sites vs FIRC (truncating for overcoverage)" )


# Relative coverage
ggplot( m2_cov, aes( tau, rel_cov, col=J ) ) +
  facet_grid( dependence ~ proptx.dependence ) +
  geom_jitter( width=0.01)+
  geom_hline( yintercept = 0 )


# Plot comparing FIRC coverage to DB-Site coverage across scenarios.
ggplot( m2_cov, aes( FIRC, `DB-SP-Sites`, col=J ) ) +
  facet_grid( dependence  ~ tau ) +
  geom_abline( slope=1, intercept=0, col="grey" ) +
  geom_hline( yintercept = 95 ) +
  geom_vline( xintercept = 95 ) +
  coord_fixed() +
  geom_point( alpha=0.75 ) +
  theme( legend.position="right", 
         legend.direction="vertical", legend.key.width=unit(1,"cm"),
         legend.box = "vertical", 
         panel.border=element_rect(size=0.3, fill=NA)) +
  labs( y="Coverage of DB-SP-Sites", x= "Coverage of FIRC" ) 
ggsave( filename="../figures/coverage_FIRC_v_DB.pdf", width=8, height=5)


