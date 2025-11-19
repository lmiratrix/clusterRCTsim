


theoretical_SE_calc <- function( params = NULL, dta = NULL, control_formula = NULL ) {
  if ( !is.null( dta ) ) {
    params2 = describe_clusterRCT( Yobs ~ Z | S.id | D.id, data= dta,
                                   control_formula = control_formula )
    
    params2 = params2[ names(params2) %in% c("ICC.2", "ICC.3",
                                             "R2.2", "R2.1",
                                             "nbar", "Jbar", "K",
                                             "p.tx") ]
    nms = setdiff( names(params2), names(params) )
    params = c( params, params2[ nms ] )
  }
  
  pp = paste( names(params), params, sep=": ")
  cli::cli_inform( "Params = {pp}" )
  
  clv = with( params, ICC.2 * (1-R2.2) / (p.tx*(1-p.tx)*K*Jbar) )
  wcl = with( params, (1 - ICC.3 - ICC.2)*(1-R2.1) / (p.tx*(1-p.tx)*K*Jbar*nbar) )
  
  cli::cli_inform( "clv = {format(clv, digits = 3)}, wcl = {format(wcl, digits = 3)} -> per = {format(clv/(clv+wcl), digits=3)}" )
  VAR = with( params, clv + wcl )
  
  sqrt(VAR)
}




analyze_sim_data <- function( sim.data, inc_MLM = FALSE ) {
  
  # Group mean center our combined covariate
  sim.data <- sim.data %>%
    group_by( D.id, S.id ) %>%
    mutate( C_mn = mean(C),
            C_gc = C - C_mn ) %>%
    ungroup()
  
  inc_agg = TRUE
  # No covariates at all
  c1 <- clusterRCT::compare_methods( Yobs ~ Z | S.id | D.id, data=sim.data,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = inc_agg,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  # Just the o.g. individual covariate
  c2 <- clusterRCT::compare_methods( Yobs ~ Z | S.id | D.id, data=sim.data,
                                     control_formula = ~ X,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = inc_agg,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  # Just the o.g. cluster covariate
  c3 <- clusterRCT::compare_methods( Yobs ~ Z | S.id | D.id, data=sim.data,
                                     control_formula = ~ W,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = inc_agg,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  # Just "C"
  c4 <- clusterRCT::compare_methods( Yobs ~ Z | S.id | D.id, data=sim.data,
                                     control_formula = ~ C,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = inc_agg,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  # C, broken into parts
  c5 <- clusterRCT::compare_methods( Yobs ~ Z | S.id | D.id, data=sim.data,
                                     control_formula = ~ C_mn + C_gc,
                                     include_MLM = inc_MLM,
                                     include_LM = TRUE,
                                     include_agg = inc_agg,
                                     include_gee = FALSE,
                                     include_method_characteristics = FALSE )
  
  
  c1 = bind_rows( raw=c1, X = c2, W = c3, 
                  C = c4, Csep = c5, 
                  .id="adjustment" )
  
  c1$tau = mean( sim.data$Y1 - sim.data$Y0 )
  c1$sdY0 = sd( sim.data$Y0 )
  
  c1
}