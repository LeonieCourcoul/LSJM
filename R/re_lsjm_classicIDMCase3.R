#' @importFrom mvtnorm dmvnorm

re_lsjm_classicIDMCase3 <- function(param, nb.e.a, Sigma.re,
                                       sharedtype, HB, Gompertz, Weibull, nb_pointsGK,
                                       alpha_y_slope, alpha_b_01, alpha_b_02, alpha_b_12, alpha_z,  gamma_z0,  beta,  beta_slope,  wk, rep_wk,
                                       sigma_epsilon,
                                       delta2_i,  Z_01_i,  Z_02_i,  Z_12_i,  X_T_i,  U_T_i,
                                       Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                       Uslope_GK_T_i,  X_GK_L_T_i,  U_GK_L_T_i,  Xslope_GK_L_T_i,  Uslope_GK_L_T_i,
                                       X_GK_0_LT_i,  U_GK_0_LT_i,  Xslope_GK_0_LT_i,  Uslope_GK_0_LT_i,
                                       X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                       Time_T_i,  Time_L_T_i,  Time_T0_i, st_T_i,  st_0_LT_i,  st_L_T_i,  st_T0_i,
                                       ck,
                                       B_T_i_02,  B_T_i_12,
                                       Bs_T_i_01,  Bs_T_i_02,  Bs_T_i_12,
                                       Bs_0_LT_i_01,  Bs_0_LT_i_02,  Bs_0_LT_i_12,
                                       Bs_L_T_i_01,
                                       Bs_T0_i_01,  Bs_T0_i_02,  left_trunc,
                                       X_base_i,  U_base_i,   y_i,  index_b_slope
){
  b_re <- matrix(param[1:nb.e.a], nrow = 1)
  b_y_slope <- as.matrix(0)
  if(sharedtype[2] || sharedtype[4] || sharedtype[6] ){
    b_y_slope <- as.matrix(b_re[index_b_slope], ncol = nb.e.a-1)
  }
  f_b_tau <- dmvnorm(x = c(b_re), mean = rep(0,length(b_re)), sigma = Sigma.re)


  log_f_Y_f_T <- re_lsjm_classicIDMCase3_cpp(sharedtype, HB, Gompertz, Weibull,
                                                nb_pointsGK,
                                                alpha_y_slope,alpha_b_01, alpha_b_02, alpha_b_12, alpha_z,  gamma_z0,  beta,  beta_slope,b_y = t(matrix(b_re, nrow = 1)),
                                                b_y_slope= t(matrix(b_y_slope, nrow = 1)),  wk, rep_wk, sigma_epsilon,
                                                delta2_i,  Z_01_i,  Z_02_i,  Z_12_i,  X_T_i,  U_T_i,
                                                Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                                Uslope_GK_T_i,  X_GK_L_T_i,  U_GK_L_T_i,  Xslope_GK_L_T_i,  Uslope_GK_L_T_i,
                                                X_GK_0_LT_i,  U_GK_0_LT_i,  Xslope_GK_0_LT_i,  Uslope_GK_0_LT_i,
                                                X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                                Time_T_i,  Time_L_T_i,  Time_T0_i, st_T_i,  st_0_LT_i,  st_L_T_i,  st_T0_i,
                                                ck,
                                                B_T_i_02,  B_T_i_12,
                                                Bs_T_i_01,  Bs_T_i_02,  Bs_T_i_12,
                                                Bs_0_LT_i_01,  Bs_0_LT_i_02,  Bs_0_LT_i_12,
                                                Bs_L_T_i_01,
                                                Bs_T0_i_01,  Bs_T0_i_02,  left_trunc,
                                                X_base_i,  U_base_i,   y_i)

  log_f_Y_f_T <- log_f_Y_f_T + log(f_b_tau)

  if(is.na(log_f_Y_f_T)){
    print(param)
    log_f_Y_f_T <- -1E09
  }
  log_f_Y_f_T

}
