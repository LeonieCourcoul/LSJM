re_lsjm_covDepIDMCase2 <- function(param, nb.e.a, nb.e.a.sigma,  Sigma.re,
                                    sharedtype, HB, W_G, nb_pointsGK,
                                    alpha_y_slope_var, alpha_z,  gamma_z0,  beta,  beta_slope, omega,  wk,
                                    delta2_i, Z_01_i, Z_02_i,  X_T_i,  U_T_i,
                                    Xslope_T_i,  Uslope_T_i, O_T_i,  W_T_i,
                                    X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                    Uslope_GK_T_i, O_GK_T_i,  W_GK_T_i,
                                    X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                    O_GK_T0_i,  W_GK_T0_i,
                                    Time_T_i,    Time_T0_i, st_T_i,    st_T0_i,
                                    B_T_i_02,
                                    Bs_T_i_01,  Bs_T_i_02,
                                    Bs_T0_i_01,  Bs_T0_i_02,  left_trunc,
                                    X_base_i,  U_base_i,   y_i, O_base_i,  W_base_i,  index_b_slope

){

  all_re <- matrix(param[1:(nb.e.a+nb.e.a.sigma)], nrow = 1)
  f_b_tau <- mvtnorm::dmvnorm(x = all_re, mean = rep(0,length(all_re)), sigma = Sigma.re)
  b_y_slope <- as.matrix(0)
  b_re <- all_re[1:nb.e.a]
  tau_re <- all_re[(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]
  if(sharedtype[2] || sharedtype[5] || sharedtype[8] ){
    b_y_slope <- as.matrix(b_re[index_b_slope], ncol = nb.e.a-1)
  }


  log_f_Y_f_T <- re_lsjm_covDepIDMCase2_cpp(  sharedtype,  HB,  W_G,
                                               nb_pointsGK ,
                                               alpha_y_slope_var,  alpha_z,  gamma_z0,  beta,  beta_slope, omega,
                                               b_y=t(matrix(b_re, nrow = 1)),  b_y_slope= t(matrix(b_y_slope, nrow = 1)), tau_re=t(matrix(tau_re, nrow = 1)),  wk,
                                               delta2_i,  Z_01_i,  Z_02_i,  X_T_i,  U_T_i,
                                               Xslope_T_i,  Uslope_T_i, O_T_i,  W_T_i,
                                              X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                               Uslope_GK_T_i, O_GK_T_i,  W_GK_T_i,
                                               X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i, O_GK_T0_i,  W_GK_T0_i,
                                               Time_T_i,   Time_T0_i, st_T_i,   st_T0_i,
                                               B_T_i_02,
                                               Bs_T_i_01,  Bs_T_i_02,
                                               Bs_T0_i_01,  Bs_T0_i_02,  left_trunc,
                                               X_base_i,  U_base_i,   y_i, O_base_i,  W_base_i
  )#50

  log_f_Y_f_T <- log_f_Y_f_T + log(f_b_tau)

  if(is.na(log_f_Y_f_T)){
    print(param)
    log_f_Y_f_T <- -1E09
  }
  log_f_Y_f_T

}
