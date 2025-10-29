#' @importFrom mvtnorm dmvnorm
re_lsjm_interintraIDMCase2 <- function(param, nb.e.a, variability_inter_visit, variability_intra_visit, Sigma.re,
                          sharedtype, HB, W_G, nb_pointsGK,
                          alpha_y_slope, alpha_inter_intra,alpha_b_01, alpha_b_02, alpha_z,  gamma_z0,  beta,  beta_slope,  wk,
                          mu.inter, sigma.epsilon.inter, mu.intra,sigma.epsilon.intra,
                          delta2_i, Z_01_i, Z_02_i,  X_T_i,  U_T_i,
                          Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                          Uslope_GK_T_i,
                          X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                          Time_T_i,    Time_T0_i, st_T_i,    st_T0_i,
                          B_T_i_02,
                          Bs_T_i_01,  Bs_T_i_02,
                          Bs_T0_i_01,  Bs_T0_i_02,  left_trunc,
                          len_visit_i,  X_base_i,  U_base_i,   y_i,  offset_ID_i, index_b_slope

){

  b_re <- matrix(param[1:nb.e.a], nrow = 1)
  b_y_slope <- as.matrix(0)
  if(sharedtype[2] || sharedtype[6] || sharedtype[10] ){
    b_y_slope <- as.matrix(b_re[index_b_slope], ncol = nb.e.a-1)
  }
  if(variability_inter_visit && variability_intra_visit){
    tau_re <- param[,(nb.e.a+1):(nb.e.a+2)]
    f_b_tau <- dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
    sigma_inter <- exp(mu.inter + tau_re[1])
    var.inter <- sigma_inter**2
    sigma_intra <- exp(mu.intra + tau_re[2])
    var.intra <- sigma_intra**2
  }
  else{
    if(variability_inter_visit){
      tau_re <- param[,(nb.e.a+1)]
      f_b_tau <- dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
      sigma_inter <- exp(mu.inter + tau_re[1])
      var.inter <- sigma_inter**2
      sigma_intra <- sigma.epsilon.intra
      var.intra <- sigma_intra**2
    }
    else{
      if(variability_intra_visit){
        tau_re <- param[,(nb.e.a+1)]
        f_b_tau <- dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
        sigma_inter <- sigma.epsilon.inter
        var.inter <- sigma_inter**2
        sigma_intra <- exp(mu.intra + tau_re[1])
        var.intra <- sigma_intra**2
      }
      else{
        f_b_tau <- dmvnorm(x = c(b_re), mean = rep(0,length(b_re)), sigma = Sigma.re)
        sigma_inter <- sigma.epsilon.inter
        var.inter <- sigma_inter**2
        sigma_intra <- sigma.epsilon.intra
        var.intra <- sigma_intra**2
      }
    }
  }
  sigma_inter_intra <- list(sigma_inter, sigma_intra, var.inter+var.intra, var.inter, var.intra, var.intra*(2*var.inter+var.intra))

  log_f_Y_f_T <- re_lsjm_interintraIDMCase2_cpp(  sharedtype,  HB,  W_G,
                                        nb_pointsGK ,  alpha_inter_intra,
                                        alpha_y_slope, alpha_b_01, alpha_b_02,  alpha_z,  gamma_z0,  beta,  beta_slope,
                                        b_y=t(matrix(b_re, nrow = 1)),  b_y_slope= t(matrix(b_y_slope, nrow = 1)),  wk,  sigma_inter_intra,
                                        delta2_i,  Z_01_i,  Z_02_i,  X_T_i,  U_T_i,
                                        Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                        Uslope_GK_T_i,
                                        X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                        Time_T_i,   Time_T0_i, st_T_i,   st_T0_i,
                                        B_T_i_02,
                                        Bs_T_i_01,  Bs_T_i_02,
                                        Bs_T0_i_01,  Bs_T0_i_02,  left_trunc,
                                        len_visit_i,  X_base_i,  U_base_i,   y_i,  offset_ID_i
  )

  log_f_Y_f_T <- log_f_Y_f_T + log(f_b_tau)

  if(is.na(log_f_Y_f_T)){
    print(param)
    log_f_Y_f_T <- -1E09
  }
  log_f_Y_f_T

}
