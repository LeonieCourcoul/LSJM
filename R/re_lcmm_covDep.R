re_lcmm_covDep <- function(param, nb.e.a, nb.e.a.sigma, Sigma.re,beta, omega, X_base_i, U_base_i, O_base_i, W_base_i, y_i){
  all_re <- matrix(param[1:(nb.e.a+nb.e.a.sigma)], nrow = 1)
  f_b_tau <- mvtnorm::dmvnorm(x = all_re, mean = rep(0,length(all_re)), sigma = Sigma.re)



  log_f_Y_f_T <- re_lcmm_covDep_cpp( beta, b_y=matrix(all_re[1:nb.e.a]), omega,  b_om=matrix(all_re[(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]),
                                     X_base_i,  U_base_i, O_base_i, W_base_i,  y_i
  )

  log_f_Y_f_T <- log_f_Y_f_T + log(f_b_tau)

  if(is.na(log_f_Y_f_T)){
    print(param)
    log_f_Y_f_T <- -1E09
  }
  log_f_Y_f_T



}
