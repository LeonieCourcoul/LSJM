#' @importFrom mvtnorm dmvnorm
re_lsmm_classic <- function(param, nb.e.a, Sigma.re,beta, X_base_i, U_base_i,y_i, sigma_epsilon){
  all_re <- matrix(param[1:(nb.e.a)], nrow = 1)
  f_b_tau <- dmvnorm(x = all_re, mean = rep(0,length(all_re)), sigma = Sigma.re)



  log_f_Y_f_T <- re_lsmm_classic_cpp( beta, b_y=all_re,
                                     X_base_i,  U_base_i,  y_i, sigma_epsilon
  )

  log_f_Y_f_T <- log_f_Y_f_T + log(f_b_tau)

  if(is.na(log_f_Y_f_T)){
    print(param)
    log_f_Y_f_T <- -1E09
  }
  log_f_Y_f_T



}
