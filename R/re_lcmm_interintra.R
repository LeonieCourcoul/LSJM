re_lcmm_interintra <- function(param, nb.e.a, variability_inter_visit, variability_intra_visit, Sigma.re,
                         beta,
                         mu.inter, sigma.epsilon.inter, mu.intra,sigma.epsilon.intra,
                         len_visit_i,  X_base_i,  U_base_i,   y_i,  offset_ID_i

){
  #browser()
  b_re <- matrix(param[1:nb.e.a], nrow = 1)
  f_b_tau <- 1
  if(variability_inter_visit && variability_intra_visit){
    tau_re <- param[,(nb.e.a+1):(nb.e.a+2)]
    f_b_tau <- mvtnorm::dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
    sigma_inter <- exp(mu.inter + tau_re[1])
    var.inter <- sigma_inter**2
    sigma_intra <- exp(mu.intra + tau_re[2])
    var.intra <- sigma_intra**2
  }
  else{
    if(variability_inter_visit){
      tau_re <- param[,(nb.e.a+1)]
      f_b_tau <- mvtnorm::dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
      sigma_inter <- exp(mu.inter + tau_re[1])
      var.inter <- sigma_inter**2
      sigma_intra <- sigma.epsilon.intra
      var.intra <- sigma.epsilon.intra**2
    }
    else{
      if(variability_intra_visit){
        tau_re <- param[,(nb.e.a+1)]
        f_b_tau <- mvtnorm::dmvnorm(x = c(b_re, tau_re), mean = rep(0,length(b_re)+length(tau_re)), sigma = Sigma.re)
        sigma_intra <- exp(mu.intra + tau_re[2])
        var.intra <- sigma_intra**2
        sigma_inter <- sigma.epsilon.inter
        var.inter <- sigma.epsilon.inter**2
      }
      else{
        sigma_inter <- sigma.epsilon.inter
        var.inter <- sigma.epsilon.inter**2
        sigma_intra <- sigma.epsilon.intra
        var.intra <- sigma.epsilon.intra**2
      }

    }
  }

  sigma_inter_intra <- list(sigma_inter, sigma_intra, var.inter+var.intra, var.inter, var.intra, var.intra*(2*var.inter+var.intra))

  log_f_Y_f_T <- re_lcmm_interintra_cpp( beta,
                                      b_y=t(matrix(b_re, nrow = 1)),   sigma_inter_intra,
                                      len_visit_i,  X_base_i,  U_base_i,   y_i,  offset_ID_i
  )

  log_f_Y_f_T <- log_f_Y_f_T + log(f_b_tau)

  if(is.na(log_f_Y_f_T)){
    print(param)
    log_f_Y_f_T <- -1E09
  }
  log_f_Y_f_T



}
