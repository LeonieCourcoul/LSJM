log_llh_lsmm_classic <- function(param, nb.e.a, nb.beta, S,Zq, X_base, offset, U_base, y.new.prog,Ind){
  #browser()
  #Manage parameter
  curseur <- 1
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  curseur <- curseur+nb.beta
  sigma_epsilon <- param[curseur]
  curseur <- curseur +1
  C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
  MatCov <- C1
  MatCov <- as.matrix(MatCov)
  random.effects <- Zq%*%t(MatCov)
  b_al <- random.effects[,1:nb.e.a]
  b_al <- matrix(b_al, ncol = nb.e.a)


  log_ll <- log_llh_lsmm_classic_cpp( S,  beta, b_al,
                                    X_base,  U_base, y.new.prog, Ind,
                                    offset, sigma_epsilon)

  ll_glob <- sum(log_ll)
  ll_glob
}
