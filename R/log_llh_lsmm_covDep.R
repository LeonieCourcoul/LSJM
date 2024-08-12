log_llh_lsmm_covDep <- function(param, nb.e.a, nb.beta, S,Zq, X_base, offset, U_base, y.new.prog,Ind,
                                nb.e.a.sigma, nb.omega ,
                                O_base , W_base,correlated_re){

  #Manage parameter
  curseur <- 1
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  curseur <- curseur+nb.beta
  omega <- param[(curseur):(curseur+nb.omega-1)]
  curseur <- curseur+nb.omega
  if(correlated_re){
    C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
    C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
    MatCov <- C1
    MatCov <- as.matrix(MatCov)
    random.effects <- Zq%*%t(MatCov)
    b_al <- random.effects[,1:nb.e.a]
    b_al <- matrix(b_al, ncol = nb.e.a)
    b_om <- random.effects[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]
    b_om <- matrix(b_om, ncol = nb.e.a.sigma)
  }
  else{
    borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
    C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
    C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
    MatCovb <- as.matrix(C1)
    MatCovSig <- as.matrix(C3)
    b_al <- Zq[,1:nb.e.a]%*%t(MatCovb)
    b_al <- matrix(b_al, ncol = nb.e.a)
    b_om <- Zq[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]%*%t(MatCovSig)
    b_om <- matrix(b_om, ncol = nb.e.a.sigma)
  }

  log_ll <- log_llh_lsmm_covDep_cpp(S, omega, b_om, beta, b_al,
                                    X_base, O_base, W_base, U_base, y.new.prog, Ind,
                                    offset)

  ll_glob <- sum(log_ll)
  ll_glob
}
