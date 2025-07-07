log_llh_lsmm_interintra <- function(param,nb.beta, Zq,
                        nb.e.a, S,
                        variability_inter_visit, variability_intra_visit,
                        correlated_re, X_base, U_base, y_new, Ind, offset, offset_ID, offset_position, len_visit)
{


  #Manage parameter
  curseur <- 1
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  curseur <- curseur+nb.beta
  if(variability_inter_visit){
    mu.inter <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.inter <- param[curseur]
    curseur <- curseur +1
  }
  if(variability_intra_visit){
    mu.intra <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.intra <- param[curseur]
    curseur <- curseur +1
  }
  ### Cholesky matrix for random effects
  if(variability_inter_visit && variability_intra_visit){
    if(correlated_re){

      C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }
    else{
      borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      C3 <- matrix(rep(0,2*nb.e.a), ncol = nb.e.a)
      C4 <- matrix(rep(0,2*nb.e.a), nrow = nb.e.a)
      Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
      Cholesky <- as.matrix(Cholesky)
    }
  }
  else{
    if(variability_inter_visit){
      if(correlated_re){
        C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
        Cholesky <- C1
      }
      else{
        borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
        C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        C3 <- matrix(rep(0,nb.e.a), ncol = 1)
        C4 <- matrix(rep(0,nb.e.a), nrow = 1)
        Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
        Cholesky <- as.matrix(Cholesky)
      }
    }
    else{
      if(variability_intra_visit){
        if(correlated_re){
          C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
          Cholesky <- C1
        }
        else{
          borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
          C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          C3 <- matrix(rep(0,nb.e.a), ncol = 1)
          C4 <- matrix(rep(0,nb.e.a), nrow = 1)
          Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
          Cholesky <- as.matrix(Cholesky)
        }
      }
    }
    if(!variability_inter_visit && !variability_intra_visit){
      C1 <- matrix(rep(0,(length(param)-curseur)**2),nrow=length(param)-curseur,ncol=length(param)-curseur)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }

  }
  # Manage random effects
  random.effects <- Zq%*%t(Cholesky)
  b_y <- random.effects[,1:nb.e.a]
  b_y <- matrix(b_y, ncol = nb.e.a)
  if(variability_inter_visit && variability_intra_visit){
    b_inter <- random.effects[,nb.e.a+1]
    sigma.inter <- exp(mu.inter + b_inter)
    var.inter <- sigma.inter**2

    b_intra <- random.effects[,nb.e.a+2]
    sigma.intra <- exp(mu.intra + b_intra)
    var.intra <- sigma.intra**2
  }
  else{
    if(variability_inter_visit){
      b_inter <- random.effects[,nb.e.a+1]
      sigma.inter <- exp(mu.inter + b_inter)
      var.inter <- sigma.inter**2
    }
    else{
      if(variability_intra_visit){
        b_intra <- random.effects[,nb.e.a+1]
        sigma.intra <- exp(mu.intra + b_intra)
        var.intra <- sigma.intra**2
      }
    }
  }
  random.effects <- Zq%*%t(Cholesky)
  b_y <- random.effects[,1:nb.e.a]
  b_y <- matrix(b_y, ncol = nb.e.a)
  if(variability_inter_visit){
    b_inter <- random.effects[,nb.e.a+1]
    sigma_inter <- exp(mu.inter + b_inter)
    var.inter <- sigma_inter**2
  }
  else{
    sigma_inter <- rep(sigma.epsilon.inter,S)
    var.inter <- rep(sigma.epsilon.inter**2,S)
  }
  if(variability_intra_visit){
    b_intra <- random.effects[,ncol(random.effects)]
    sigma_intra <- exp(mu.intra + b_intra)
    var.intra <- sigma_intra**2
  }
  else{
    sigma_intra <- rep(sigma.epsilon.intra,S)
    var.intra <- rep(sigma.epsilon.intra**2,S)
  }
  sigma_inter_intra <- list(sigma_inter, sigma_intra, var.inter+var.intra, var.inter, var.intra, var.intra*(2*var.inter+var.intra))

  log_ll <- log_llh_lsmm_interintra_cpp(S,  sigma_inter_intra, X_base, U_base,
                                        y_new, offset_ID, beta, b_y,
                                        offset, offset_position,
                                        len_visit,
                                        Ind )


  ll_glob <- sum(log_ll)
  ll_glob
}
