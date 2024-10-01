lsmm_covDep <- function(formFixed, formRandom, formGroup,
                        formFixedVar, formRandomVar,correlated_re ,
                        data.long, idVar, list.long, time.prog1,
                        S1 , S2,
                        nproc , clustertype, maxiter, print.info ,
                        file, epsa, epsb, epsd, binit){

  X_base <- list.long$X
  X_base <- as.matrix(X_base)
  U_base <- list.long$U
  U_base <- as.matrix(U_base)
  nb.e.a <- ncol(U_base)
  list.var <- data.manag.sigma(formGroup,formFixedVar, formRandomVar,data.long)
  O_base <- list.var$X
  O_base <- as.matrix(O_base)
  W_base <- list.var$U
  W_base <- as.matrix(W_base)
  nb.omega <- ncol(O_base)
  nb.e.a.sigma <- ncol(W_base)
  offset <- list.long$offset
  Ind <- list.long$I
  binit_initial <- binit
  list.init.long <- initial.long(formFixed, formRandom, idVar, data.long,
                                 ncol(X_base), nproc = nproc)
  sigma_epsilon <- list.init.long$sigma
  mu.log.sigma <- log(sigma_epsilon)
  cov_mat <- diag(nb.e.a)
  cholesky_b <- cov_mat[lower.tri(cov_mat, diag = T)]
  #cholesky_b <- list.init.long$long_model$cholesky
  priorMean.beta <- list.init.long$priorMean.beta
  names_param <- c()
  binit <- priorMean.beta
  names_param <- c(names_param, paste(colnames(X_base),"Y",sep = "_"))
  binit <- c(binit,mu.log.sigma,rep(0,nb.omega-1))
  names_param <- c(names_param, paste(colnames(O_base),"Var",sep = "_"))
  if(correlated_re){
    binit <- c(binit,
               cholesky_b,
               rep(0, nb.e.a*nb.e.a.sigma),
               rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma))
    nb.chol <- length(c(cholesky_b,
                        rep(0, nb.e.a*nb.e.a.sigma),
                        rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma)))
    for(i in 1:length(c(cholesky_b,
                        rep(0, nb.e.a*nb.e.a.sigma),
                        rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma)))){
      names_param <- c(names_param, paste("chol", i, sep = "_"))
    }
  }
  else{
    binit <- c(binit,
               cholesky_b,
               rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma))
    for(i in 1:length(cholesky_b)){
      names_param <- c(names_param, paste("chol_b", i, sep = "_"))
    }
    for(i in 1:length(rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma))){
      names_param <- c(names_param, paste("chol_tau", i, sep = "_"))
    }
    nb.chol <- length(c(cholesky_b,
                        rep(0, choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma)))
  }
  nb.beta <- length(priorMean.beta)

  Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a+nb.e.a.sigma)
  Zq <- apply(Zq1, 2, qnorm)

  if(!is.null(binit_initial)){
    if(length(binit_initial) == length(binit)){
      binit <- binit_initial
    }
    else{
      stop("binit has not the correct number of arguments")
    }
  }

  message(paste("First estimation with ", S1, " QMC draws"))
  estimation1 <- marqLevAlg(binit, fn = log_llh_lsmm_covDep, minimize = FALSE,
                           nb.e.a = nb.e.a, nb.beta = nb.beta, S = S1,Zq = Zq, X_base = X_base, offset = offset,
                           U_base = U_base, y.new.prog = list.long$y.new, Ind = Ind,
                           nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega,
                           O_base = O_base, W_base=W_base, correlated_re = correlated_re,
                           nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                           file = file, blinding = TRUE, epsa = epsa, epsb = epsb, epsd = epsd)



  estimation2 <- NULL
  info_conv_step2 <- NULL
  if(!is.null(S2)){
    Zq2 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a+nb.e.a.sigma)
    Zq <- apply(Zq2, 2, qnorm)


    message(paste("Second estimation with ", S2, " QMC draws"))
    estimation2 <- marqLevAlg(estimation1$b, fn = log_llh_lsmm_covDep, minimize = FALSE,
                              nb.e.a = nb.e.a, nb.beta = nb.beta, S = S2,Zq = Zq, X_base = X_base, offset = offset,
                              U_base = U_base, y.new.prog = list.long$y.new, Ind = Ind,
                              nb.e.a.sigma = nb.e.a.sigma, nb.omega = nb.omega,
                              O_base = O_base, W_base=W_base,correlated_re = correlated_re,
                              nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                              file = file, blinding = TRUE, epsa = 10000000, epsb = 10000000, epsd = 0.99999)




    var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
    var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
    sd.param <- sqrt(diag(var_trans))
    param_est <-  estimation2$b
    #Delta-method
    if(correlated_re){
      curseur <- length(estimation2$b) - nb.chol + 1
      C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
      C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
      C1 <- as.matrix(C1)
      Index.C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
      Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+nb.e.a.sigma,2)+nb.e.a+nb.e.a.sigma)
      Index.C1 <- as.matrix(Index.C1)

      MatCov <- C1%*%t(C1)
      param_est <- c(param_est,unique(c(t(MatCov))))
      vec_name1 <-  c(paste(colnames(U_base),"_Location",sep = "_"),paste(colnames(W_base),"_Scale",sep = "_"))
      mat_name1 <- outer(vec_name1,vec_name1, paste, sep = "*cov*")
      names_param <- c(names_param, c(mat_name1[lower.tri(mat_name1, diag= T)]))

      #for(i in 1:length(unique(c(t(MatCov))))){
      #  names_param <- c(names_param, paste("cov", i, sep = "_"))
      #}
      var_trans <- matrix(rep(0,length(estimation2$b)**2),nrow=length(estimation2$b),ncol=length(estimation2$b))
      var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
      trig.cov <- var_trans[curseur:length(estimation2$b),curseur:length(estimation2$b)]
      trig.cov <- trig.cov+t(trig.cov)
      diag(trig.cov) <- diag(trig.cov)/2
      cov.cholesky <- trig.cov
      Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))
      element.chol <- estimation2$b[curseur:length(estimation2$b)]
      for(i in 1:ncol(C1)){
        for(j in i:ncol(C1)){
          resultat <- 0
          k <- i
          m <- j
          for(t in 1:min(i,j,ncol(Cov.delta))){
            for(s in 1:min(k,m,ncol(Cov.delta))){
              resultat <- resultat +
                C1[j,t]*C1[m,s]*cov.cholesky[Index.C1[i,t],Index.C1[k,s]] +
                C1[j,t]*C1[k,s]*cov.cholesky[Index.C1[i,t],Index.C1[m,s]] +
                C1[i,t]*C1[m,s]*cov.cholesky[Index.C1[j,t],Index.C1[k,s]] +
                C1[i,t]*C1[k,s]*cov.cholesky[Index.C1[j,t],Index.C1[m,s]]
            }
          }
          sd.param <- c(sd.param,sqrt(resultat))
        }
      }
    }
    else{
      curseur <- length(estimation2$b) - nb.chol + 1
      borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:borne1]
      C1 <- as.matrix(C1)
      Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
      Index.C1 <- as.matrix(Index.C1)
      borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
      C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- estimation2$b[(borne1+1):borne3]
      C3 <- as.matrix(C3)

      Index.C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      Index.C3[lower.tri(Index.C3, diag=T)] <- 1:(choose(nb.e.a.sigma,2)+nb.e.a.sigma)
      Index.C3 <- as.matrix(Index.C3)

      MatCovb <- C1%*%t(C1)
      MatCovSig <- C3%*%t(C3)
      param_est <- c(param_est,unique(c(t(MatCovb))),unique(c(t(MatCovSig))))
      vec_name1 <-  paste(colnames(U_base),"Location",sep = "_")
      mat_name1 <- outer(vec_name1,vec_name1, paste, sep = "*cov*")
      vec_name2 <-  paste(colnames(W_base),"Scale",sep = "_")
      mat_name2 <- outer(vec_name2,vec_name2, paste, sep = "*cov*")
      names_param <- c(names_param, c(mat_name1[upper.tri(mat_name1, diag= T)]))
      names_param <- c(names_param, c(mat_name2[upper.tri(mat_name2, diag= T)]))


      #for(i in 1:length(unique(c(t(MatCovb))))){
      #  names_param <- c(names_param, paste("covB", i, sep = "_"))
      #}
      #for(i in 1:length(unique(c(t(MatCovSig))))){
      #  names_param <- c(names_param, paste("covSig", i, sep = "_"))
      #}
      var_trans <- matrix(rep(0,length(estimation2$b)**2),nrow=length(estimation2$b),ncol=length(estimation2$b))
      var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
      trig.cov <- var_trans[curseur:borne1,curseur:borne1]
      trig.cov <- trig.cov+t(trig.cov)
      diag(trig.cov) <- diag(trig.cov)/2
      cov.cholesky <- trig.cov
      Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))

      element.chol <- estimation2$b[curseur:length(estimation2$b)]
      for(i in 1:ncol(C1)){
        for(j in i:ncol(C1)){
          resultat <- 0
          k <- i
          m <- j
          for(t in 1:min(i,j,ncol(Cov.delta))){
            for(s in 1:min(k,m,ncol(Cov.delta))){
              resultat <- resultat +
                C1[j,t]*C1[m,s]*cov.cholesky[Index.C1[i,t],Index.C1[k,s]] +
                C1[j,t]*C1[k,s]*cov.cholesky[Index.C1[i,t],Index.C1[m,s]] +
                C1[i,t]*C1[m,s]*cov.cholesky[Index.C1[j,t],Index.C1[k,s]] +
                C1[i,t]*C1[k,s]*cov.cholesky[Index.C1[j,t],Index.C1[m,s]]
            }
          }
          sd.param <- c(sd.param,sqrt(resultat))
        }
      }
      trig.cov <- var_trans[(borne1+1):borne3,(borne1+1):borne3]
      trig.cov <- trig.cov+t(trig.cov)
      diag(trig.cov) <- diag(trig.cov)/2
      cov.cholesky <- trig.cov
      Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))

      element.chol <- estimation2$b[(borne1+1):borne3]
      for(i in 1:ncol(C3)){
        for(j in i:ncol(C3)){
          resultat <- 0
          k <- i
          m <- j
          for(t in 1:min(i,j,ncol(Cov.delta))){
            for(s in 1:min(k,m,ncol(Cov.delta))){
              resultat <- resultat +
                C3[j,t]*C3[m,s]*cov.cholesky[Index.C3[i,t],Index.C3[k,s]] +
                C3[j,t]*C3[k,s]*cov.cholesky[Index.C3[i,t],Index.C3[m,s]] +
                C3[i,t]*C3[m,s]*cov.cholesky[Index.C3[j,t],Index.C3[k,s]] +
                C3[i,t]*C3[k,s]*cov.cholesky[Index.C3[j,t],Index.C3[m,s]]
            }
          }
          sd.param <- c(sd.param,sqrt(resultat))
        }
      }
    }
  info_conv_step2 <-  list(conv = estimation2$istop, niter = estimation2$ni,
                         convcrit = c(estimation2$ca, estimation2$cb, estimation2$rdm))
  }
  else{
    var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
    var_trans[upper.tri(var_trans, diag=T)] <- estimation1$v
    sd.param <- sqrt(diag(var_trans))
    param_est <-  estimation1$b

    #Delta-method
    if(correlated_re){
      curseur <- length(estimation1$b) - nb.chol + 1
      C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
      C1[lower.tri(C1, diag=T)] <- estimation1$b[curseur:length(estimation1$b)]
      C1 <- as.matrix(C1)
      Index.C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
      Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+nb.e.a.sigma,2)+nb.e.a+nb.e.a.sigma)
      Index.C1 <- as.matrix(Index.C1)

      MatCov <- C1%*%t(C1)
      param_est <- c(param_est,unique(c(t(MatCov))))
      vec_name1 <-  c(paste(colnames(U_base),"_Location",sep = "_"),paste(colnames(W_base),"_Scale",sep = "_"))
      mat_name1 <- outer(vec_name1,vec_name1, paste, sep = "*cov*")
      names_param <- c(names_param, c(mat_name1[lower.tri(mat_name1, diag= T)]))

      var_trans <- matrix(rep(0,length(estimation1$b)**2),nrow=length(estimation1$b),ncol=length(estimation1$b))
      var_trans[upper.tri(var_trans, diag=T)] <- estimation1$v
      trig.cov <- var_trans[curseur:length(estimation1$b),curseur:length(estimation1$b)]
      trig.cov <- trig.cov+t(trig.cov)
      diag(trig.cov) <- diag(trig.cov)/2
      cov.cholesky <- trig.cov
      Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))
      element.chol <- estimation1$b[curseur:length(estimation1$b)]
      for(i in 1:ncol(C1)){
        for(j in i:ncol(C1)){
          resultat <- 0
          k <- i
          m <- j
          for(t in 1:min(i,j,ncol(Cov.delta))){
            for(s in 1:min(k,m,ncol(Cov.delta))){
              resultat <- resultat +
                C1[j,t]*C1[m,s]*cov.cholesky[Index.C1[i,t],Index.C1[k,s]] +
                C1[j,t]*C1[k,s]*cov.cholesky[Index.C1[i,t],Index.C1[m,s]] +
                C1[i,t]*C1[m,s]*cov.cholesky[Index.C1[j,t],Index.C1[k,s]] +
                C1[i,t]*C1[k,s]*cov.cholesky[Index.C1[j,t],Index.C1[m,s]]
            }
          }
          sd.param <- c(sd.param,sqrt(resultat))
        }
      }
    }
    else{
      curseur <- length(estimation1$b) - nb.chol + 1
      borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- estimation1$b[curseur:borne1]
      C1 <- as.matrix(C1)
      Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
      Index.C1 <- as.matrix(Index.C1)
      borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
      C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      C3[lower.tri(C3, diag=T)] <- estimation1$b[(borne1+1):borne3]
      C3 <- as.matrix(C3)

      Index.C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
      Index.C3[lower.tri(Index.C3, diag=T)] <- 1:(choose(nb.e.a.sigma,2)+nb.e.a.sigma)
      Index.C3 <- as.matrix(Index.C3)

      MatCovb <- C1%*%t(C1)
      MatCovSig <- C3%*%t(C3)
      param_est <- c(param_est,unique(c(t(MatCovb))),unique(c(t(MatCovSig))))
      vec_name1 <-  paste(colnames(U_base),"Location",sep = "_")
      mat_name1 <- outer(vec_name1,vec_name1, paste, sep = "*cov*")
      vec_name2 <-  paste(colnames(W_base),"Scale",sep = "_")
      mat_name2 <- outer(vec_name2,vec_name2, paste, sep = "*cov*")
      names_param <- c(names_param, c(mat_name1[upper.tri(mat_name1, diag= T)]))
      names_param <- c(names_param, c(mat_name2[upper.tri(mat_name2, diag= T)]))


      var_trans <- matrix(rep(0,length(estimation1$b)**2),nrow=length(estimation1$b),ncol=length(estimation1$b))
      var_trans[upper.tri(var_trans, diag=T)] <- estimation1$v
      trig.cov <- var_trans[curseur:borne1,curseur:borne1]
      trig.cov <- trig.cov+t(trig.cov)
      diag(trig.cov) <- diag(trig.cov)/2
      cov.cholesky <- trig.cov
      Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))

      element.chol <- estimation1$b[curseur:length(estimation1$b)]
      for(i in 1:ncol(C1)){
        for(j in i:ncol(C1)){
          resultat <- 0
          k <- i
          m <- j
          for(t in 1:min(i,j,ncol(Cov.delta))){
            for(s in 1:min(k,m,ncol(Cov.delta))){
              resultat <- resultat +
                C1[j,t]*C1[m,s]*cov.cholesky[Index.C1[i,t],Index.C1[k,s]] +
                C1[j,t]*C1[k,s]*cov.cholesky[Index.C1[i,t],Index.C1[m,s]] +
                C1[i,t]*C1[m,s]*cov.cholesky[Index.C1[j,t],Index.C1[k,s]] +
                C1[i,t]*C1[k,s]*cov.cholesky[Index.C1[j,t],Index.C1[m,s]]
            }
          }
          sd.param <- c(sd.param,sqrt(resultat))
        }
      }
      trig.cov <- var_trans[(borne1+1):borne3,(borne1+1):borne3]
      trig.cov <- trig.cov+t(trig.cov)
      diag(trig.cov) <- diag(trig.cov)/2
      cov.cholesky <- trig.cov
      Cov.delta <- matrix(NA, ncol = ncol(cov.cholesky), nrow = ncol(cov.cholesky))

      element.chol <- estimation1$b[(borne1+1):borne3]
      for(i in 1:ncol(C3)){
        for(j in i:ncol(C3)){
          resultat <- 0
          k <- i
          m <- j
          for(t in 1:min(i,j,ncol(Cov.delta))){
            for(s in 1:min(k,m,ncol(Cov.delta))){
              resultat <- resultat +
                C3[j,t]*C3[m,s]*cov.cholesky[Index.C3[i,t],Index.C3[k,s]] +
                C3[j,t]*C3[k,s]*cov.cholesky[Index.C3[i,t],Index.C3[m,s]] +
                C3[i,t]*C3[m,s]*cov.cholesky[Index.C3[j,t],Index.C3[k,s]] +
                C3[i,t]*C3[k,s]*cov.cholesky[Index.C3[j,t],Index.C3[m,s]]
            }
          }
          sd.param <- c(sd.param,sqrt(resultat))
        }
      }
    }
  }


  table.res <- cbind(param_est, sd.param)
  table.res <- as.data.frame(table.res)
  colnames(table.res) <- c("Estimation", "SE")
  rownames(table.res) <- names_param

  time.prog2 <- Sys.time()
  time.prog.fin <- difftime(time.prog2, time.prog1)

  result <- list("table.res" = table.res,
                 "result_step1" = estimation1,
                 "result_step2" = estimation2,
                 "info_conv_step1" = list(conv = estimation1$istop, niter = estimation1$ni,
                                          convcrit = c(estimation1$ca, estimation1$cb, estimation1$rdm)),
                 "info_conv_step2" = info_conv_step2,
                 "time.computation" = time.prog.fin,
                 "control" = list(formFixed = formFixed,
                                formRandom = formRandom,
                                formGroup = formGroup,
                                formVar = "cov-dependent",
                                formFixedVar = formFixedVar,
                                formRandomVar = formRandomVar,
                                correlated_re = correlated_re,
                                data.long = data.long,
                                idVar = idVar,
                                X_base = X_base,
                                U_base = U_base,
                                O_base = O_base,
                                W_base = W_base,
                                y.new.prog = list.long$y.new,
                                Ind = Ind,
                                nb.beta = nb.beta,
                                nb.e.a = nb.e.a,
                                nb.omega = nb.omega,
                                nb.e.a.sigma = nb.e.a.sigma,
                                nb.chol = nb.chol,
                                S1 = S1, S2 = S2,
                                nproc = nproc,
                                clustertype = clustertype,
                                maxiter = maxiter,
                                print.info = print.info,
                                file = file,
                                epsa = epsa,
                                epsb = epsb,
                                epsd = epsd
                                ))

  class(result) <- c("lsmm_covDep")
  return(result)
}
