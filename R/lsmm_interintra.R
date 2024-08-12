lsmm_interintra <- function(formFixed, formRandom, formGroup,
                var_inter, var_intra, formGroupVisit,correlated_re ,
                data.long, idVar, list.long, time.prog1,
                S1 , S2,
                nproc , clustertype, maxiter, print.info ,
                file, epsa, epsb, epsd, binit){

  data.long.court <- data.long[!duplicated(data.long[,c("id",all.vars(formGroupVisit))]),]
  binit_initial <- binit
  list.init.long <- initial.long(formFixed, formRandom, idVar, data.long.court,
                                 ncol(list.long$X), nproc = 1)
  sigma_epsilon <- list.init.long$sigma
  mu.log.sigma <- log(sigma_epsilon)
  cholesky_b <- list.init.long$long_model$cholesky
  priorMean.beta <- list.init.long$priorMean.beta
  X_base <- list.long$X; U_base <- list.long$U
  y.new <- list.long$y.new
  ID.visit <- data.long[all.vars(formGroupVisit)][,1]; offset <- list.long$offset
  Ind <- list.long$I

  offset_ID <- c()
  len_visit <- c(0)
  for(oo in 1:Ind){
    ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
    offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
    len_visit <- c(len_visit,length(unique(ID.visit_i)))
    frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
    offset_ID <- rbind(offset_ID, frame_offset_ID_i)
  }
  offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))



  nb.e.a <- ncol(U_base)
  nb.beta <-  length(priorMean.beta)
  binit <- c(priorMean.beta, mu.log.sigma/2,mu.log.sigma/2)
  names_param <- c()
  names_param <- c(names_param, paste(colnames(X_base),"Y",sep = "_"))
  if(var_inter){
    names_param <- c(names_param, "mu_inter")
  }
  else{
    names_param <- c(names_param, "sigma_inter")
  }
  if(var_intra){
    names_param <- c(names_param, "mu_intra")
  }
  else{
    names_param <- c(names_param, "sigma_intra")
  }
  if(var_inter && var_intra){
    if(correlated_re){
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- cholesky_b
      C2 <- matrix(c(0.1, 0,0, 0.1),nrow=2,ncol=2, byrow = TRUE)
      C3 <- matrix(rep(0,2*nb.e.a), ncol = nb.e.a)
      C4 <- matrix(rep(0,2*nb.e.a), nrow = nb.e.a)
      Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
      binit <- c(binit,Cholesky[lower.tri(Cholesky, diag = T)])
      for(i in 1:length(c(Cholesky[lower.tri(Cholesky, diag = T)]))){
        names_param <- c(names_param, paste("chol", i, sep = "_"))
      }
      nb.chol <- length(c(Cholesky[lower.tri(Cholesky, diag = T)]))
    }
    else{
      binit <- c(binit,cholesky_b,
                      0.1,0.05,0.1)
      for(i in 1:length(c(cholesky_b))){
        names_param <- c(names_param, paste("cholB", i, sep = "_"))
      }
      names_param <- c(names_param, "cholSig_1", "cholSig_2", "cholSig_3")
      nb.chol <- length(c(cholesky_b,
                   0.1,0.05,0.1))
    }
  }
  else{
    if(var_inter){
      if(correlated_re){
        C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
        C1[lower.tri(C1, diag=T)] <- cholesky_b
        C2 <- matrix(c(0.1),nrow=1,ncol=1, byrow = TRUE)
        C3 <- matrix(rep(0,nb.e.a), ncol = nb.e.a)
        C4 <- matrix(rep(0,nb.e.a), nrow = nb.e.a)
        Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
        binit <- c(binit,Cholesky[lower.tri(Cholesky, diag = T)])
        for(i in 1:length(c(Cholesky[lower.tri(Cholesky, diag = T)]))){
          names_param <- c(names_param, paste("chol", i, sep = "_"))
        }
        nb.chol <- length(Cholesky[lower.tri(Cholesky, diag = T)])
      }
      else{
        binit <- c(binit,cholesky_b,
                        0.1)
        nb.chol <- length(c(cholesky_b,0.1) )
        for(i in 1:length(c(cholesky_b))){
          names_param <- c(names_param, paste("cholB", i, sep = "_"))
        }
        names_param <- c(names_param, "cholInter")
      }
    }
    else{
      if(var_intra){
        if(correlated_re){
          C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
          C1[lower.tri(C1, diag=T)] <- cholesky_b
          C2 <- matrix(c(0.1),nrow=1,ncol=1, byrow = TRUE)
          C3 <- matrix(rep(0,nb.e.a), ncol = nb.e.a)
          C4 <- matrix(rep(0,nb.e.a), nrow = nb.e.a)
          Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
          binit <- c(binit,Cholesky[lower.tri(Cholesky, diag = T)])
          nb.chol <- length(Cholesky[lower.tri(Cholesky, diag = T)])
          for(i in 1:length(c(Cholesky[lower.tri(Cholesky, diag = T)]))){
            names_param <- c(names_param, paste("chol", i, sep = "_"))
          }
        }
        else{
          binit <- c(binit,cholesky_b,
                     0.1)
          nb.chol <- length(c(cholesky_b,0.1) )
          for(i in 1:length(c(cholesky_b))){
            names_param <- c(names_param, paste("cholB", i, sep = "_"))
          }
          names_param <- c(names_param, "cholIntra")
        }
      }
      else{
        binit <- c(binit,cholesky_b)
        nb.chol <- length(c(cholesky_b))
        for(i in 1:length(c(cholesky_b))){
          names_param <- c(names_param, paste("cholB", i, sep = "_"))
        }
      }

      }

  }

  if(var_inter && var_intra){
    Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a+2)
    Zq <- apply(Zq1, 2, qnorm)
  }
  else{
    if(var_inter || var_intra){
      Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a+1)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a)
      Zq <- apply(Zq1, 2, qnorm)
    }
  }

  if(!is.null(binit_initial)){
    if(length(binit_initial) == length(binit)){
      binit <- binit_initial
    }
    else{
      stop("binit has not the correct number of arguments")
    }
  }

  message("First estimation")
  estimation1 <- marqLevAlg(binit, fn = log_llh_lsmm_interintra, minimize = FALSE,
                           nb.beta = nb.beta, Zq=Zq,
                           nb.e.a = nb.e.a, S = S1,
                           variability_inter_visit = var_inter, variability_intra_visit = var_intra,
                           correlated_re = correlated_re, X_base = X_base, U_base = U_base, y_new = y.new,
                           Ind = Ind, offset = offset, offset_ID = offset_ID, offset_position = offset_position, len_visit = len_visit,
                           nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                           file = file, blinding = FALSE, epsa = epsa, epsb = epsb, epsd = epsd)


  if(var_inter && var_intra){
    Zq1 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a+2)
    Zq <- apply(Zq1, 2, qnorm)
  }
  else{
    if(var_inter || var_intra){
      Zq1 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a+1)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      Zq1 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a)
      Zq <- apply(Zq1, 2, qnorm)
    }
  }

  message("Second estimation")
  estimation2 <- marqLevAlg(estimation1$b, fn = log_llh_lsmm_interintra, minimize = FALSE,
                           nb.beta = nb.beta, Zq=Zq,
                           nb.e.a = nb.e.a, S = S2,
                           variability_inter_visit = var_inter, variability_intra_visit = var_intra,
                           correlated_re = correlated_re, X_base = X_base, U_base = U_base, y_new = y.new,
                           Ind = Ind, offset = offset, offset_ID = offset_ID, offset_position = offset_position, len_visit = len_visit,
                           nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                           file = file, blinding = FALSE, epsa = epsa, epsb = epsb, epsd = epsd)

  var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
  var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
  sd.param <- sqrt(diag(var_trans))
  param_est <-  estimation2$b

  #Delta-method

  if(correlated_re){
    if(var_inter && var_intra){
      curseur <- length(estimation2$b) - nb.chol + 1
      C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
      C1 <- as.matrix(C1)
      Index.C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
      Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+2,2)+nb.e.a+2)
      Index.C1 <- as.matrix(Index.C1)

      MatCov <- C1%*%t(C1)
      param_est <- c(param_est,unique(c(t(MatCov))))
      for(i in 1:length(unique(c(t(MatCov))))){
        names_param <- c(names_param, paste("cov", i, sep = "_"))
      }
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
      if(var_inter || var_intra){
        curseur <- length(estimation2$b) - nb.chol + 1
        C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
        C1 <- as.matrix(C1)
        Index.C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
        Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+1,2)+nb.e.a+1)
        Index.C1 <- as.matrix(Index.C1)

        MatCov <- C1%*%t(C1)
        param_est <- c(param_est,unique(c(t(MatCov))))
        for(i in 1:length(unique(c(t(MatCov))))){
          names_param <- c(names_param, paste("cov", i, sep = "_"))
        }
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
        C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
        C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
        C1 <- as.matrix(C1)
        Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
        Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
        Index.C1 <- as.matrix(Index.C1)

        MatCov <- C1%*%t(C1)
        param_est <- c(param_est,unique(c(t(MatCov))))
        for(i in 1:length(unique(c(t(MatCov))))){
          names_param <- c(names_param, paste("covB", i, sep = "_"))
        }
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
    MatCovb <- C1%*%t(C1)
    param_est <- c(param_est,unique(c(t(MatCovb))))
    for(i in 1:length(unique(c(t(MatCovb))))){
      names_param <- c(names_param, paste("covB", i, sep = "_"))
    }

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

    if(var_inter && var_intra){
      borne3 <- borne1 + choose(n = 2, k = 2) + 2
      C3 <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
      C3[lower.tri(C3, diag=T)] <- estimation2$b[(borne1+1):borne3]
      C3 <- as.matrix(C3)

      Index.C3 <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
      Index.C3[lower.tri(Index.C3, diag=T)] <- 1:(choose(2,2)+2)
      Index.C3 <- as.matrix(Index.C3)

      MatCovSig <- C3%*%t(C3)
      param_est <- c(param_est,unique(c(t(MatCovSig))))
      for(i in 1:length(unique(c(t(MatCovSig))))){
        names_param <- c(names_param, paste("covSig", i, sep = "_"))
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
                 "info_conv_step2" = list(conv = estimation2$istop, niter = estimation2$ni,
                                          convcrit = c(estimation2$ca, estimation2$cb, estimation2$rdm), time = time.prog.fin),
                 "control" = list(formFixed = formFixed,
                                  formRandom = formRandom,
                                  formGroup = formGroup,
                                  formVar = "inter-intra",
                                  var_inter = var_inter,
                                  var_intra = var_intra,
                                  correlated_re = correlated_re,
                                  formGroupVisit = formGroupVisit,
                                  data.long = data.long,
                                  data.long.court = data.long.court,
                                  idVar = idVar,
                                  X_base = X_base,
                                  U_base = U_base,
                                  y.new.prog = list.long$y.new,
                                  Ind = Ind,
                                  nb.beta = nb.beta,
                                  nb.e.a = nb.e.a,
                                  offset = offset,
                                  offset_ID = offset_ID,
                                  offset_position = offset_position,
                                  len_visit = len_visit,
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

  class(result) <- c("lsmm_interintra")
  return(result)
}


