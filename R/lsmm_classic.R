lsmm_classic <- function(formFixed, formRandom, formGroup,
                        data.long, idVar, list.long, time.prog1,
                        S1 , S2,
                        nproc , clustertype, maxiter, print.info ,
                        file, epsa, epsb, epsd, binit){
  #browser()
  X_base <- list.long$X
  X_base <- as.matrix(X_base)
  U_base <- list.long$U
  U_base <- as.matrix(U_base)
  nb.e.a <- ncol(U_base)
  offset <- list.long$offset
  Ind <- list.long$I
  binit_initial <- binit
  list.init.long <- initial.long(formFixed, formRandom, idVar, data.long,
                                 ncol(X_base), nproc = nproc)
  sigma_epsilon <- list.init.long$sigma
  cov_mat <- diag(nb.e.a)
  cholesky_b <- cov_mat[lower.tri(cov_mat, diag = T)]
  priorMean.beta <- list.init.long$priorMean.beta
  names_param <- c()
  binit <- priorMean.beta
  names_param <- c(names_param, paste(colnames(X_base),"Y",sep = "_"))
  binit <- c(binit, sigma_epsilon)
  names_param <- c(names_param, "sigma")
  binit <- c(binit,
             cholesky_b)
  for(i in 1:length(cholesky_b)){
    names_param <- c(names_param, paste("chol_b", i, sep = "_"))
  }
  nb.chol <- length(cholesky_b)
  nb.beta <- length(priorMean.beta)

  Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a)
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
  estimation1 <- marqLevAlg::marqLevAlg(binit, fn = log_llh_lsmm_classic, minimize = FALSE,
                           nb.e.a = nb.e.a, nb.beta = nb.beta, S = S1,Zq = Zq, X_base = X_base, offset = offset,
                           U_base = U_base, y.new.prog = list.long$y.new, Ind = Ind,
                           nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                           file = file, blinding = TRUE, epsa = epsa, epsb = epsb, epsd = epsd)

  estimation2 <- NULL
  info_conv_step2 <- NULL

  if(!is.null(S2)){
    Zq2 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a)
    Zq <- apply(Zq2, 2, qnorm)

    message(paste("Second estimation with ", S2, " QMC draws"))
    estimation2 <- marqLevAlg::marqLevAlg(estimation1$b, fn = log_llh_lsmm_classic, minimize = FALSE,
                                          nb.e.a = nb.e.a, nb.beta = nb.beta, S = S2,Zq = Zq, X_base = X_base, offset = offset,
                                          U_base = U_base, y.new.prog = list.long$y.new, Ind = Ind,
                                          nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                                          file = file, blinding = TRUE, epsa = 10000000, epsb = 10000000, epsd = 0.99999)




    var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
    var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
    sd.param <- sqrt(diag(var_trans))
    param_est <-  estimation2$b

    #Delta-method
    curseur <- length(estimation2$b) - nb.chol + 1
    C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
    C1 <- as.matrix(C1)

    Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
    Index.C1 <- as.matrix(Index.C1)

    MatCov <- C1%*%t(C1)
    param_est <- c(param_est,unique(c(t(MatCov))))
    vec_name1 <-  paste(colnames(U_base),"Location",sep = "_")
    mat_name1 <- outer(vec_name1,vec_name1, paste, sep = "*cov*")
    names_param <- c(names_param, c(mat_name1[upper.tri(mat_name1, diag= T)]))
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

    info_conv_step2 <- list(conv = estimation2$istop, niter = estimation2$ni,
                            convcrit = c(estimation2$ca, estimation2$cb, estimation2$rdm))
  }
  else{
    var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
    var_trans[upper.tri(var_trans, diag=T)] <- estimation1$v
    sd.param <- sqrt(diag(var_trans))
    param_est <-  estimation1$b

    #Delta-method
    curseur <- length(estimation1$b) - nb.chol + 1
    C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    C1[lower.tri(C1, diag=T)] <- estimation1$b[curseur:length(estimation1$b)]
    C1 <- as.matrix(C1)

    Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
    Index.C1 <- as.matrix(Index.C1)

    MatCov <- C1%*%t(C1)
    param_est <- c(param_est,unique(c(t(MatCov))))
    vec_name1 <-  paste(colnames(U_base),"Location",sep = "_")
    mat_name1 <- outer(vec_name1,vec_name1, paste, sep = "*cov*")
    names_param <- c(names_param, c(mat_name1[upper.tri(mat_name1, diag= T)]))

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
                                  formVar = "classic",
                                  data.long = data.long,
                                  idVar = idVar,
                                  X_base = X_base,
                                  U_base = U_base,
                                  y.new.prog = list.long$y.new,
                                  Ind = Ind,
                                  nb.beta = nb.beta,
                                  nb.e.a = nb.e.a,
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

  class(result) <- c("lsmm_classic")
  return(result)

}
