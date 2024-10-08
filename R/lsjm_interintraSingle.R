lsjm_interintraSingle <- function(Objectlsmm, Time, deltas, hazard_baseline_01, nb.knots.splines,
                              formSurv_01, nb_pointsGK, sharedtype_01, formSlopeFixed, formSlopeRandom,
                              index_beta_slope , index_b_slope,timeVar,
                              S1, S2, binit,nproc , clustertype, maxiter,
                              print.info, file , epsa , epsb , epsd){



  time.prog1 <- Sys.time()

  data.long <- as.data.frame(Objectlsmm$control$data.long)
  formGroup <- Objectlsmm$control$formGroup; formFixed <- Objectlsmm$control$formFixed; formRandom <- Objectlsmm$control$formRandom; correlated_re <- Objectlsmm$control$correlated_re
  #timeVar <- Objectlsmm$control$timeVar
  variability_inter_visit <- Objectlsmm$control$var_inter; variability_intra_visit <- Objectlsmm$control$var_intra; formGroupVisit <- Objectlsmm$control$formGroupVisit; Ind <- Objectlsmm$control$Ind
  Time_T <- Time[["Time_T"]]
  data.long$Time_T <- data.long[all.vars(Time_T)][,1]
  if(!is.null(Time[["Time_T0"]])){
    left_trunc <- TRUE
    Time_T0 <- Time[["Time_T0"]]
    data.long$Time_T0 <- data.long[all.vars(Time_T0)][,1]
  }
  else{
    left_trunc <- FALSE
  }

  delta1 <- deltas[["delta1"]]
  data.long$delta1 <- data.long[all.vars(delta1)][,1]
  knots_01 <- NULL

  ## Survival initialisation
  message("Survival initialisation")
  data.id <- data.long[!duplicated(data.long$id),]
  data.id <- as.data.frame(data.id)
  data.id$Time_T[which(data.id$Time_T == 0)] <- 1e-20
  if(hazard_baseline_01 == "Splines"){
    pp <- seq(0,1, length.out = nb.knots.splines[1]+2)
    pp <- utils::tail(utils::head(pp,-1),-1)
    vec.time.01 <- data.id[which(data.id$delta1==1),"Time_T"]
    kn <- quantile(vec.time.01, pp, names = FALSE)
    kn <- kn[kn<max(data.id$Time_T)]
    knots_01 <- sort(c(rep(range(data.id$Time_T,0), 4L), kn))
  }

  if(length(all.vars(formSurv_01))==0){
    Surv_01 <- Surv(Time_T, delta1) ~ 1
  }
  else{
    Surv_01 <- as.formula(paste("Surv(Time_T,delta1) ~ ", paste(all.vars(formSurv_01), collapse="+")))
  }

  if(hazard_baseline_01 == "Exponential" ){
    mod_surv <- survreg(Surv_01, data = data.id, dist = "exponential")
    alpha_01 <- mod_surv$coefficients
    alpha_01[1] <- -alpha_01[1]
  }
  else{
    if(hazard_baseline_01 == "Weibull"){
      mod_surv <- survreg(Surv_01, data = data.id, dist = "weibull")
      alpha_01 <- mod_surv$coefficients
      alpha_01[1] <- -mod_surv$coefficients[1]/sqrt(mod_surv$scale)
      shape_01 <- sqrt(1/mod_surv$scale)
    }
    else{
      if(hazard_baseline_01 == "Gompertz"){
        mod_surv <- flexsurvreg(Surv_01, data = data.id, dist = "gompertz")
        gompertz.1_01 <- sqrt(mod_surv$res[2,1])
        gompertz.2_01 <- mod_surv$res[1,1]
        alpha_01 <- NULL
        if(length(all.vars(formSurv_01)) >=1){
          alpha_01 <- mod_surv$res[3:(3+length(all.vars(formSurv_01))-1),1]
        }
      }
      else{
        if(hazard_baseline_01 == "Splines"){
          list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = nb_pointsGK)
          B <- splines::splineDesign(knots_01, data.id$Time_T, ord = 4L)
          Bs <- splines::splineDesign(knots_01, c(t(list.GK_T$st)), ord = 4L)
          opt_splines_01 <- optim(rep(0,ncol(B)), fn2,event = data.id$delta1,W2 = B,P = list.GK_T$P,wk = list.GK_T$wk,
                                  Time = data.id$Time_T,W2s = Bs,id.GK = list.GK_T$id.GK, method="BFGS", hessian = T)
          tmp_model <- coxph(Surv_01,
                             data = data.id,
                             x = TRUE)


          alpha_01 <- tmp_model$coefficients

        }
      }

    }


  }



  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); B_T_01 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Time_T0 = c(0)
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0);
  list.long <- data.manag.long(formGroup,formFixed, formRandom,data.long)
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  ID.visit <- data.long[all.vars(formGroupVisit)][,1]; offset <- list.long$offset

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

  list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = nb_pointsGK)
  st_T <- list.GK_T$st
  if(left_trunc){
    list.GK_T0 <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T0, k = nb_pointsGK)
    st_T0 <- list.GK_T0$st
    Time_T0 <- data.id$Time_T0
  }
  if(("value" %in% sharedtype_01)){
    list.data_T <- data.time(data.id, data.id$Time_T, formFixed, formRandom,timeVar)
    list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),formFixed, formRandom,timeVar)
    X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
    X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
    if(left_trunc){
      list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                   formFixed, formRandom,timeVar)
      X_GK_T0 <- list.data.GK_T0$Xtime
      U_GK_T0 <- list.data.GK_T0$Utime
    }
  }

  if(("slope" %in% sharedtype_01)){
    list.data_T <- data.time(data.id, data.id$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
    list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), formSlopeFixed, formSlopeRandom, timeVar)
    Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
    Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
    if(left_trunc){
      list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                   formSlopeFixed, formSlopeRandom,timeVar)
      Xslope_GK_T0 <- list.data.GK_T0$Xtime
      Uslope_GK_T0 <- list.data.GK_T0$Utime
    }
  }

  list.surv <- data.manag.surv(formGroup, formSurv_01, data.long)
  Z_01 <- list.surv$Z
  if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
  if(hazard_baseline_01 == "Splines"){
    Z_01 <- as.matrix(Z_01[,-1])
    B_T_01 <- splineDesign(knots_01, data.id$Time_T, ord = 4L)
    Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
    if(left_trunc){
      Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
    }
  }

  Matrices <- list( "delta1" = data.id$delta1,"Z_01" = Z_01, "Time_T" = data.id$Time_T,
                    "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                    "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                    "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01),
                    "Bs_T_01" = Bs_T_01,  "Time_T0" = Time_T0,
                    "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                    "Bs_T0_01" = Bs_T0_01, "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
  )

  nb.beta <-  ncol(X_base)
  nb.alpha <- c(ncol(Z_01))
  nb.e.a <- ncol(U_base)
  name_ZO1 <- colnames(Z_01)

  binit_CR <- c()
  names.param <- c()
  # 01
  if(hazard_baseline_01 == "Weibull"){
    binit_CR <- c(binit_CR, shape_01)
    names.param <- c(names.param, 'shape_01')
  }
  else{
    if(hazard_baseline_01 == "Gompertz"){
      binit_CR <- c(binit_CR, gompertz.1_01, gompertz.2_01)
      names.param <- c(names.param, 'gompertz.1_01', 'gompertz.2_01')
    }
    else{
      if(hazard_baseline_01 == "Splines"){
        binit_CR <- c(binit_CR,opt_splines_01$par)
        for(i in 1:length(opt_splines_01$par)){
          names.param <- c(names.param, paste("splines01", i, sep = "_"))
        }
      }

    }
  }
  binit_CR <- c(binit_CR, alpha_01)
  if(!is.null(alpha_01)){
    names.param <- c(names.param, paste(name_ZO1,"01",sep = "_"))
  }
  if("random effects" %in% sharedtype_01){
    binit_CR <- c(binit_CR, rep(0,nb.e.a))
    names.param <- c(names.param, paste("re",colnames(U_base),"01",sep = "_"))
  }
  if("value" %in% sharedtype_01){
    binit_CR <- c(binit_CR, 0)
    names.param <- c(names.param, 'value 01')
  }
  if("slope" %in% sharedtype_01){
    binit_CR <- c(binit_CR, 0)
    names.param <- c(names.param, 'slope 01')
  }
  if("variability inter" %in% sharedtype_01){
    binit_CR <- c(binit_CR,0)
    names.param <- c(names.param, 'variability inter 01')
  }
  if("variability intra" %in% sharedtype_01){
    binit_CR <- c(binit_CR,0)
    names.param <- c(names.param, 'variability intra 01')
  }
  if(is.null(Objectlsmm$result_step2)){
    binit_CR <- c(binit_CR, Objectlsmm$result_step1$b)
  }
  else{
    binit_CR <- c(binit_CR, Objectlsmm$result_step2$b)
  }

  if(variability_inter_visit && variability_intra_visit){
    Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a+2)
    Zq <- apply(Zq1, 2, qnorm)
  }
  else{
    if(variability_inter_visit || variability_intra_visit){
      Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a+1)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      Zq1 <- spacefillr::generate_sobol_owen_set(S1,  nb.e.a)
      Zq <- apply(Zq1, 2, qnorm)
    }
  }

  if(is.null(binit)){
    binit <- binit_CR
  }
  else{
    if(length(binit) != length(binit_CR)){
      stop("binit has not the correct number of arguments")
    }
  }

  message(paste("First estimation with ", S1, " QMC draws"))

  estimation1 <- marqLevAlg(binit, fn = logR_llh_lsjm_interintraSingle, minimize = FALSE,

                            hazard_baseline_01 = hazard_baseline_01, sharedtype_01 = sharedtype_01,
                            ord.splines = nb.knots.splines + 2, nb.beta = nb.beta, Zq = Zq, nb_pointsGK = nb_pointsGK,
                            nb.e.a = nb.e.a, S = S1, wk = gaussKronrod()$wk, rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk)), sk_GK = gaussKronrod()$sk, nb.alpha = nb.alpha,
                            variability_inter_visit = variability_inter_visit, variability_intra_visit = variability_intra_visit,
                            correlated_re = correlated_re, Matrices = Matrices, left_trunc = left_trunc,
                            knots.hazard_baseline.splines_01 = knots_01,
                            index_beta_slope = index_beta_slope, index_b_slope = index_b_slope, Ind = Ind,
                            nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                            file = file, blinding = FALSE, epsa = epsa, epsb = epsb, epsd = epsd)


  estimation2 <- NULL
  info_conv_step2 <- NULL

  if(!is.null(S2)){
    message(paste("Second estimation with ", S2, " QMC draws"))
    if(variability_inter_visit && variability_intra_visit){
      Zq1 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a+2)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      if(variability_inter_visit || variability_intra_visit){
        Zq1 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a+1)
        Zq <- apply(Zq1, 2, qnorm)
      }
      else{
        Zq1 <- spacefillr::generate_sobol_owen_set(S2,  nb.e.a)
        Zq <- apply(Zq1, 2, qnorm)
      }
    }

    estimation2 <- marqLevAlg(estimation1$b, fn = logR_llh_lsjm_interintraSingle, minimize = FALSE,

                              hazard_baseline_01 = hazard_baseline_01, sharedtype_01 = sharedtype_01,
                              ord.splines = nb.knots.splines + 2, nb.beta = nb.beta, Zq = Zq,
                              nb.e.a = nb.e.a, S = S2, wk = gaussKronrod()$wk, rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk)), sk_GK = gaussKronrod()$sk,nb_pointsGK = nb_pointsGK,
                              nb.alpha = nb.alpha,
                              variability_inter_visit = variability_inter_visit, variability_intra_visit = variability_intra_visit,
                              correlated_re = correlated_re, Matrices = Matrices, left_trunc = left_trunc,
                              knots.hazard_baseline.splines_01 = knots_01,
                              index_beta_slope = index_beta_slope, index_b_slope = index_b_slope, Ind = Ind,

                              nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                              file = file, blinding = FALSE, epsa = 10000000, epsb = 10000000, epsd = 0.99999)

    var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
    var_trans[upper.tri(var_trans, diag=T)] <- estimation2$v
    sd.param <- sqrt(diag(var_trans))
    param_est <-  estimation2$b

    #Delta-method
    nb.chol <- Objectlsmm$control$nb.chol
    if(correlated_re){
      if(variability_inter_visit && variability_intra_visit){
        curseur <- length(estimation2$b) - nb.chol + 1
        C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
        C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
        C1 <- as.matrix(C1)
        Index.C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
        Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+2,2)+nb.e.a+2)
        Index.C1 <- as.matrix(Index.C1)

        MatCov <- C1%*%t(C1)
        param_est <- c(param_est,unique(c(t(MatCov))))
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
        if(variability_inter_visit || variability_intra_visit){
          curseur <- length(estimation2$b) - nb.chol + 1
          C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- estimation2$b[curseur:length(estimation2$b)]
          C1 <- as.matrix(C1)
          Index.C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
          Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+1,2)+nb.e.a+1)
          Index.C1 <- as.matrix(Index.C1)

          MatCov <- C1%*%t(C1)
          param_est <- c(param_est,unique(c(t(MatCov))))
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

      if(variability_inter_visit && variability_intra_visit){
        borne3 <- borne1 + choose(n = 2, k = 2) + 2
        C3 <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
        C3[lower.tri(C3, diag=T)] <- estimation2$b[(borne1+1):borne3]
        C3 <- as.matrix(C3)

        Index.C3 <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
        Index.C3[lower.tri(Index.C3, diag=T)] <- 1:(choose(2,2)+2)
        Index.C3 <- as.matrix(Index.C3)

        MatCovSig <- C3%*%t(C3)
        param_est <- c(param_est,unique(c(t(MatCovSig))))


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

    info_conv_step2 <- list(conv = estimation2$istop, niter = estimation2$ni,
                            convcrit = c(estimation2$ca, estimation2$cb, estimation2$rdm))
  }
  else{
    var_trans <- matrix(rep(0,length(binit)**2),nrow=length(binit),ncol=length(binit))
    var_trans[upper.tri(var_trans, diag=T)] <- estimation1$v
    sd.param <- sqrt(diag(var_trans))
    param_est <-  estimation1$b

    #Delta-method
    nb.chol <- Objectlsmm$control$nb.chol
    if(correlated_re){
      if(variability_inter_visit && variability_intra_visit){
        curseur <- length(estimation1$b) - nb.chol + 1
        C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
        C1[lower.tri(C1, diag=T)] <- estimation1$b[curseur:length(estimation1$b)]
        C1 <- as.matrix(C1)
        Index.C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
        Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+2,2)+nb.e.a+2)
        Index.C1 <- as.matrix(Index.C1)

        MatCov <- C1%*%t(C1)
        param_est <- c(param_est,unique(c(t(MatCov))))
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
        if(variability_inter_visit || variability_intra_visit){
          curseur <- length(estimation1$b) - nb.chol + 1
          C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- estimation1$b[curseur:length(estimation1$b)]
          C1 <- as.matrix(C1)
          Index.C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
          Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a+1,2)+nb.e.a+1)
          Index.C1 <- as.matrix(Index.C1)

          MatCov <- C1%*%t(C1)
          param_est <- c(param_est,unique(c(t(MatCov))))
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
          C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
          C1[lower.tri(C1, diag=T)] <- estimation1$b[curseur:length(estimation1$b)]
          C1 <- as.matrix(C1)
          Index.C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
          Index.C1[lower.tri(Index.C1, diag=T)] <- 1:(choose(nb.e.a,2)+nb.e.a)
          Index.C1 <- as.matrix(Index.C1)

          MatCov <- C1%*%t(C1)
          param_est <- c(param_est,unique(c(t(MatCov))))
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
      MatCovb <- C1%*%t(C1)
      param_est <- c(param_est,unique(c(t(MatCovb))))


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

      if(variability_inter_visit && variability_intra_visit){
        borne3 <- borne1 + choose(n = 2, k = 2) + 2
        C3 <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
        C3[lower.tri(C3, diag=T)] <- estimation1$b[(borne1+1):borne3]
        C3 <- as.matrix(C3)

        Index.C3 <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
        Index.C3[lower.tri(Index.C3, diag=T)] <- 1:(choose(2,2)+2)
        Index.C3 <- as.matrix(Index.C3)

        MatCovSig <- C3%*%t(C3)
        param_est <- c(param_est,unique(c(t(MatCovSig))))


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
  }

  table.res <- cbind(param_est, sd.param)
  table.res <- as.data.frame(table.res)
  colnames(table.res) <- c("Estimation", "SE")
  rownames(table.res) <- c(names.param, rownames(Objectlsmm$table.res))

  time.prog2 <- Sys.time()
  time.prog.fin <- difftime(time.prog2, time.prog1)

  result <- list("table.res" = table.res,
                 "result_step1" = estimation1,
                 "result_step2" = estimation2,
                 "info_conv_step1" = list(conv = estimation1$istop, niter = estimation1$ni,
                                          convcrit = c(estimation1$ca, estimation1$cb, estimation1$rdm)),
                 "info_conv_step2" = info_conv_step2,
                 "time.computation" = time.prog.fin,
                 "control" = list(Objectlsmm = Objectlsmm,
                                  Time = Time,
                                  deltas = deltas,
                                  hazard_baseline_01 = hazard_baseline_01,
                                  nb.knots.splines = nb.knots.splines,
                                  formSurv_01 = formSurv_01,
                                  nb_pointsGK = nb_pointsGK,
                                  sharedtype_01 = sharedtype_01,
                                  formSlopeFixed = formSlopeFixed,
                                  formSlopeRandom = formSlopeRandom,
                                  index_b_slope = index_b_slope,
                                  index_beta_slope = index_beta_slope,
                                  knots.hazard_baseline.splines_01= knots_01,
                                  left_trunc = left_trunc,
                                  nb.alpha = nb.alpha,
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

  class(result) <- c("lsjm_interintraSingle")
  return(result)

}
