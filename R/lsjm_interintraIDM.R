#' @importFrom stats quantile qnorm as.formula optim
#' @importFrom utils head tail
#' @importFrom survival Surv survreg coxph
#' @importFrom flexsurv flexsurvreg
#' @importFrom splines splineDesign
#' @importFrom spacefillr generate_sobol_owen_set
#' @importFrom marqLevAlg marqLevAlg
lsjm_interintraIDM <- function(Objectlsmm, Time, deltas, hazard_baseline_01, hazard_baseline_02, hazard_baseline_12, nb.knots.splines,
                               formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                               formSlopeFixed, formSlopeRandom,
                               index_beta_slope , index_b_slope,timeVar,
                               S1, S2, binit,nproc , clustertype, maxiter,
                               print.info, file , epsa , epsb , epsd){

  time.prog1 <- Sys.time()

  data.long <- as.data.frame(Objectlsmm$control$data.long)
  formGroup <- Objectlsmm$control$formGroup; formFixed <- Objectlsmm$control$formFixed; formRandom <- Objectlsmm$control$formRandom; correlated_re <- Objectlsmm$control$correlated_re
  variability_inter_visit <- Objectlsmm$control$var_inter; variability_intra_visit <- Objectlsmm$control$var_intra; formGroupVisit <- Objectlsmm$control$formGroupVisit; Ind <- Objectlsmm$control$Ind
 # timeVar <- Objectlsmm$control$timeVar
  Time_T <- Time[["Time_T"]]
  data.long$Time_T <- data.long[all.vars(Time_T)][,1]
  Time_R <- Time[["Time_R"]]
  data.long$Time_R <- data.long[all.vars(Time_R)][,1]
  Time_L <- Time[["Time_L"]]
  data.long$Time_L <- data.long[all.vars(Time_L)][,1]
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
  delta2 <- deltas[["delta2"]]
  data.long$delta2 <- data.long[all.vars(delta2)][,1]
  knots_01 <- NULL
  knots_12 <- NULL
  knots_02 <- NULL

  ## Sub-groups creation for IDM
  data.long$Case <- ifelse(data.long$delta1 == 1 & data.long$Time_R == data.long$Time_L, "Case1bis",
                           ifelse(data.long$delta1 == 1, "Case1",
                                  ifelse(data.long$Time_L == data.long$Time_T, "Case2", "Case3")))
  id.Case1 <- unique(data.long$id[which(data.long$Case == "Case1")])
  id.Case1bis <- unique(data.long$id[which(data.long$Case == "Case1bis")])
  id.Case2 <- unique(data.long$id[which(data.long$Case == "Case2")])
  id.Case3 <- unique(data.long$id[which(data.long$Case == "Case3")])
  nbCase1 <- length(id.Case1); nbCase1bis <- length(id.Case1bis); nbCase2 <- length(id.Case2); nbCase3 <- length(id.Case3)

  ## Survival initialisation
  message("Survival initialisation")

  data.id <- data.long[!duplicated(data.long$id),]
  data.id <- as.data.frame(data.id)
  data.id$Time_R[which(data.id$Time_R == 0)] <- 1e-20
  data.id$Time_L[which(data.id$Time_L == 0)] <- 1e-20
  data.id$Time_T[which(data.id$Time_T == 0)] <- 1e-20

  if(hazard_baseline_01 == "Splines"){
    pp <- seq(0,1, length.out = nb.knots.splines[1]+2)
    pp <- tail(head(pp,-1),-1)
    vec.time.01 <- data.id[which(data.id$delta1==1),"Time_R"]
    kn <- quantile(vec.time.01, pp, names = FALSE)
    kn <- kn[kn<max(data.id$Time_T)]
    knots_01 <- sort(c(rep(range(data.id$Time_T,0), 4L), kn))
  }
  if(hazard_baseline_02 == "Splines"){
    pp <- seq(0,1, length.out =nb.knots.splines[2]+2)
    pp <- tail(head(pp,-1),-1)
    vec.time.02 <- data.id[which(data.id$delta1==0 & data.id$delta2==1),"Time_T"]
    kn <- quantile(vec.time.02, pp, names = FALSE)
    kn <- kn[kn<max(data.id$Time_T)]
    knots_02 <- sort(c(rep(range(data.id$Time_T,0), 4L), kn))
  }
  if(hazard_baseline_12 == "Splines"){
    pp <- seq(0,1, length.out = nb.knots.splines[3]+2)
    pp <- tail(head(pp,-1),-1)
    vec.time.12 <- data.id[which(data.id$delta1==1 & data.id$delta2==1), "Time_T"]
    kn <- quantile(vec.time.12, pp, names = FALSE)
    kn <- kn[kn<max(data.id$Time_T)]
    knots_12 <- sort(c(rep(range(data.id$Time_T,0), 4L), kn))
  }

  if(length(all.vars(formSurv_01))==0){
    Surv_01 <- Surv(Time_R, delta1) ~ 1
  }
  else{
    Surv_01 <- as.formula(paste("Surv(Time_R,delta1) ~ ", paste(all.vars(formSurv_01), collapse="+")))
  }
  if(length(all.vars(formSurv_02))==0){
    Surv_02 <- Surv(Time_T, delta2) ~ 1
  }
  else{
    Surv_02 <- as.formula(paste("Surv(Time_T,delta2) ~ ", paste(all.vars(formSurv_02), collapse="+")))
  }
  if(length(all.vars(formSurv_12))==0){
    Surv_12 <- Surv(Time_T, delta2) ~ 1
  }
  else{
    Surv_12 <- as.formula(paste("Surv(Time_T,delta2) ~ ", paste(all.vars(formSurv_12), collapse="+")))
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
          B <- splineDesign(knots_01, data.id$Time_T, ord = 4L)
          Bs <- splineDesign(knots_01, c(t(list.GK_T$st)), ord = 4L)
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

  if(hazard_baseline_02 == "Exponential" ){
    mod_surv <- survreg(Surv_02, data = data.id, dist = "exponential")
    alpha_02 <- mod_surv$coefficients
    alpha_02[1] <- -alpha_02[1]
  }
  else{
    if(hazard_baseline_02 == "Weibull"){
      mod_surv <- survreg(Surv_02, data = data.id, dist = "weibull")
      alpha_02 <- mod_surv$coefficients
      alpha_02[1] <- -mod_surv$coefficients[1]/sqrt(mod_surv$scale)
      shape_02 <- sqrt(1/mod_surv$scale)
    }
    else{
      if(hazard_baseline_02 == "Gompertz"){
        mod_surv <- flexsurvreg(Surv_02, data = data.id, dist = "gompertz")
        gompertz.1_02 <- sqrt(mod_surv$res[2,1])
        gompertz.2_02 <- mod_surv$res[1,1]
        alpha_02 <- NULL
        if(length(all.vars(formSurv_02)) >=1){
          alpha_02 <- mod_surv$res[3:(3+length(all.vars(formSurv_02))-1),1]
        }
      }
      else{
        if(hazard_baseline_02 == "Splines"){
          list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = nb_pointsGK)
          B <- splineDesign(knots_02, data.id$Time_T, ord = 4L)
          Bs <- splineDesign(knots_02, c(t(list.GK_T$st)), ord = 4L)
          opt_splines_02 <- optim(rep(0,ncol(B)), fn2,event = data.id$delta2,W2 = B,P = list.GK_T$P,wk = list.GK_T$wk,
                                  Time = data.id$Time_T,W2s = Bs,id.GK = list.GK_T$id.GK, method="BFGS", hessian = T)
          tmp_model <- coxph(Surv_02,
                             data = data.id,
                             x = TRUE)


          alpha_02 <- tmp_model$coefficients

        }
      }
    }
  }

  if(hazard_baseline_12 == "Exponential" ){
    mod_surv <- survreg(Surv_12, data = data.id, dist = "exponential")
    alpha_12 <- mod_surv$coefficients
    alpha_12[1] <- -alpha_12[1]
  }
  else{
    if(hazard_baseline_12 == "Weibull"){
      mod_surv <- survreg(Surv_12, data = data.id, dist = "weibull")
      alpha_12 <- mod_surv$coefficients
      alpha_12[1] <- -mod_surv$coefficients[1]/sqrt(mod_surv$scale)
      shape_12 <- sqrt(1/mod_surv$scale)
    }
    else{
      if(hazard_baseline_12 == "Gompertz"){
        mod_surv <- flexsurvreg(Surv_12, data = data.id, dist = "gompertz")
        gompertz.1_12 <- sqrt(mod_surv$res[2,1])
        gompertz.2_12 <- mod_surv$res[1,1]
        alpha_12 <- NULL
        if(length(all.vars(formSurv_12)) >=1){
          alpha_12 <- mod_surv$res[3:(3+length(all.vars(formSurv_02))-1),1]
        }
      }
      else{
        if(hazard_baseline_12 == "Splines"){
          list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = nb_pointsGK)
          B <- splineDesign(knots_12, data.id$Time_T, ord = 4L)
          Bs <- splineDesign(knots_12, c(t(list.GK_T$st)), ord = 4L)
          opt_splines_12 <- optim(rep(0,ncol(B)), fn2,event = data.id$delta2,W2 = B,P = list.GK_T$P,wk = list.GK_T$wk,
                                  Time = data.id$Time_T,W2s = Bs,id.GK = list.GK_T$id.GK, method="BFGS", hessian = T)
          tmp_model <- coxph(Surv_12,
                             data = data.id,
                             x = TRUE)


          alpha_12 <- tmp_model$coefficients

        }
      }
    }
  }


  binit_CI <- c()
  # 01
  if(hazard_baseline_01 == "Weibull"){
    binit_CI <- c(binit_CI, shape_01)
  }
  else{
    if(hazard_baseline_01 == "Gompertz"){
      binit_CI <- c(binit_CI, gompertz.1_01, gompertz.2_01)
    }
    else{
      if(hazard_baseline_01 == "Splines"){
        binit_CI <- c(binit_CI,opt_splines_01$par)
      }
    }
  }
  binit_CI <- c(binit_CI, alpha_01)
  if("random effects" %in% sharedtype_01){
    binit_CI <- c(binit_CI, rep(0,nb.e.a))
  }
  if("value" %in% sharedtype_01){
    binit_CI <- c(binit_CI, 0)
  }
  if("slope" %in% sharedtype_01){
    binit_CI <- c(binit_CI, 0)
  }
  if("variability inter" %in% sharedtype_01){
    binit_CI <- c(binit_CI,0)
  }
  if("variability intra" %in% sharedtype_01){
    binit_CI <- c(binit_CI,0)
  }
  # 02
  if(hazard_baseline_02 == "Weibull"){
    binit_CI <- c(binit_CI, shape_02)
  }
  else{
    if(hazard_baseline_02 == "Gompertz"){
      binit_CI <- c(binit_CI, gompertz.1_02, gompertz.2_02)
    }
    else{
      if(hazard_baseline_02 == "Splines"){
        binit_CI <- c(binit_CI,opt_splines_02$par)
      }

    }
  }
  binit_CI <- c(binit_CI, alpha_02)
  if("random effects" %in% sharedtype_02){
    binit_CI <- c(binit_CI, rep(0,nb.e.a))
  }
  if("value" %in% sharedtype_02){
    binit_CI <- c(binit_CI, 0)
  }
  if("slope" %in% sharedtype_02){
    binit_CI <- c(binit_CI, 0)
  }
  if("variability inter" %in% sharedtype_02){
    binit_CI <- c(binit_CI,0)
  }
  if("variability intra" %in% sharedtype_02){
    binit_CI <- c(binit_CI,0)
  }
  # 12
  if(hazard_baseline_12 == "Weibull"){
    binit_CI <- c(binit_CI, shape_12)
  }
  else{
    if(hazard_baseline_12 == "Gompertz"){
      binit_CI <- c(binit_CI, gompertz.1_12, gompertz.2_12)
    }
    else{
      if(hazard_baseline_12 == "Splines"){
        binit_CI <- c(binit_CI,opt_splines_12$par)
      }

    }
  }
  binit_CI <- c(binit_CI, alpha_12)
  if("random effects" %in% sharedtype_12){
    binit_CI <- c(binit_CI, rep(0,nb.e.a))
  }
  if("value" %in% sharedtype_12){
    binit_CI <- c(binit_CI, 0)
  }
  if("slope" %in% sharedtype_12){
    binit_CI <- c(binit_CI, 0)
  }
  if("variability inter" %in% sharedtype_12){
    binit_CI <- c(binit_CI,0)
  }
  if("variability intra" %in% sharedtype_12){
    binit_CI <- c(binit_CI,0)
  }

  if(is.null(Objectlsmm$result_step2)){
    binit_CI <- c(binit_CI, Objectlsmm$result_step1$b)
  }
  else{
    binit_CI <- c(binit_CI, Objectlsmm$result_step2$b)
  }



  # Initialisation avec le mileu de l'intervalle
  estimation.noCI <- NULL
  if(is.null(binit) &&  ((nbCase1 + nbCase3)/(nbCase1 + nbCase1bis + nbCase2 + nbCase3) > 0.10)){#nombre de cas 1 et 3 supérieur à 10% => seuil à discuter
    data.long_noCI <- data.long
    data.long_noCI$Time_L_initnoCI <- NA
    data.long_noCI$Time_R_initnoCI <- NA

    ## On calcule le temps du milieu d'intervalle pour les personnes devenues démentes et on prends le temps de décès ou de censure pour les autres
    data.long_noCI$Time_L_initnoCI[which(data.long_noCI$delta1 == 1)] <- (data.long_noCI$Time_L[which(data.long_noCI$delta1 == 1)] + data.long_noCI$Time_R[which(data.long_noCI$delta1 == 1)])/2
    data.long_noCI$Time_R_initnoCI[which(data.long_noCI$delta1 == 1)] <- (data.long_noCI$Time_L[which(data.long_noCI$delta1 == 1)] + data.long_noCI$Time_R[which(data.long_noCI$delta1 == 1)])/2
    data.long_noCI$Time_L_initnoCI[which(data.long_noCI$delta1 == 0)] <- data.long_noCI$Time_T[which(data.long_noCI$delta1 == 0)]
    data.long_noCI$Time_R_initnoCI[which(data.long_noCI$delta1 == 0)] <- data.long_noCI$Time_T[which(data.long_noCI$delta1 == 0)]
    #On créer les nouveaux sous-groupe pour l'initialisation sans CI
    data.long_noCI$Case <- ifelse(data.long_noCI$delta1 == 1 & data.long_noCI$Time_R_initnoCI == data.long_noCI$Time_L_initnoCI, "Case1bis",
                                  ifelse(data.long_noCI$delta1 == 1, "Case1",
                                         ifelse(data.long_noCI$Time_L_initnoCI == data.long_noCI$Time_T, "Case2", "Case3")))
    id.Case1_noCI <- unique(data.long_noCI$id[which(data.long_noCI$Case == "Case1")])
    id.Case1bis_noCI <- unique(data.long_noCI$id[which(data.long_noCI$Case == "Case1bis")])
    id.Case2_noCI <- unique(data.long_noCI$id[which(data.long_noCI$Case == "Case2")])
    id.Case3_noCI <- unique(data.long_noCI$id[which(data.long_noCI$Case == "Case3")])
    nbCase1_noCI <- length(id.Case1_noCI); nbCase1bis_noCI <- length(id.Case1bis_noCI); nbCase2_noCI <- length(id.Case2_noCI); nbCase3_noCI <- length(id.Case3_noCI)
    if(nbCase1_noCI + nbCase3_noCI!=0){
      stop("Something is wrong.")
    }

    data.id <- data.long_noCI[!duplicated(data.long_noCI$id),]
    data.id <- as.data.frame(data.id)
    data.id$Time_R_initnoCI[which(data.id$Time_R_initnoCI == 0)] <- 1e-20
    data.id$Time_L_initnoCI[which(data.id$Time_L_initnoCI == 0)] <- 1e-20
    data.id$Time_T[which(data.id$Time_T == 0)] <- 1e-20

    # Matrices of data for Case1bis
    Case1bis <- NULL
    st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
    X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); B_T_01 = as.matrix(0); B_T_02 = as.matrix(0); B_T_12 = as.matrix(0);
    Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L = as.matrix(0); U_GK_L = as.matrix(0);
    Xslope_GK_L = as.matrix(0); Uslope_GK_L = as.matrix(0);
    st_L = as.matrix(0); Bs_L_01 = as.matrix(0); Bs_L_02 = as.matrix(0); Bs_L_12 = as.matrix(0); B_L_01 = as.matrix(0); B_L_02 = as.matrix(0); B_L_12 = as.matrix(0);
    X_L = as.matrix(0); U_L = as.matrix(0); Xslope_L = as.matrix(0); Uslope_L = as.matrix(0); Time_T0 = c(0)
    st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
    Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)
    if(length(id.Case1bis_noCI)>0){
      data.long.Case1bis <- data.long_noCI[which(data.long_noCI$id %in% id.Case1bis_noCI),]
      data.id.Case1bis <- data.long.Case1bis[!duplicated(data.long.Case1bis$id),]
      list.long.Case1bis <- data.manag.long(formGroup,formFixed, formRandom,data.long.Case1bis)
      X_base <- list.long.Case1bis$X; U_base <- list.long.Case1bis$U; y.new <- list.long.Case1bis$y.new
      ID.visit <- data.long.Case1bis[all.vars(formGroupVisit)][,1]; offset <- list.long.Case1bis$offset

      offset_ID <- c()
      len_visit <- c(0)
      for(oo in 1:nbCase1bis_noCI){
        ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
        offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
        len_visit <- c(len_visit,length(unique(ID.visit_i)))
        frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
        offset_ID <- rbind(offset_ID, frame_offset_ID_i)
      }
      offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))

      list.GK_T <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_T, k = nb_pointsGK)
      list.GK_L <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_L_initnoCI, k = nb_pointsGK)
      st_T <- list.GK_T$st
      st_L <- list.GK_L$st
      if(left_trunc){
        list.GK_T0 <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_T0, k = nb_pointsGK)
        st_T0 <- list.GK_T0$st
        Time_T0 <- data.id.Case1bis$Time_T0
      }

      if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
        list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, formFixed, formRandom,timeVar)
        list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L_initnoCI, formFixed, formRandom,timeVar)
        list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),formFixed, formRandom,timeVar)
        list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)),formFixed, formRandom,timeVar)
        X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
        X_L <- list.data_L$Xtime; U_L <- list.data_L$Utime
        X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
        X_GK_L <- list.data.GK_L$Xtime; U_GK_L <- list.data.GK_L$Utime

        if(left_trunc){
          list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                       formFixed, formRandom,timeVar)
          X_GK_T0 <- list.data.GK_T0$Xtime
          U_GK_T0 <- list.data.GK_T0$Utime
        }
      }

      if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
        list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
        list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L_initnoCI, formSlopeFixed, formSlopeRandom, timeVar)
        list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), formSlopeFixed, formSlopeRandom, timeVar)
        list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)), formSlopeFixed, formSlopeRandom, timeVar)
        Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
        Xslope_L <- list.data_L$Xtime; Uslope_L <- list.data_L$Utime
        Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
        Xslope_GK_L <- list.data.GK_L$Xtime; Uslope_GK_L <- list.data.GK_L$Utime

        if(left_trunc){
          list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                       formSlopeFixed, formSlopeRandom,timeVar)
          Xslope_GK_T0 <- list.data.GK_T0$Xtime
          Uslope_GK_T0 <- list.data.GK_T0$Utime
        }
      }

      list.surv <- data.manag.surv(formGroup, formSurv_01, data.long.Case1bis)
      Z_01 <- list.surv$Z
      if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
      list.surv <- data.manag.surv(formGroup, formSurv_02, data.long.Case1bis)
      Z_02 <- list.surv$Z
      if(hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
      list.surv <- data.manag.surv(formGroup, formSurv_12, data.long.Case1bis)
      Z_12 <- list.surv$Z
      if(hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
      if(hazard_baseline_01 == "Splines"){
        Z_01 <- as.matrix(Z_01[,-1])
        B_T_01 <- splineDesign(knots_01, data.id.Case1bis$Time_T, ord = 4L)
        B_L_01 <- splineDesign(knots_01, data.id.Case1bis$Time_L_initnoCI, ord = 4L)
        Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
        Bs_L_01 <- splineDesign(knots_01, c(t(st_L)), ord = 4L)
        if(left_trunc){
          Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
        }
      }
      if(hazard_baseline_02 == "Splines"){
        Z_02 <- as.matrix(Z_02[,-1])
        B_T_02 <- splineDesign(knots_02, data.id.Case1bis$Time_T, ord = 4L)
        B_L_02 <- splineDesign(knots_02, data.id.Case1bis$Time_L_initnoCI, ord = 4L)
        Bs_T_02 <- splineDesign(knots_02, c(t(st_T)), ord = 4L)
        Bs_L_02 <- splineDesign(knots_02, c(t(st_L)), ord = 4L)
        if(left_trunc){
          Bs_T0_02 <- splineDesign(knots_02, c(t(st_T0)), ord = 4L)
        }
      }
      if(hazard_baseline_12 == "Splines"){
        Z_12 <- as.matrix(Z_12[,-1])
        B_T_12 <- splineDesign(knots_12, data.id.Case1bis$Time_T, ord = 4L)
        B_L_12 <- splineDesign(knots_12, data.id.Case1bis$Time_L_initnoCI, ord = 4L)
        Bs_T_12 <- splineDesign(knots_12, c(t(st_T)), ord = 4L)
        Bs_L_12 <- splineDesign(knots_12, c(t(st_L)), ord = 4L)
      }

      Case1bis <- list( "delta2" = data.id.Case1bis$delta2, "Z_01" = Z_01, "Z_02" = Z_02, "Z_12" = Z_12, "Time_T" = data.id.Case1bis$Time_T,
                        "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                        "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                        "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01), "B_T_02" = t(B_T_02), "B_T_12" = t(B_T_12),
                        "Bs_T_01" = Bs_T_01, "Bs_T_02" = Bs_T_02, "Bs_T_12" = Bs_T_12, "X_GK_L" = X_GK_L, "U_GK_L" = U_GK_L,
                        "Xslope_GK_L" = Xslope_GK_L, "Uslope_GK_L" = Uslope_GK_L, "Time_L" = data.id.Case1bis$Time_L_initnoCI,
                        "st_L" = t(st_L), "Bs_L_01" = Bs_L_01, "Bs_L_02" = Bs_L_02, "Bs_L_12" = Bs_L_12, "B_L_01" = t(B_L_01), "B_L_02" = t(B_L_02), "B_L_12" = t(B_L_12),
                        "X_L" = X_L, "U_L" = t(U_L), "Xslope_L" = Xslope_L, "Uslope_L" = t(Uslope_L), "Time_T0" = Time_T0,
                        "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                        "Bs_T0_01" = Bs_T0_01, "Bs_T0_02" = Bs_T0_02, "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
                        #  ,"random_effects" = re_Case1bis, "var_random_effects" = var_re_Case1bis
      )

      name_ZO1 <- colnames(Z_01)
      name_ZO2 <- colnames(Z_02)
      name_Z12 <- colnames(Z_12)
    }

    Case2 <- NULL
    st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
    X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
    B_T_01 = as.matrix(0); B_T_02 = as.matrix(0);
    Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Time_T0 = c(0);
    st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
    Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)
    if(length(id.Case2_noCI)>0){
      data.long.Case2 <- data.long_noCI[which(data.long_noCI$id %in% id.Case2_noCI),]
      data.id.Case2 <- data.long.Case2[!duplicated(data.long.Case2$id),]
      list.long.Case2 <- data.manag.long(formGroup,formFixed, formRandom,data.long.Case2)
      X_base <- list.long.Case2$X; U_base <- list.long.Case2$U; y.new <- list.long.Case2$y.new
      ID.visit <- data.long.Case2[all.vars(formGroupVisit)][,1]; offset <- list.long.Case2$offset
      offset_ID <- c()
      len_visit <- c(0)
      for(oo in 1:nbCase2_noCI){
        ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
        offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
        len_visit <- c(len_visit,length(unique(ID.visit_i)))
        frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
        offset_ID <- rbind(offset_ID, frame_offset_ID_i)
      }
      offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))

      list.GK_T <- data.GaussKronrod(data.id.Case2, a = 0, b = data.id.Case2$Time_T, k = nb_pointsGK)
      st_T <- list.GK_T$st
      if(left_trunc){
        list.GK_T0 <- data.GaussKronrod(data.id.Case2, a = 0, b = data.id.Case2$Time_T0, k = nb_pointsGK)
        st_T0 <- list.GK_T0$st
        Time_T0 <- data.id.Case2$Time_T0
      }
      if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) ){
        list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, formFixed, formRandom,timeVar)
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

      if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) ){
        list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
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

      list.surv <- data.manag.surv(formGroup, formSurv_01, data.long.Case2)
      Z_01 <- list.surv$Z
      if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
      list.surv <- data.manag.surv(formGroup, formSurv_02, data.long.Case2)
      Z_02 <- list.surv$Z
      if(hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
      if(hazard_baseline_01 == "Splines"){
        Z_01 <- as.matrix(Z_01[,-1])
        B_T_01 <- splineDesign(knots_01, data.id.Case2$Time_T, ord = 4L)
        Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
        if(left_trunc){
          Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
        }
      }
      if(hazard_baseline_02 == "Splines"){
        Z_02 <- as.matrix(Z_02[,-1])
        B_T_02 <- splineDesign(knots_02, data.id.Case2$Time_T, ord = 4L)
        Bs_T_02 <- splineDesign(knots_02, c(t(st_T)), ord = 4L)
        if(left_trunc){
          Bs_T0_02 <- splineDesign(knots_02, c(t(st_T0)), ord = 4L)
        }
      }
      Case2 <- list( "delta2" = data.id.Case2$delta2, "Z_01" = Z_01, "Z_02" = Z_02, "Z_12" = Z_12, "Time_T" = data.id.Case2$Time_T,
                     "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                     "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                     "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01), "B_T_02" = t(B_T_02),
                     "Bs_T_01" = Bs_T_01, "Bs_T_02" = Bs_T_02, "Time_T0" = Time_T0,
                     "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                     "Bs_T0_01" = Bs_T0_01, "Bs_T0_02" = Bs_T0_02, "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
                     #,"random_effects" = re_Case2, "var_random_effects" = var_re_Case2
      )
      name_ZO1 <- colnames(Z_01)
      name_ZO2 <- colnames(Z_02)
      name_Z12 <- colnames(Z_12)

    }

    nb.beta <-  ncol(X_base)
    nb.alpha <- c(ncol(Z_01), ncol(Z_02), ncol(Z_12))
    nb.e.a <- ncol(U_base)

    binit_noCI <- c()
    names.param <- c()
    # 01
    if(hazard_baseline_01 == "Weibull"){
      binit_noCI <- c(binit_noCI, shape_01)
      names.param <- c(names.param, 'shape_01')
    }
    else{
      if(hazard_baseline_01 == "Gompertz"){
        binit_noCI <- c(binit_noCI, gompertz.1_01, gompertz.2_01)
        names.param <- c(names.param, 'gompertz.1_01', 'gompertz.2_01')
      }
      else{
        if(hazard_baseline_01 == "Splines"){
          binit_noCI <- c(binit_noCI,opt_splines_01$par)
          for(i in 1:length(opt_splines_01$par)){
            names.param <- c(names.param, paste("splines01", i, sep = "_"))
          }
        }

      }
    }
    binit_noCI <- c(binit_noCI, alpha_01)
    if(!is.null(alpha_01)){
      names.param <- c(names.param, paste(name_ZO1,"",sep = "_"))
    }
    if("random effects" %in% sharedtype_01){
      binit_noCI <- c(binit_noCI, rep(0,nb.e.a))
      names.param <- c(names.param, paste("re",colnames(U_base),"01",sep = "_"))
    }
    if("value" %in% sharedtype_01){
      binit_noCI <- c(binit_noCI, 0)
      names.param <- c(names.param, 'value 01')
    }
    if("slope" %in% sharedtype_01){
      binit_noCI <- c(binit_noCI, 0)
      names.param <- c(names.param, 'slope 01')
    }
    if("variability inter" %in% sharedtype_01){
      binit_noCI <- c(binit_noCI,0)
      names.param <- c(names.param, 'variability inter 01')
    }
    if("variability intra" %in% sharedtype_01){
      binit_noCI <- c(binit_noCI,0)
      names.param <- c(names.param, 'variability intra 01')
    }
    # 02
    if(hazard_baseline_02 == "Weibull"){
      binit_noCI <- c(binit_noCI, shape_02)
      names.param <- c(names.param, 'shape_02')
    }
    else{
      if(hazard_baseline_02 == "Gompertz"){
        binit_noCI <- c(binit_noCI, gompertz.1_02, gompertz.2_02)
        names.param <- c(names.param, 'gompertz.1_02', 'gompertz.2_02')
      }
      else{
        if(hazard_baseline_02 == "Splines"){
          binit_noCI <- c(binit_noCI,opt_splines_02$par)
          for(i in 1:length(opt_splines_02$par)){
            names.param <- c(names.param, paste("splines02", i, sep = "_"))
          }
        }

      }
    }
    binit_noCI <- c(binit_noCI, alpha_02)
    if(!is.null(alpha_02)){
      names.param <- c(names.param, paste(name_ZO2,"",sep = "_"))
    }
    if("random effects" %in% sharedtype_02){
      binit_noCI <- c(binit_noCI, rep(0,nb.e.a))
      names.param <- c(names.param, paste("re",colnames(U_base),"02",sep = "_"))
    }
    if("value" %in% sharedtype_02){
      binit_noCI <- c(binit_noCI, 0)
      names.param <- c(names.param, 'value 02')
    }
    if("slope" %in% sharedtype_02){
      binit_noCI <- c(binit_noCI, 0)
      names.param <- c(names.param, 'slope 02')
    }
    if("variability inter" %in% sharedtype_02){
      binit_noCI <- c(binit_noCI,0)
      names.param <- c(names.param, 'variability inter 02')
    }
    if("variability intra" %in% sharedtype_02){
      binit_noCI <- c(binit_noCI,0)
      names.param <- c(names.param, 'variability intra 02')
    }
    # 12
    if(hazard_baseline_12 == "Weibull"){
      binit_noCI <- c(binit_noCI, shape_12)
      names.param <- c(names.param, 'shape_12')
    }
    else{
      if(hazard_baseline_12 == "Gompertz"){
        binit_noCI <- c(binit_noCI, gompertz.1_12, gompertz.2_12)
        names.param <- c(names.param, 'gompertz.1_12', 'gompertz.2_12')
      }
      else{
        if(hazard_baseline_12 == "Splines"){
          binit_noCI <- c(binit_noCI,opt_splines_12$par)
          for(i in 1:length(opt_splines_12$par)){
            names.param <- c(names.param, paste("splines12", i, sep = "_"))
          }
        }

      }
    }
    binit_noCI <- c(binit_noCI, alpha_12)
    if(!is.null(alpha_12)){
      names.param <- c(names.param, paste(name_Z12,"",sep = "_"))
    }
    if("random effects" %in% sharedtype_12){
      binit_noCI <- c(binit_noCI, rep(0,nb.e.a))
      names.param <- c(names.param, paste("re",colnames(U_base),"12",sep = "_"))
    }
    if("value" %in% sharedtype_12){
      binit_noCI <- c(binit_noCI, 0)
      names.param <- c(names.param, 'value 12')
    }
    if("slope" %in% sharedtype_12){
      binit_noCI <- c(binit_noCI, 0)
      names.param <- c(names.param, 'slope 12')
    }
    if("variability inter" %in% sharedtype_12){
      binit_noCI <- c(binit_noCI,0)
      names.param <- c(names.param, 'variability inter 12')
    }
    if("variability intra" %in% sharedtype_12){
      binit_noCI <- c(binit_noCI,0)
      names.param <- c(names.param, 'variability intra 12')
    }

    if(is.null(Objectlsmm$result_step2)){
      binit_noCI <- c(binit_noCI, Objectlsmm$result_step1$b)
    }
    else{
      binit_noCI <- c(binit_noCI, Objectlsmm$result_step2$b)
    }

    if(variability_inter_visit && variability_intra_visit){
      Zq1 <- generate_sobol_owen_set(S1,  nb.e.a+2)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      if(variability_inter_visit || variability_intra_visit){
        Zq1 <- generate_sobol_owen_set(S1,  nb.e.a+1)
        Zq <- apply(Zq1, 2, qnorm)
      }
      else{
        Zq1 <- generate_sobol_owen_set(S1,  nb.e.a)
        Zq <- apply(Zq1, 2, qnorm)
      }
    }

    binit_noCI <- round(binit_noCI, 5)
    estimation.noCI <- marqLevAlg(binit_noCI, fn = logR_llh_lsjm_interintraIDM, minimize = FALSE,

                                        hazard_baseline_01 = hazard_baseline_01, sharedtype_01 = sharedtype_01,
                                        hazard_baseline_02 = hazard_baseline_02, sharedtype_02 = sharedtype_02,
                                        hazard_baseline_12 = hazard_baseline_12, sharedtype_12 = sharedtype_12,
                                        ord.splines = nb.knots.splines + 2, nb.beta = nb.beta, Zq = Zq, nb_pointsGK = nb_pointsGK,
                                        nb.e.a = nb.e.a, S = S1, wk = gaussKronrod()$wk, rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk)), sk_GK = gaussKronrod()$sk, nb.alpha = nb.alpha,
                                        variability_inter_visit = variability_inter_visit, variability_intra_visit = variability_intra_visit,
                                        correlated_re = correlated_re, Case1 = NULL, Case1bis = Case1bis, Case2 = Case2, Case3 = NULL,
                                        nbCase1 = nbCase1_noCI, nbCase1bis = nbCase1bis_noCI, nbCase2 = nbCase2_noCI, nbCase3 = nbCase3_noCI, left_trunc = left_trunc,
                                        knots.hazard_baseline.splines_01 = knots_01,
                                        knots.hazard_baseline.splines_02 = knots_02,
                                        knots.hazard_baseline.splines_12 = knots_12,
                                        index_beta_slope = index_beta_slope, index_b_slope = index_b_slope,
                                        nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                                        file = file, blinding = FALSE, epsa = epsa, epsb = epsb, epsd = 1.9)


  }


  message("First estimation")
  Case1 <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0)
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0) ; B_T_12 = as.matrix(0)
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L_R = as.matrix(0); U_GK_L_R = as.matrix(0)
  Xslope_GK_L_R = as.matrix(0); Uslope_GK_L_R = as.matrix(0);
  st_L_R = as.matrix(0); Bs_L_R_01 = as.matrix(0); Bs_L_R_02 = as.matrix(0); Bs_L_R_12 = as.matrix(0); Time_T0 = c(0)
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0); X_0_LR = as.matrix(0); U_0_LR = as.matrix(0); Xslope_0_LR = as.matrix(0);
  Uslope_0_LR = as.matrix(0); Bs_0_LR_01 <- as.matrix(0);Bs_0_LR_02 <- as.matrix(0); Bs_0_LR_12 <- as.matrix(0)
  if(length(id.Case1)>0){

    data.long.Case1 <- data.long[which(data.long$id %in% id.Case1),]
    data.id.Case1 <- data.long.Case1[!duplicated(data.long.Case1$id),]
    list.long.Case1 <- data.manag.long(formGroup,formFixed, formRandom,data.long.Case1)
    X_base <- list.long.Case1$X; U_base <- list.long.Case1$U; y.new <- list.long.Case1$y.new
    ID.visit <- data.long.Case1[all.vars(formGroupVisit)][,1]; offset <- list.long.Case1$offset

    offset_ID <- c()
    len_visit <- c(0)
    for(oo in 1:nbCase1){
      ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      len_visit <- c(len_visit,length(unique(ID.visit_i)))
      frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
      offset_ID <- rbind(offset_ID, frame_offset_ID_i)
    }
    offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))

    list.GK_T <- data.GaussKronrod(data.id.Case1, a = 0, b = data.id.Case1$Time_T, k = nb_pointsGK)
    list.GK_L_R <- data.GaussKronrod(data.id.Case1, a = data.id.Case1$Time_L, b = data.id.Case1$Time_R, k = nb_pointsGK)
    st_T <- list.GK_T$st
    st_L_R <- list.GK_L_R$st
    if(left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case1, a = 0, b = data.id.Case1$Time_T0, k = nb_pointsGK)
      st_T0 <- list.GK_T0$st
      Time_T0 <- data.id.Case1$Time_T0
    }
    if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
      list.data_T <- data.time(data.id.Case1, data.id.Case1$Time_T, formFixed, formRandom,timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),formFixed, formRandom,timeVar)
      list.data.GK_L_R <- data.time(list.GK_L_R$data.id2, c(t(st_L_R)),formFixed, formRandom,timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
      X_GK_L_R <- list.data.GK_L_R$Xtime; U_GK_L_R <- list.data.GK_L_R$Utime

      if(left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     formFixed, formRandom,timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
      list.data_T <- data.time(data.id.Case1, data.id.Case1$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), formSlopeFixed, formSlopeRandom, timeVar)
      list.data.GK_L_R <- data.time(list.GK_L_R$data.id2, c(t(st_L_R)), formSlopeFixed, formSlopeRandom, timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
      Xslope_GK_L_R <- list.data.GK_L_R$Xtime; Uslope_GK_L_R <- list.data.GK_L_R$Utime

      if(left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     formSlopeFixed, formSlopeRandom,timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }
    list.surv <- data.manag.surv(formGroup, formSurv_01, data.long.Case1)
    Z_01 <- list.surv$Z
    if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_02, data.long.Case1)
    Z_02 <- list.surv$Z
    if(hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_12, data.long.Case1)
    Z_12 <- list.surv$Z
    if(hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
    if(hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(knots_01, data.id.Case1$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
      Bs_L_R_01 <- splineDesign(knots_01, c(t(st_L_R)), ord = 4L)
      if(left_trunc){
        Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(knots_02, data.id.Case1$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(knots_02, c(t(st_T)), ord = 4L)
      Bs_L_R_02 <- splineDesign(knots_02, c(t(st_L_R)), ord = 4L)
      if(left_trunc){
        Bs_T0_02 <- splineDesign(knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_12 == "Splines"){
      Z_12 <- as.matrix(Z_12[,-1])
      B_T_12 <- splineDesign(knots_12, data.id.Case1$Time_T, ord = 4L)
      Bs_T_12 <- splineDesign(knots_12, c(t(st_T)), ord = 4L)
      Bs_L_R_12 <- splineDesign(knots_12, c(t(st_L_R)), ord = 4L)
    }
    ## Pour l'intégrale (à optmiser plus tard)
    print("go integrale Case1")
    st_0_LR <- c()
    X_0_LR <- c()
    U_0_LR <- c()
    if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
      X_0_LR <- c()
      U_0_LR <- c()
    }
    if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
      Xslope_0_LR <- c()
      Uslope_0_LR <- c()
    }

    if(hazard_baseline_01== "Splines"){
      Bs_0_LR_01 <- c();
    }
    if(hazard_baseline_02== "Splines"){
      Bs_0_LR_02 <- c();
    }
    if(hazard_baseline_12== "Splines"){
      Bs_0_LR_12 <- c();
    }
    for(id.integrale in 1:nbCase1){
      data.id.integrale <- data.id.Case1[id.integrale,]
      st_L_R_i <- st_L_R[id.integrale,]
      for(st.integrale in st_L_R_i){
        list.GK_0_stLR <- data.GaussKronrod(data.id.integrale, a = 0, b = st.integrale, k = nb_pointsGK)
        st_0_stLR_i <- list.GK_0_stLR$st
        st_0_LR <- rbind(st_0_LR, st_0_stLR_i)
        if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
          list.data.GK_0_stLR <- data.time(list.GK_0_stLR$data.id2, c(t(st_0_stLR_i)),formFixed, formRandom,timeVar)
          X_0_stLR_i <- list.data.GK_0_stLR$Xtime; U_0_stLR_i <- list.data.GK_0_stLR$Utime
          X_0_LR <- rbind(X_0_LR,X_0_stLR_i); U_0_LR <- rbind(U_0_LR,U_0_stLR_i)
        }
        if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
          list.data.GK_0_stLR <- data.time(list.GK_0_stLR$data.id2, c(t(st_0_stLR_i)),formSlopeFixed, formSlopeRandom,timeVar)
          Xslope_0_stLR_i <- list.data.GK_0_stLR$Xtime; Uslope_0_stLR_i <- list.data.GK_0_stLR$Utime
          Xslope_0_LR <- rbind(Xslope_0_LR,Xslope_0_stLR_i); Uslope_0_LR <- rbind(Uslope_0_LR,Uslope_0_stLR_i)
        }
        if(hazard_baseline_01 == "Splines"){
          Bs_0_LR_01 <- rbind(Bs_0_LR_01,splineDesign(knots_01, c(t(st_0_stLR_i)), ord = 4L))
        }
        if(hazard_baseline_02 == "Splines"){
          Bs_0_LR_02 <- rbind(Bs_0_LR_02,splineDesign(knots_02, c(t(st_0_stLR_i)), ord = 4L))
        }
        if(hazard_baseline_12 == "Splines"){
          Bs_0_LR_12 <- rbind(Bs_0_LR_12,splineDesign(knots_12, c(t(st_0_stLR_i)), ord = 4L))
        }
      }
    }
    Case1 <- list( "delta2" = data.id.Case1$delta2, "Z_01" = Z_01, "Z_02" = Z_02, "Z_12" = Z_12, "Time_T" = data.id.Case1$Time_T,
                   "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                   "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                   "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01), "B_T_02" = t(B_T_02), "B_T_12" = t(B_T_12),
                   "Bs_T_01" = Bs_T_01, "Bs_T_02" = Bs_T_02, "Bs_T_12" = Bs_T_12, "X_GK_L_R" = X_GK_L_R, "U_GK_L_R" = U_GK_L_R,
                   "Xslope_GK_L_R" = Xslope_GK_L_R, "Uslope_GK_L_R" = Uslope_GK_L_R, "Time_L_R" = data.id.Case1$Time_R-data.id.Case1$Time_L,
                   "st_L_R" = t(st_L_R), "Bs_L_R_01" = Bs_L_R_01, "Bs_L_R_02" = Bs_L_R_02, "Bs_L_R_12" = Bs_L_R_12, "Time_T0" = Time_T0,
                   "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                   "Bs_T0_01" = Bs_T0_01, "Bs_T0_02" = Bs_T0_02, "data.longCase1" = data.long.Case1, "id.Case1" = id.Case1,
                   "st_0_LR" = st_0_LR, "X_0_LR" = X_0_LR, "U_0_LR" = U_0_LR, "Xslope_0_LR" = Xslope_0_LR, "Uslope_0_LR" = Uslope_0_LR,
                   "Bs_0_LR_01" = Bs_0_LR_01, "Bs_0_LR_02" = Bs_0_LR_02, "Bs_0_LR_12" = Bs_0_LR_12, "Time_L" = data.id.Case1$Time_L,
                   "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
                   # ,"random_effects" = re_Case1, "var_random_effects" = var_re_Case1
    )
  }

  Case1bis <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); B_T_01 = as.matrix(0); B_T_02 = as.matrix(0); B_T_12 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L = as.matrix(0); U_GK_L = as.matrix(0);
  Xslope_GK_L = as.matrix(0); Uslope_GK_L = as.matrix(0);
  st_L = as.matrix(0); Bs_L_01 = as.matrix(0); Bs_L_02 = as.matrix(0); Bs_L_12 = as.matrix(0); B_L_01 = as.matrix(0); B_L_02 = as.matrix(0); B_L_12 = as.matrix(0);
  X_L = as.matrix(0); U_L = as.matrix(0); Xslope_L = as.matrix(0); Uslope_L = as.matrix(0); Time_T0 = c(0)
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)
  if(length(id.Case1bis)>0){
    data.long.Case1bis <- data.long[which(data.long$id %in% id.Case1bis),]
    data.id.Case1bis <- data.long.Case1bis[!duplicated(data.long.Case1bis$id),]
    list.long.Case1bis <- data.manag.long(formGroup,formFixed, formRandom,data.long.Case1bis)
    X_base <- list.long.Case1bis$X; U_base <- list.long.Case1bis$U; y.new <- list.long.Case1bis$y.new
    ID.visit <- data.long.Case1bis[all.vars(formGroupVisit)][,1]; offset <- list.long.Case1bis$offset

    offset_ID <- c()
    len_visit <- c(0)
    for(oo in 1:nbCase1bis){
      ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      len_visit <- c(len_visit,length(unique(ID.visit_i)))
      frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
      offset_ID <- rbind(offset_ID, frame_offset_ID_i)
    }
    offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))



    list.GK_T <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_T, k = nb_pointsGK)
    list.GK_L <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_L, k = nb_pointsGK)
    st_T <- list.GK_T$st
    st_L <- list.GK_L$st
    if(left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_T0, k = nb_pointsGK)
      st_T0 <- list.GK_T0$st
      Time_T0 <- data.id.Case1bis$Time_T0
    }
    if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
      list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, formFixed, formRandom,timeVar)
      list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L, formFixed, formRandom,timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),formFixed, formRandom,timeVar)
      list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)),formFixed, formRandom,timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_L <- list.data_L$Xtime; U_L <- list.data_L$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
      X_GK_L <- list.data.GK_L$Xtime; U_GK_L <- list.data.GK_L$Utime

      if(left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     formFixed, formRandom,timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
      list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
      list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L, formSlopeFixed, formSlopeRandom, timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), formSlopeFixed, formSlopeRandom, timeVar)
      list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)), formSlopeFixed, formSlopeRandom, timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_L <- list.data_L$Xtime; Uslope_L <- list.data_L$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
      Xslope_GK_L <- list.data.GK_L$Xtime; Uslope_GK_L <- list.data.GK_L$Utime

      if(left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     formSlopeFixed, formSlopeRandom,timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    list.surv <- data.manag.surv(formGroup, formSurv_01, data.long.Case1bis)
    Z_01 <- list.surv$Z
    if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_02, data.long.Case1bis)
    Z_02 <- list.surv$Z
    if(hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_12, data.long.Case1bis)
    Z_12 <- list.surv$Z
    if(hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
    if(hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(knots_01, data.id.Case1bis$Time_T, ord = 4L)
      B_L_01 <- splineDesign(knots_01, data.id.Case1bis$Time_L, ord = 4L)
      Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
      Bs_L_01 <- splineDesign(knots_01, c(t(st_L)), ord = 4L)
      if(left_trunc){
        Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(knots_02, data.id.Case1bis$Time_T, ord = 4L)
      B_L_02 <- splineDesign(knots_02, data.id.Case1bis$Time_L, ord = 4L)
      Bs_T_02 <- splineDesign(knots_02, c(t(st_T)), ord = 4L)
      Bs_L_02 <- splineDesign(knots_02, c(t(st_L)), ord = 4L)
      if(left_trunc){
        Bs_T0_02 <- splineDesign(knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_12 == "Splines"){
      Z_12 <- as.matrix(Z_12[,-1])
      B_T_12 <- splineDesign(knots_12, data.id.Case1bis$Time_T, ord = 4L)
      B_L_12 <- splineDesign(knots_12, data.id.Case1bis$Time_L, ord = 4L)
      Bs_T_12 <- splineDesign(knots_12, c(t(st_T)), ord = 4L)
      Bs_L_12 <- splineDesign(knots_12, c(t(st_L)), ord = 4L)
    }

    Case1bis <- list( "delta2" = data.id.Case1bis$delta2, "Z_01" = Z_01, "Z_02" = Z_02, "Z_12" = Z_12, "Time_T" = data.id.Case1bis$Time_T,
                      "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                      "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                      "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01), "B_T_02" = t(B_T_02), "B_T_12" = t(B_T_12),
                      "Bs_T_01" = Bs_T_01, "Bs_T_02" = Bs_T_02, "Bs_T_12" = Bs_T_12, "X_GK_L" = X_GK_L, "U_GK_L" = U_GK_L,
                      "Xslope_GK_L" = Xslope_GK_L, "Uslope_GK_L" = Uslope_GK_L, "Time_L" = data.id.Case1bis$Time_L,
                      "st_L" = t(st_L), "Bs_L_01" = Bs_L_01, "Bs_L_02" = Bs_L_02, "Bs_L_12" = Bs_L_12, "B_L_01" = t(B_L_01), "B_L_02" = t(B_L_02), "B_L_12" = t(B_L_12),
                      "X_L" = X_L, "U_L" = t(U_L), "Xslope_L" = Xslope_L, "Uslope_L" = t(Uslope_L), "Time_T0" = Time_T0,
                      "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                      "Bs_T0_01" = Bs_T0_01, "Bs_T0_02" = Bs_T0_02, "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
                      #,"random_effects" = re_Case1bis, "var_random_effects" = var_re_Case1bis
    )

  }

  Case2 <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)
  if(length(id.Case2)>0){
    data.long.Case2 <- data.long[which(data.long$id %in% id.Case2),]
    data.id.Case2 <- data.long.Case2[!duplicated(data.long.Case2$id),]
    list.long.Case2 <- data.manag.long(formGroup,formFixed, formRandom,data.long.Case2)
    X_base <- list.long.Case2$X; U_base <- list.long.Case2$U; y.new <- list.long.Case2$y.new
    ID.visit <- data.long.Case2[all.vars(formGroupVisit)][,1]; offset <- list.long.Case2$offset

    offset_ID <- c()
    len_visit <- c(0)
    for(oo in 1:nbCase2){
      ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      len_visit <- c(len_visit,length(unique(ID.visit_i)))
      frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
      offset_ID <- rbind(offset_ID, frame_offset_ID_i)
    }
    offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))

    list.GK_T <- data.GaussKronrod(data.id.Case2, a = 0, b = data.id.Case2$Time_T, k = nb_pointsGK)
    st_T <- list.GK_T$st
    if(left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case2, a = 0, b = data.id.Case2$Time_T0, k = nb_pointsGK)
      st_T0 <- list.GK_T0$st
      Time_T0 <- data.id.Case2$Time_T0
    }
    if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) ){
      list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, formFixed, formRandom,timeVar)
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

    if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) ){
      list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
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

    list.surv <- data.manag.surv(formGroup, formSurv_01, data.long.Case2)
    Z_01 <- list.surv$Z
    if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_02, data.long.Case2)
    Z_02 <- list.surv$Z
    if(hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    if(hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(knots_01, data.id.Case2$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
      if(left_trunc){
        Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(knots_02, data.id.Case2$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(knots_02, c(t(st_T)), ord = 4L)
      if(left_trunc){
        Bs_T0_02 <- splineDesign(knots_02, c(t(st_T0)), ord = 4L)
      }
    }


    Case2 <- list( "delta2" = data.id.Case2$delta2, "Z_01" = Z_01, "Z_02" = Z_02, "Z_12" = Z_12, "Time_T" = data.id.Case2$Time_T,
                   "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                   "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                   "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01), "B_T_02" = t(B_T_02),
                   "Bs_T_01" = Bs_T_01, "Bs_T_02" = Bs_T_02, "Time_T0" = Time_T0,
                   "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                   "Bs_T0_01" = Bs_T0_01, "Bs_T0_02" = Bs_T0_02, "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
                   #, "random_effects" = re_Case2, "var_random_effects" = var_re_Case2
    )

  }

  Case3 <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); B_T_01 = as.matrix(0); B_T_02 = as.matrix(0); B_T_12 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L_T = as.matrix(0); U_GK_L_T = as.matrix(0);
  Xslope_GK_L_T = as.matrix(0); Uslope_GK_L_T = as.matrix(0); Time_L_T = c(0);
  st_L_T = as.matrix(0); B_GK_L_T_01=as.matrix(0); B_GK_L_T_02 = as.matrix(0); B_GK_L_T_12 = as.matrix(0); Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0); X_0_LT = as.matrix(0); U_0_LT = as.matrix(0); Xslope_0_LT = as.matrix(0);
  Uslope_0_LT = as.matrix(0); Bs_0_LT_01 <- as.matrix(0);Bs_0_LT_02 <- as.matrix(0); Bs_0_LT_12 <- as.matrix(0)
  if(length(id.Case3)>0){
    data.long.Case3 <- data.long[which(data.long$id %in% id.Case3),]
    data.id.Case3 <- data.long.Case3[!duplicated(data.long.Case3$id),]
    list.long.Case3 <- data.manag.long(formGroup,formFixed, formRandom,data.long.Case3)
    X_base <- list.long.Case3$X; U_base <- list.long.Case3$U; y.new <- list.long.Case3$y.new
    ID.visit <- data.long.Case3[all.vars(formGroupVisit)][,1]; offset <- list.long.Case3$offset

    offset_ID <- c()
    len_visit <- c(0)
    for(oo in 1:nbCase3){
      ID.visit_i <- ID.visit[offset[oo]:(offset[oo+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      len_visit <- c(len_visit,length(unique(ID.visit_i)))
      frame_offset_ID_i <- cbind(offset_ID_i, rep(oo, length(offset_ID_i)))
      offset_ID <- rbind(offset_ID, frame_offset_ID_i)
    }
    offset_position <- as.vector(c(1, 1 + cumsum(tapply(offset_ID[,2], offset_ID[,2], length))))

    list.GK_T <- data.GaussKronrod(data.id.Case3, a = 0, b = data.id.Case3$Time_T, k = nb_pointsGK)
    list.GK_L_T <- data.GaussKronrod(data.id.Case3, a = data.id.Case3$Time_L, b = data.id.Case3$Time_T, k = nb_pointsGK)
    st_T <- list.GK_T$st
    st_L_T <- list.GK_L_T$st
    if(left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case3, a = 0, b = data.id.Case3$Time_T0, k = nb_pointsGK)
      st_T0 <- list.GK_T0$st
      Time_T0 <- data.id.Case3$Time_T0
    }
    if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
      list.data_T <- data.time(data.id.Case3, data.id.Case3$Time_T, formFixed, formRandom,timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),formFixed, formRandom,timeVar)
      list.data.GK_L_T <- data.time(list.GK_L_T$data.id2, c(t(st_L_T)),formFixed, formRandom,timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
      X_GK_L_T <- list.data.GK_L_T$Xtime; U_GK_L_T <- list.data.GK_L_T$Utime

      if(left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     formFixed, formRandom,timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
      list.data_T <- data.time(data.id.Case3, data.id.Case3$Time_T, formSlopeFixed, formSlopeRandom, timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), formSlopeFixed, formSlopeRandom, timeVar)
      list.data.GK_L_T <- data.time(list.GK_L_T$data.id2, c(t(st_L_T)), formSlopeFixed, formSlopeRandom, timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
      Xslope_GK_L_T <- list.data.GK_L_T$Xtime; Uslope_GK_L_T <- list.data.GK_L_T$Utime

      if(left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     formSlopeFixed, formSlopeRandom,timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    list.surv <- data.manag.surv(formGroup, formSurv_01, data.long.Case3)
    Z_01 <- list.surv$Z
    if(hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_02, data.long.Case3)
    Z_02 <- list.surv$Z
    if(hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    list.surv <- data.manag.surv(formGroup, formSurv_12, data.long.Case3)
    Z_12 <- list.surv$Z
    if(hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
    if(hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(knots_01, data.id.Case3$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(knots_01, c(t(st_T)), ord = 4L)
      B_GK_L_T_01 <- splineDesign(knots_01, c(t(st_L_T)), ord = 4L)
      if(left_trunc){
        Bs_T0_01 <- splineDesign(knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(knots_02, data.id.Case3$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(knots_02, c(t(st_T)), ord = 4L)
      B_GK_L_T_02 <- splineDesign(knots_02, c(t(st_L_T)), ord = 4L)
      if(left_trunc){
        Bs_T0_02 <- splineDesign(knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    if(hazard_baseline_12 == "Splines"){
      Z_12 <- as.matrix(Z_12[,-1])
      B_T_12 <- splineDesign(knots_12, data.id.Case3$Time_T, ord = 4L)
      Bs_T_12 <- splineDesign(knots_12, c(t(st_T)), ord = 4L)
      B_GK_L_T_12 <- splineDesign(knots_12, c(t(st_L_T)), ord = 4L)
    }

    ## Pour l'intégrale (à optmiser plus tard)
    print("go integrale Case3")
    st_0_LT <- c()
    if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
      X_0_LT <- c()
      U_0_LT <- c()
    }
    if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
      Xslope_0_LT <- c()
      Uslope_0_LT <- c()
    }

    if(hazard_baseline_01 == "Splines"){
      Bs_0_LT_01 <- c();
    }
    if(hazard_baseline_02== "Splines"){
      Bs_0_LT_02 <- c();
    }
    if(hazard_baseline_12== "Splines"){
      Bs_0_LT_12 <- c();
    }

    for(id.integrale in 1:nbCase3){
      # print(id.integrale)
      data.id.integrale <- data.id.Case3[id.integrale,]
      st_L_T_i <- st_L_T[id.integrale,]
      for(st.integrale in st_L_T_i){
        list.GK_0_stLT <- data.GaussKronrod(data.id.integrale, a = 0, b = st.integrale, k = nb_pointsGK)
        st_0_stLT_i <- list.GK_0_stLT$st
        st_0_LT <- rbind(st_0_LT, st_0_stLT_i)
        if(("value" %in% sharedtype_01) || ("value" %in% sharedtype_02) || ("value" %in% sharedtype_12)){
          list.data.GK_0_stLT <- data.time(list.GK_0_stLT$data.id2, c(t(st_0_stLT_i)),formFixed, formRandom,timeVar)
          X_0_stLT_i <- list.data.GK_0_stLT$Xtime; U_0_stLT_i <- list.data.GK_0_stLT$Utime
          X_0_LT <- rbind(X_0_LT,X_0_stLT_i); U_0_LT <- rbind(U_0_LT,U_0_stLT_i)
        }
        if(("slope" %in% sharedtype_01) || ("slope" %in% sharedtype_02) || ("slope" %in% sharedtype_12)){
          list.data.GK_0_stLT <- data.time(list.GK_0_stLT$data.id2, c(t(st_0_stLT_i)),formSlopeFixed, formSlopeRandom,timeVar)
          Xslope_0_stLT_i <- list.data.GK_0_stLT$Xtime; Uslope_0_stLT_i <- list.data.GK_0_stLT$Utime
          Xslope_0_LT <- rbind(Xslope_0_LT,Xslope_0_stLT_i); Uslope_0_LT <- rbind(Uslope_0_LT,Uslope_0_stLT_i)
        }
        if(hazard_baseline_01 == "Splines"){
          Bs_0_LT_01 <- rbind(Bs_0_LT_01,splineDesign(knots_01, c(t(st_0_stLT_i)), ord = 4L))
        }
        if(hazard_baseline_02 == "Splines"){
          Bs_0_LT_02 <- rbind(Bs_0_LT_02,splineDesign(knots_02, c(t(st_0_stLT_i)), ord = 4L))
        }
        if(hazard_baseline_12 == "Splines"){
          Bs_0_LT_12 <- rbind(Bs_0_LT_12,splineDesign(knots_12, c(t(st_0_stLT_i)), ord = 4L))
        }
      }
    }


    Case3 <- list( "delta2" = data.id.Case3$delta2, "Z_01" = Z_01, "Z_02" = Z_02, "Z_12" = Z_12, "Time_T" = data.id.Case3$Time_T,
                   "st_T" = t(st_T), "X_GK_T" = X_GK_T, "U_GK_T" = U_GK_T, "Xslope_GK_T" = Xslope_GK_T, "Uslope_GK_T" = Uslope_GK_T,
                   "X_T" = X_T, "U_T" = t(U_T), "Xslope_T" = Xslope_T, "Uslope_T" = t(Uslope_T), "X_base" = X_base, "U_base" = U_base,
                   "y.new" = y.new, "ID.visit" = ID.visit, "offset" = offset, "B_T_01" = t(B_T_01), "B_T_02" = t(B_T_02), "B_T_12" = t(B_T_12),
                   "Bs_T_01" = Bs_T_01, "Bs_T_02" = Bs_T_02, "Bs_T_12" = Bs_T_12, "X_GK_L_T" = X_GK_L_T, "U_GK_L_T" = U_GK_L_T,
                   "Xslope_GK_L_T" = Xslope_GK_L_T, "Uslope_GK_L_T" = Uslope_GK_L_T, "Time_L_T" = data.id.Case3$Time_T-data.id.Case3$Time_L,
                   "st_L_T" = t(st_L_T), "Bs_L_T_01" = B_GK_L_T_01, "Bs_L_T_02" = B_GK_L_T_02, "Bs_L_T_12" = B_GK_L_T_12, "Time_T0" = Time_T0,
                   "st_T0" = t(st_T0), "X_GK_T0" = X_GK_T0, "U_GK_T0" = U_GK_T0, "Xslope_GK_T0" = Xslope_GK_T0, "Uslope_GK_T0" = Uslope_GK_T0,
                   "Bs_T0_01" = Bs_T0_01, "Bs_T0_02" = Bs_T0_02,
                   "st_0_LT" = st_0_LT, "X_0_LT" = X_0_LT, "U_0_LT" = U_0_LT, "Xslope_0_LT" = Xslope_0_LT, "Uslope_0_LT" = Uslope_0_LT,
                   "Bs_0_LT_01" = Bs_0_LT_01, "Bs_0_LT_02" = Bs_0_LT_02, "Bs_0_LT_12" = Bs_0_LT_12, "Time_L" = data.id.Case3$Time_L,
                   "offset_ID" = offset_ID, "len_visit" = len_visit, "offset_position" = offset_position
    )



  }

  names.param <- c()
  name_ZO1 <- colnames(Z_01)
  name_ZO2 <- colnames(Z_02)
  name_Z12 <- colnames(Z_12)
  # 01
  if(hazard_baseline_01 == "Weibull"){
    names.param <- c(names.param, 'shape_01')
  }
  else{
    if(hazard_baseline_01 == "Gompertz"){
      names.param <- c(names.param, 'gompertz.1_01', 'gompertz.2_01')
    }
    else{
      if(hazard_baseline_01 == "Splines"){
        for(i in 1:length(opt_splines_01$par)){
          names.param <- c(names.param, paste("splines01", i, sep = "_"))
        }
      }

    }
  }
  if(!is.null(alpha_01)){
    names.param <- c(names.param, paste(name_ZO1,"01",sep = "_"))
  }
  if("random effects" %in% sharedtype_01){
    names.param <- c(names.param, paste("re",colnames(U_base),"01",sep = "_"))
  }
  if("value" %in% sharedtype_01){
    names.param <- c(names.param, 'value 01')
  }
  if("slope" %in% sharedtype_01){
    names.param <- c(names.param, 'slope 01')
  }
  if("variability inter" %in% sharedtype_01){
    names.param <- c(names.param, 'variability inter 01')
  }
  if("variability intra" %in% sharedtype_01){
    names.param <- c(names.param, 'variability intra 01')
  }
  # 02
  if(hazard_baseline_02 == "Weibull"){
    names.param <- c(names.param, 'shape_02')
  }
  else{
    if(hazard_baseline_02 == "Gompertz"){
      names.param <- c(names.param, 'gompertz.1_02', 'gompertz.2_02')
    }
    else{
      if(hazard_baseline_02 == "Splines"){
        for(i in 1:length(opt_splines_02$par)){
          names.param <- c(names.param, paste("splines02", i, sep = "_"))
        }
      }

    }
  }
  if(!is.null(alpha_02)){
    names.param <- c(names.param, paste(name_ZO2,"02",sep = "_"))
  }
  if("random effects" %in% sharedtype_02){
    names.param <- c(names.param, paste("re",colnames(U_base),"02",sep = "_"))
  }
  if("value" %in% sharedtype_02){
    names.param <- c(names.param, 'value 02')
  }
  if("slope" %in% sharedtype_02){
    names.param <- c(names.param, 'slope 02')
  }
  if("variability inter" %in% sharedtype_02){
    names.param <- c(names.param, 'variability inter 02')
  }
  if("variability intra" %in% sharedtype_02){
    names.param <- c(names.param, 'variability intra 02')
  }
  # 12
  if(hazard_baseline_12 == "Weibull"){
    names.param <- c(names.param, 'shape_12')
  }
  else{
    if(hazard_baseline_12 == "Gompertz"){
      names.param <- c(names.param, 'gompertz.1_12', 'gompertz.2_12')
    }
    else{
      if(hazard_baseline_12 == "Splines"){
        for(i in 1:length(opt_splines_12$par)){
          names.param <- c(names.param, paste("splines12", i, sep = "_"))
        }
      }

    }
  }
  if(!is.null(alpha_12)){
    names.param <- c(names.param, paste(name_Z12,"12",sep = "_"))
  }
  if("random effects" %in% sharedtype_12){
    names.param <- c(names.param, paste("re",colnames(U_base),"12",sep = "_"))
  }
  if("value" %in% sharedtype_12){
    names.param <- c(names.param, 'value 12')
  }
  if("slope" %in% sharedtype_12){
    names.param <- c(names.param, 'slope 12')
  }
  if("variability inter" %in% sharedtype_12){
    names.param <- c(names.param, 'variability inter 12')
  }
  if("variability intra" %in% sharedtype_12){
    names.param <- c(names.param, 'variability intra 12')
  }


  if(is.null(binit)){
    if(is.null(estimation.noCI)){
      binit <- binit_CI
    }
    else{
      binit <- estimation.noCI$b
    }
  }
  else{
    if(length(binit) != length(binit_CI)){
      stop("binit has not the correct number of arguments")
    }
  }

  nb.beta <-  ncol(X_base)
  nb.alpha <- c(ncol(Z_01), ncol(Z_02), ncol(Z_12))
  nb.e.a <- ncol(U_base)

  if(variability_inter_visit && variability_intra_visit){
    Zq1 <- generate_sobol_owen_set(S1,  nb.e.a+2)
    Zq <- apply(Zq1, 2, qnorm)
  }
  else{
    if(variability_inter_visit || variability_intra_visit){
      Zq1 <- generate_sobol_owen_set(S1,  nb.e.a+1)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      Zq1 <- generate_sobol_owen_set(S1,  nb.e.a)
      Zq <- apply(Zq1, 2, qnorm)
    }
  }
  message(paste("First estimation with ", S1, " QMC draws"))
  estimation1 <- marqLevAlg(binit, fn = logR_llh_lsjm_interintraIDM, minimize = FALSE,
                                 hazard_baseline_01 = hazard_baseline_01, sharedtype_01 = sharedtype_01,
                                 hazard_baseline_02 = hazard_baseline_02, sharedtype_02 = sharedtype_02,
                                 hazard_baseline_12 = hazard_baseline_12, sharedtype_12 = sharedtype_12,
                                 ord.splines = nb.knots.splines + 2, nb.beta = nb.beta, Zq = Zq, nb_pointsGK = nb_pointsGK,
                                 nb.e.a = nb.e.a, S = S1, wk = gaussKronrod()$wk, rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk)), sk_GK = gaussKronrod()$sk, nb.alpha = nb.alpha,
                                 variability_inter_visit = variability_inter_visit, variability_intra_visit = variability_intra_visit,
                                 correlated_re = correlated_re, Case1 = Case1, Case1bis = Case1bis, Case2 = Case2, Case3 = Case3,
                                 nbCase1 = nbCase1, nbCase1bis = nbCase1bis, nbCase2 = nbCase2, nbCase3 = nbCase3, left_trunc = left_trunc,
                                 knots.hazard_baseline.splines_01 = knots_01,
                                 knots.hazard_baseline.splines_02 = knots_02,
                                 knots.hazard_baseline.splines_12 = knots_12,
                                 index_beta_slope = index_beta_slope,index_b_slope = index_b_slope,
                                 nproc = nproc, clustertype = clustertype, maxiter = maxiter, print.info = print.info,
                                 file = file, blinding = FALSE, epsa = epsa, epsb = epsb, epsd = epsd)


  estimation2 <- NULL
  info_conv_step2 <- NULL

  if(!is.null(S2) & estimation1$istop == 1){
    message(paste("Second estimation with ", S2, " QMC draws"))
    if(variability_inter_visit && variability_intra_visit){
      Zq1 <- generate_sobol_owen_set(S2,  nb.e.a+2)
      Zq <- apply(Zq1, 2, qnorm)
    }
    else{
      if(variability_inter_visit || variability_intra_visit){
        Zq1 <- generate_sobol_owen_set(S2,  nb.e.a+1)
        Zq <- apply(Zq1, 2, qnorm)
      }
      else{
        Zq1 <- generate_sobol_owen_set(S2,  nb.e.a)
        Zq <- apply(Zq1, 2, qnorm)
      }
    }

    estimation2 <- marqLevAlg(estimation1$b, fn = logR_llh_lsjm_interintraIDM, minimize = FALSE,

                              hazard_baseline_01 = hazard_baseline_01, sharedtype_01 = sharedtype_01,
                              hazard_baseline_02 = hazard_baseline_02, sharedtype_02 = sharedtype_02,
                              hazard_baseline_12 = hazard_baseline_12, sharedtype_12 = sharedtype_12,
                              ord.splines = nb.knots.splines + 2, nb.beta = nb.beta, Zq = Zq,
                              nb.e.a = nb.e.a, S = S2, wk = gaussKronrod()$wk, rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk)), sk_GK = gaussKronrod()$sk,nb_pointsGK = nb_pointsGK,
                              nb.alpha = nb.alpha,
                              variability_inter_visit = variability_inter_visit, variability_intra_visit = variability_intra_visit,
                              correlated_re = correlated_re, Case1 = Case1, Case1bis = Case1bis, Case2 = Case2, Case3 = Case3,
                              nbCase1 = nbCase1, nbCase1bis = nbCase1bis, nbCase2 = nbCase2, nbCase3 = nbCase3, left_trunc = left_trunc,
                              knots.hazard_baseline.splines_01 = knots_01,
                              knots.hazard_baseline.splines_02 = knots_02,
                              knots.hazard_baseline.splines_12 = knots_12,
                              index_beta_slope = index_beta_slope,index_b_slope = index_b_slope,

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
                                  hazard_baseline_02 = hazard_baseline_02,
                                  hazard_baseline_12 = hazard_baseline_12,
                                  nb.knots.splines = nb.knots.splines,
                                  formSurv_01 = formSurv_01,
                                  formSurv_02 = formSurv_02,
                                  formSurv_12 = formSurv_12,
                                  nb_pointsGK = nb_pointsGK,
                                  sharedtype_01 = sharedtype_01,
                                  sharedtype_02 = sharedtype_02,
                                  sharedtype_12 = sharedtype_12,
                                  formSlopeFixed = formSlopeFixed,
                                  formSlopeRandom = formSlopeRandom,
                                  index_b_slope = index_b_slope,
                                  index_beta_slope = index_beta_slope,
                                  knots.hazard_baseline.splines_01= knots_01,
                                  knots.hazard_baseline.splines_02 = knots_02,
                                  knots.hazard_baseline.splines_12 = knots_12,
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

  class(result) <- c("lsjm_interintraIDM")
  return(result)


}


































