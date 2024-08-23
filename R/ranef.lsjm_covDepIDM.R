#' @rdname ranef
#' @export
#'

ranef.lsjm_covDepIDM <- function(object,...){

  x <- object
  param <- x$result_step2$b
  cv.Pred <- c()
  x$control$nproc <- 1


  shape_01 <- 0; shape_02 <- 0; shape_12 <- 0
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0; Gompertz.1_02 <- 0; Gompertz.2_02 <- 0; Gompertz.1_12 <- 0; Gompertz.2_12 <- 0
  alpha.current_01 <- 0; alpha.current_02 <- 0; alpha.current_12 <- 0; alpha.slope_01 <- 0; alpha.slope_02 <- 0;alpha.slope_12 <- 0
  alpha_01 <- c(0); alpha_02 <- c(0); alpha_12<- c(0);
  gamma_01 <- c(0); gamma_02 <- c(0); gamma_12 <- c(0);
  beta_slope <- c(0);  omega <- c(0); alpha.var_01 <- 0; alpha.var_02 <- 0; alpha.var_12 <- 0

  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(x$control$hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]**2
    curseur <- curseur + 1
  }

  if(x$control$hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]**2
    Gompertz.2_01 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[1]+2 + 2
  }
  ### Covariables :
  nb.alpha_01 <- x$control$nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("current value" %in% x$control$sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% x$control$sharedtype_01){
    alpha.var_01 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Risque 02
  if(x$control$hazard_baseline_02 == "Weibull"){
    shape_02 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline_02 == "Gompertz"){
    Gompertz.1_02 <- param[curseur]**2
    Gompertz.2_02 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_02 == "Splines"){
    gamma_02 <- param[(curseur):(curseur+x$control$nb.knots.splines[2]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[2]+2+ 2
  }
  ### Covariables :
  nb.alpha_02 <- x$control$nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("current value" %in% x$control$sharedtype_02){
    alpha.current_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_02){
    alpha.slope_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% x$control$sharedtype_02){
    alpha.var_02 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Risque 12
  if(x$control$hazard_baseline_12 == "Weibull"){
    shape_12 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline_12 == "Gompertz"){
    Gompertz.1_12 <- param[curseur]**2
    Gompertz.2_12 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_12 == "Splines"){
    gamma_12 <- param[(curseur):(curseur+x$control$nb.knots.splines[3]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[3]+2 + 2
  }
  ### Covariables :
  nb.alpha_12 <-  x$control$nb.alpha[3]
  if(nb.alpha_12 >=1){
    alpha_12 <- param[(curseur):(curseur+nb.alpha_12-1)]
    curseur <- curseur+nb.alpha_12
  }
  ### Association
  if("current value" %in%  x$control$sharedtype_12){
    alpha.current_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in%  x$control$sharedtype_12){
    alpha.slope_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in%  x$control$sharedtype_12){
    alpha.var_12 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
  if( "slope" %in%  x$control$sharedtype_01 || "slope" %in%  x$control$sharedtype_02 || "slope" %in%  x$control$sharedtype_12){
    beta_slope <- beta[x$control$index_beta_slope]
  }
  curseur <- curseur+x$control$Objectlsmm$control$nb.beta
  omega <- param[(curseur):(curseur+x$control$Objectlsmm$control$nb.omega-1)]
  curseur <- curseur+x$control$Objectlsmm$control$nb.omega
  if(x$control$Objectlsmm$control$correlated_re){
    C1 <- matrix(rep(0,(x$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)**2),nrow=x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma,ncol=x$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)
    C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
    Cholesky <- C1
    Cholesky <- as.matrix(Cholesky)
  }
  else{
    borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
    C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    borne3 <- borne1 + choose(n = x$control$Objectlsmm$control$nb.e.a.sigma, k = 2) + x$control$Objectlsmm$control$nb.e.a.sigma
    C3 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a.sigma)**2),nrow=x$control$Objectlsmm$control$nb.e.a.sigma,ncol=x$control$Objectlsmm$control$nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
    C2 <- matrix(0, ncol = x$control$Objectlsmm$control$nb.e.a.sigma, nrow = x$control$Objectlsmm$control$nb.e.a)

    C4 <- matrix(0, ncol = x$control$Objectlsmm$control$nb.e.a, nrow = x$control$Objectlsmm$control$nb.e.a.sigma)
    Cholesky <- rbind(cbind(C1,C2),cbind(C4,C3))
    Cholesky <- as.matrix(Cholesky)
  }

  MatCov <- Cholesky%*%t(Cholesky)

  sharedtype <- c("current value" %in% x$control$sharedtype_01, "slope" %in% x$control$sharedtype_01, "variability" %in% x$control$sharedtype_01,
                  "current value" %in% x$control$sharedtype_02, "slope" %in% x$control$sharedtype_02, "variability" %in% x$control$sharedtype_02,
                  "current value" %in% x$control$sharedtype_12, "slope" %in% x$control$sharedtype_12, "variability" %in% x$control$sharedtype_12)
  HB <- list(x$control$hazard_baseline_01, x$control$hazard_baseline_02, x$control$hazard_baseline_12, x$control$left_trunc)
  W_G <- c(shape_01, shape_02, shape_12, Gompertz.1_01, Gompertz.2_01, Gompertz.1_02, Gompertz.2_02, Gompertz.1_12, Gompertz.2_12)
  alpha_y_slope_var <- c(alpha.current_01,alpha.current_02,alpha.current_12, alpha.slope_01,alpha.slope_02,alpha.slope_12,
                         alpha.var_01,alpha.var_02,alpha.var_12)
  alpha_z <- list(alpha_01, alpha_02, alpha_12)
  gamma_z0 <- list(gamma_01, gamma_02, gamma_12)
  wk = gaussKronrod()$wk
  rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk))
  sk_GK <- gaussKronrod()$sk
  # Différencier chaque cas, pour chaque cas créer les matrices puis faire une fonction à maximiser (dedans un appel à C++ très proche des fonctions déjà existantes)
  # => il faut appeler marqlevalg individu par individu !
  #1. On associe chaque individu à un cas
  data.long <-x$control$Objectlsmm$control$data.long
  data.long <- x$control$Objectlsmm$control$data.long
  Time_T <- x$control$Time[["Time_T"]]
  data.long$Time_T <- data.long[all.vars(Time_T)][,1]
  Time_R <- x$control$Time[["Time_R"]]
  data.long$Time_R <- data.long[all.vars(Time_R)][,1]
  Time_L <- x$control$Time[["Time_L"]]
  data.long$Time_L <- data.long[all.vars(Time_L)][,1]
  if(!is.null(x$control$Time[["Time_T0"]])){
    left_trunc <- TRUE
    Time_T0 <- x$control$Time[["Time_T0"]]
    data.long$Time_T0 <- data.long[all.vars(Time_T0)][,1]
  }
  else{
    left_trunc <- FALSE
  }
  delta1 <- x$control$deltas[["delta1"]]
  data.long$delta1 <- data.long[all.vars(delta1)][,1]
  delta2 <- x$control$deltas[["delta2"]]
  data.long$delta2 <- data.long[all.vars(delta2)][,1]
  data.long$Case <- ifelse(data.long$delta1 == 1 & data.long$Time_R == data.long$Time_L, "Case1bis",
                           ifelse(data.long$delta1 == 1, "Case1",
                                  ifelse(data.long$Time_L == data.long$Time_T, "Case2", "Case3")))
  id.Case1 <- unique(data.long$id[which(data.long$Case == "Case1")])
  id.Case1bis <- unique(data.long$id[which(data.long$Case == "Case1bis")])
  id.Case2 <- unique(data.long$id[which(data.long$Case == "Case2")])
  id.Case3 <- unique(data.long$id[which(data.long$Case == "Case3")])
  nbCase1 <- length(id.Case1); nbCase1bis <- length(id.Case1bis); nbCase2 <- length(id.Case2); nbCase3 <- length(id.Case3)
  random.effects.Predictions <- matrix(NA, nrow = nbCase1+nbCase1bis+nbCase2+nbCase3, ncol = x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)

  #2. Cas 1
  message("Case1")
  Case1 <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0)
  O_GK_T = as.matrix(0); W_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  O_T = as.matrix(0); W_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0) ; B_T_12 = as.matrix(0)
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L_R = as.matrix(0); U_GK_L_R = as.matrix(0)
  Xslope_GK_L_R = as.matrix(0); Uslope_GK_L_R = as.matrix(0);
  O_GK_L_R = as.matrix(0); W_GK_L_R = as.matrix(0)
  st_L_R = as.matrix(0); Bs_L_R_01 = as.matrix(0); Bs_L_R_02 = as.matrix(0); Bs_L_R_12 = as.matrix(0); Time_T0 = c(0)
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  O_GK_T0 = as.matrix(0); W_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0); X_0_LR = as.matrix(0); U_0_LR = as.matrix(0); Xslope_0_LR = as.matrix(0);
  Uslope_0_LR = as.matrix(0);
  O_0_LR = as.matrix(0); W_0_LR = as.matrix(0);
  Bs_0_LR_01 <- as.matrix(0);Bs_0_LR_02 <- as.matrix(0); Bs_0_LR_12 <- as.matrix(0)
  if(length(id.Case1)>0){

    data.long.Case1 <- data.long[which(data.long$id %in% id.Case1),]
    data.id.Case1 <- data.long.Case1[!duplicated(data.long.Case1$id),]
    list.long.Case1 <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long.Case1)
    X_base <- list.long.Case1$X; U_base <- list.long.Case1$U; y.new <- list.long.Case1$y.new
    time.measures <- data.long.Case1[,x$control$Objectlsmm$control$timeVar]
    offset <- list.long.Case1$offset

    list.var <- data.manag.sigma(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,data.long.Case1)
    O_base <- list.var$X
    O_base <- as.matrix(O_base)
    W_base <- list.var$U
    W_base <- as.matrix(W_base)


    list.GK_T <- data.GaussKronrod(data.id.Case1, a = 0, b = data.id.Case1$Time_T, k = x$control$nb_pointsGK)
    list.GK_L_R <- data.GaussKronrod(data.id.Case1, a = data.id.Case1$Time_L, b = data.id.Case1$Time_R, k = x$control$nb_pointsGK)
    st_T <- list.GK_T$st
    st_L_R <- list.GK_L_R$st
    if(x$control$left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case1, a = 0, b = data.id.Case1$Time_T0, k = x$control$nb_pointsGK)
      st_T0 <- list.GK_T0$st
      Time_T0 <- data.id.Case1$Time_T0
    }
    if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case1, data.id.Case1$Time_T, x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_L_R <- data.time(list.GK_L_R$data.id2, c(t(st_L_R)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
      X_GK_L_R <- list.data.GK_L_R$Xtime; U_GK_L_R <- list.data.GK_L_R$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case1, data.id.Case1$Time_T, x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_L_R <- data.time(list.GK_L_R$data.id2, c(t(st_L_R)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
      Xslope_GK_L_R <- list.data.GK_L_R$Xtime; Uslope_GK_L_R <- list.data.GK_L_R$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case1, data.id.Case1$Time_T, x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_L_R <- data.time(list.GK_L_R$data.id2, c(t(st_L_R)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      O_T <- list.data_T$Xtime; W_T <- list.data_T$Utime
      O_GK_T <- list.data.GK_T$Xtime; W_GK_T <- list.data.GK_T$Utime
      O_GK_L_R <- list.data.GK_L_R$Xtime; W_GK_L_R <- list.data.GK_L_R$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
        O_GK_T0 <- list.data.GK_T0$Xtime
        W_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long.Case1)
    Z_01 <- list.surv$Z
    if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_02, data.long.Case1)
    Z_02 <- list.surv$Z
    if(x$control$hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_12, data.long.Case1)
    Z_12 <- list.surv$Z
    if(x$control$hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
    if(x$control$hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(x$control$knots_01, data.id.Case1$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(x$control$knots_01, c(t(st_T)), ord = 4L)
      Bs_L_R_01 <- splineDesign(x$control$knots_01, c(t(st_L_R)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_01 <- splineDesign(x$control$knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(x$control$knots_02, data.id.Case1$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(x$control$knots_02, c(t(st_T)), ord = 4L)
      Bs_L_R_02 <- splineDesign(x$control$knots_02, c(t(st_L_R)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_02 <- splineDesign(x$control$knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_12 == "Splines"){
      Z_12 <- as.matrix(Z_12[,-1])
      B_T_12 <- splineDesign(x$control$knots_12, data.id.Case1$Time_T, ord = 4L)
      Bs_T_12 <- splineDesign(x$control$knots_12, c(t(st_T)), ord = 4L)
      Bs_L_R_12 <- splineDesign(x$control$knots_12, c(t(st_L_R)), ord = 4L)
    }

    ## Pour l'intégrale (à optmiser plus tard)
    #print("go integrale Case1")
    st_0_LR <- c()
    X_GK_0_LR <- c()
    U_GK_0_LR <- c()
    if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
      X_GK_0_LR <- c()
      U_GK_0_LR <- c()
    }
    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
      Xslope_GK_0_LR <- c()
      Uslope_GK_0_LR <- c()
    }
    if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
      O_GK_0_LR <- c()
      W_GK_0_LR <- c()
    }

    if(x$control$hazard_baseline_01== "Splines"){
      Bs_0_LR_01 <- c();
    }
    if(x$control$hazard_baseline_02== "Splines"){
      Bs_0_LR_02 <- c();
    }
    if(x$control$hazard_baseline_12== "Splines"){
      Bs_0_LR_12 <- c();
    }
    for(id.integrale in 1:nbCase1){
      data.id.integrale <- data.id.Case1[id.integrale,]
      st_L_R_i <- st_L_R[id.integrale,]
      for(st.integrale in st_L_R_i){
        list.GK_0_stLR <- data.GaussKronrod(data.id.integrale, a = 0, b = st.integrale, k = x$control$nb_pointsGK)
        st_0_stLR_i <- list.GK_0_stLR$st
        st_0_LR <- rbind(st_0_LR, st_0_stLR_i)
        if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
          list.data.GK_0_stLR <- data.time(list.GK_0_stLR$data.id2, c(t(st_0_stLR_i)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
          X_0_stLR_i <- list.data.GK_0_stLR$Xtime; U_0_stLR_i <- list.data.GK_0_stLR$Utime
          X_GK_0_LR <- rbind(X_GK_0_LR,X_0_stLR_i); U_GK_0_LR <- rbind(U_GK_0_LR,U_0_stLR_i)
        }
        if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
          list.data.GK_0_stLR <- data.time(list.GK_0_stLR$data.id2, c(t(st_0_stLR_i)),x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
          Xslope_0_stLR_i <- list.data.GK_0_stLR$Xtime; Uslope_0_stLR_i <- list.data.GK_0_stLR$Utime
          Xslope_GK_0_LR <- rbind(Xslope_GK_0_LR,Xslope_0_stLR_i); Uslope_GK_0_LR <- rbind(Uslope_GK_0_LR,Uslope_0_stLR_i)
        }
        if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
          list.data.GK_0_stLR <- data.time(list.GK_0_stLR$data.id2, c(t(st_0_stLR_i)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
          O_0_stLR_i <- list.data.GK_0_stLR$Xtime; W_0_stLR_i <- list.data.GK_0_stLR$Utime
          O_GK_0_LR <- rbind(O_GK_0_LR,X_0_stLR_i); W_GK_0_LR <- rbind(W_GK_0_LR,U_0_stLR_i)
        }
        if(x$control$hazard_baseline_01 == "Splines"){
          Bs_0_LR_01 <- rbind(Bs_0_LR_01,splineDesign(x$control$knots_01, c(t(st_0_stLR_i)), ord = 4L))
        }
        if(x$control$hazard_baseline_02 == "Splines"){
          Bs_0_LR_02 <- rbind(Bs_0_LR_02,splineDesign(x$control$knots_02, c(t(st_0_stLR_i)), ord = 4L))
        }
        if(x$control$hazard_baseline_12 == "Splines"){
          Bs_0_LR_12 <- rbind(Bs_0_LR_12,splineDesign(x$control$knots_12, c(t(st_0_stLR_i)), ord = 4L))
        }
      }



    }

    X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
    O_GK_T0_i <- as.matrix(0); W_GK_T0_i <- as.matrix(0); O_T_i <- c(0); W_T_i <- c(0);
    X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ; X_GK_L_R_i<- as.matrix(0) ; U_GK_L_R_i<- as.matrix(0) ; X_GK_0_LR_i<- as.matrix(0); U_GK_0_LR_i<- as.matrix(0);
    Xslope_T_i <- c(0); Uslope_T_i <- c(0); Xslope_GK_0_LR_i <- as.matrix(0); Uslope_GK_0_LR_i <- as.matrix(0);
    Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0); Xslope_GK_L_R_i <- as.matrix(0); Uslope_GK_L_R_i <- as.matrix(0)
    Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0)
    O_GK_T_i <- as.matrix(0); W_GK_T_i<- as.matrix(0) ; O_GK_L_R_i<- as.matrix(0) ; W_GK_L_R_i<- as.matrix(0) ; O_GK_0_LR_i<- as.matrix(0); W_GK_0_LR_i<- as.matrix(0);
    B_T_i_12 <- c(0);
    Bs_T_i_12 <- as.matrix(0);
    Bs_0_LR_i_01<- as.matrix(0);  Bs_0_LR_i_02<- as.matrix(0);  Bs_0_LR_i_12<- as.matrix(0);
    Bs_L_R_i_01<- as.matrix(0);
    Bs_T0_i_01<- as.matrix(0);  Bs_T0_i_02<- as.matrix(0);
    st_T_i <- c(0); st_T0_i <- c(0)

    binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a + x$control$Objectlsmm$control$nb.e.a.sigma)
    for(id.Case1_boucle in 1:nbCase1){
      delta2_i <- data.id.Case1$delta2[id.Case1_boucle]
      Z_01_i <- Z_01[id.Case1_boucle,]
      Z_02_i <- Z_02[id.Case1_boucle,]
      Z_12_i <- Z_12[id.Case1_boucle,]
      if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
        X_T_i <- X_T[id.Case1_boucle,]
        U_T_i <- U_T[id.Case1_boucle,]
        X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        X_GK_L_R_i <- as.matrix(X_GK_L_R[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        U_GK_L_R_i <- as.matrix(U_GK_L_R[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        X_GK_0_LR_i <- as.matrix(X_GK_0_LR[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),])
        U_GK_0_LR_i <- as.matrix(U_GK_0_LR[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),])
        if(x$control$left_trunc){
          X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
          U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        }
      }
      if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
        Xslope_T_i <- Xslope_T[id.Case1_boucle,]
        Uslope_T_i <- Uslope_T[id.Case1_boucle,]
        Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        Xslope_GK_L_R_i <- as.matrix(Xslope_GK_L_R[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        Uslope_GK_L_R_i <- as.matrix(Uslope_GK_L_R[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        Xslope_GK_0_LR_i <- as.matrix(Xslope_GK_0_LR[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),])
        Uslope_GK_0_LR_i <- as.matrix(Uslope_GK_0_LR[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),])
        if(x$control$left_trunc){
          Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
          Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        }
      }
      if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
        O_T_i <- O_T[id.Case1_boucle,]
        W_T_i <- W_T[id.Case1_boucle,]
        O_GK_T_i <- as.matrix(O_GK_T[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        W_GK_T_i <- as.matrix(W_GK_T[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        O_GK_L_R_i <- as.matrix(O_GK_L_R[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        W_GK_L_R_i <- as.matrix(W_GK_L_R[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        O_GK_0_LR_i <- as.matrix(O_GK_0_LR[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),])
        W_GK_0_LR_i <- as.matrix(W_GK_0_LR[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),])
        if(x$control$left_trunc){
          O_GK_T0_i <- as.matrix(O_GK_T0[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
          W_GK_T0_i <- as.matrix(W_GK_T0[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),])
        }
      }
      Time_T_i <- data.id.Case1$Time_T[id.Case1_boucle]
      Time_L_R_i <- data.id.Case1$Time_R[id.Case1_boucle]-data.id.Case1$Time_L[id.Case1_boucle]
      if(x$control$left_trunc){
        Time_T0_i <- data.id.Case1$Time_T0[id.Case1_boucle]
        st_T0_i <- st_T0[id.Case1_boucle,]
      }
      st_T_i <- st_T[id.Case1_boucle,]
      st_0_LR_i <- st_0_LR[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),]
      st_L_R_i <- st_L_R[id.Case1_boucle,]
      ck <- ((sk_GK+1)/4)*Time_L_R_i + data.id.Case1$Time_L[id.Case1_boucle]/2
      if(x$control$hazard_baseline_01 == "Splines"){
        Bs_L_R_i_01 <- Bs_L_R_01[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),]
        Bs_0_LR_i_01 <- Bs_0_LR_01[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),]
        }
      }
      if(x$control$hazard_baseline_02 == "Splines"){
        Bs_0_LR_i_02 <- Bs_0_LR_02[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),]
        }
      }
      if(x$control$hazard_baseline_02 == "Splines"){
        B_T_i_12 <- B_T_12[id.Case1_boucle,]
        Bs_T_i_12 <- Bs_T_12[(x$control$nb_pointsGK*(id.Case1_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1_boucle)),]
        Bs_0_LR_i_12 <- Bs_0_LR_12[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1_boucle),]

      }

      X_base_i <- X_base[offset[id.Case1_boucle]:(offset[id.Case1_boucle+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[id.Case1_boucle+1]-offset[id.Case1_boucle])
      X_base_i <- unique(X_base_i)
      U_i <- U_base[offset[id.Case1_boucle]:(offset[id.Case1_boucle+1]-1),]
      U_i <- matrix(U_i, nrow = offset[id.Case1_boucle+1]-offset[id.Case1_boucle])
      U_base_i <- unique(U_i)
      y_i <- y.new[offset[id.Case1_boucle]:(offset[id.Case1_boucle+1]-1)]

      W_base_i <- W_base[offset[id.Case1_boucle]:(offset[id.Case1_boucle+1]-1),]
      W_base_i <- matrix(W_base_i, nrow = offset[id.Case1_boucle+1]-offset[id.Case1_boucle])
      W_base_i <- unique(W_base_i)
      O_base_i <- O_base[offset[id.Case1_boucle]:(offset[id.Case1_boucle+1]-1),]
      O_base_i <- matrix(O_base_i, nrow = offset[id.Case1_boucle+1]-offset[id.Case1_boucle])
      O_base_i <- unique(O_base_i)

      random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase1, minimize = FALSE, nb.e.a = x$control$Objectlsmm$control$nb.e.a,
                                     nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma, Sigma.re = MatCov,
                                     sharedtype = sharedtype, HB = HB, W_G = W_G , nb_pointsGK = x$control$nb_pointsGK,
                                     alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,
                                     wk = wk, rep_wk = rep_wk,
                                     delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i, Z_12_i=Z_12_i,
                                     X_T_i=X_T_i,  U_T_i=U_T_i,Xslope_T_i =Xslope_T_i ,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                     X_GK_T_i = X_GK_T_i,  U_GK_T_i = U_GK_T_i,  Xslope_GK_T_i= Xslope_GK_T_i, Uslope_GK_T_i = Uslope_GK_T_i, O_GK_T_i = O_GK_T_i,  W_GK_T_i = W_GK_T_i,
                                     X_GK_L_R_i = X_GK_L_R_i,  U_GK_L_R_i = U_GK_L_R_i,  Xslope_GK_L_R_i = Xslope_GK_L_R_i,  Uslope_GK_L_R_i = Uslope_GK_L_R_i, O_GK_L_R_i = O_GK_L_R_i,  W_GK_L_R_i = W_GK_L_R_i,
                                     X_GK_0_LR_i = X_GK_0_LR_i,  U_GK_0_LR_i = U_GK_0_LR_i,  Xslope_GK_0_LR_i = Xslope_GK_0_LR_i,  Uslope_GK_0_LR_i = Uslope_GK_0_LR_i, O_GK_0_LR_i = O_GK_0_LR_i,  W_GK_0_LR_i = W_GK_0_LR_i,
                                     X_GK_T0_i = X_GK_T0_i,  U_GK_T0_i = U_GK_T0_i,  Xslope_GK_T0_i = Xslope_GK_T0_i,  Uslope_GK_T0_i = Uslope_GK_T0_i, O_GK_T0_i = O_GK_T0_i,  W_GK_T0_i = W_GK_T0_i,
                                     Time_T_i = Time_T_i,  Time_L_R_i = Time_L_R_i,  Time_T0_i = Time_T0_i,
                                     st_T_i = st_T_i,  st_0_LR_i = st_0_LR_i,  st_L_R_i = st_L_R_i,  st_T0_i = st_T0_i, ck = ck,
                                     B_T_i_12 = B_T_i_12,
                                     Bs_T_i_12 = Bs_T_i_12,
                                     Bs_0_LR_i_01 = Bs_0_LR_i_01,  Bs_0_LR_i_02 = Bs_0_LR_i_02,  Bs_0_LR_i_12 = Bs_0_LR_i_12,
                                     Bs_L_R_i_01 = Bs_L_R_i_01,
                                     Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02 = Bs_T0_i_02, left_trunc= x$control$left_trunc,
                                     X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i, O_base_i = O_base_i, W_base_i = W_base_i,
                                     index_b_slope = index_b_slope,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

      while(random.effects_i$istop !=1){
        binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
        random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase1, minimize = FALSE, nb.e.a = x$control$Objectlsmm$control$nb.e.a,
                                       nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma, Sigma.re = MatCov,
                                       sharedtype = sharedtype, HB = HB, W_G = W_G , nb_pointsGK = x$control$nb_pointsGK,
                                       alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,omega = omega,
                                       wk = wk, rep_wk = rep_wk,
                                       delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i, Z_12_i=Z_12_i,
                                       X_T_i=X_T_i,  U_T_i=U_T_i,Xslope_T_i =Xslope_T_i ,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                       X_GK_T_i = X_GK_T_i,  U_GK_T_i = U_GK_T_i,  Xslope_GK_T_i= Xslope_GK_T_i, Uslope_GK_T_i = Uslope_GK_T_i, O_GK_T_i = O_GK_T_i,  W_GK_T_i = W_GK_T_i,
                                       X_GK_L_R_i = X_GK_L_R_i,  U_GK_L_R_i = U_GK_L_R_i,  Xslope_GK_L_R_i = Xslope_GK_L_R_i,  Uslope_GK_L_R_i = Uslope_GK_L_R_i, O_GK_L_R_i = O_GK_L_R_i,  W_GK_L_R_i = W_GK_L_R_i,
                                       X_GK_0_LR_i = X_GK_0_LR_i,  U_GK_0_LR_i = U_GK_0_LR_i,  Xslope_GK_0_LR_i = Xslope_GK_0_LR_i,  Uslope_GK_0_LR_i = Uslope_GK_0_LR_i, O_GK_0_LR_i = O_GK_0_LR_i,  W_GK_0_LR_i = W_GK_0_LR_i,
                                       X_GK_T0_i = X_GK_T0_i,  U_GK_T0_i = U_GK_T0_i,  Xslope_GK_T0_i = Xslope_GK_T0_i,  Uslope_GK_T0_i = Uslope_GK_T0_i, O_GK_T0_i = O_GK_T0_i,  W_GK_T0_i = W_GK_T0_i,
                                       Time_T_i = Time_T_i,  Time_L_R_i = Time_L_R_i,  Time_T0_i = Time_T0_i,
                                       st_T_i = st_T_i,  st_0_LR_i = st_0_LR_i,  st_L_R_i = st_L_R_i,  st_T0_i = st_T0_i, ck = ck,
                                       B_T_i_12 = B_T_i_12,
                                       Bs_T_i_12 = Bs_T_i_12,
                                       Bs_0_LR_i_01 = Bs_0_LR_i_01,  Bs_0_LR_i_02 = Bs_0_LR_i_02,  Bs_0_LR_i_12 = Bs_0_LR_i_12,
                                       Bs_L_R_i_01 = Bs_L_R_i_01,
                                       Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02 = Bs_T0_i_02, left_trunc= x$control$left_trunc,
                                       X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i, O_base_i = O_base_i, W_base_i = W_base_i,
                                       index_b_slope = index_b_slope,
                                       nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                       file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
        binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a + x$control$Objectlsmm$control$nb.e.a.sigma)
      }

      random.effects.Predictions[id.Case1_boucle,] <- c(data.id.Case1$id[id.Case1_boucle] ,random.effects_i$b)
      CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
      time.measures_i <- time.measures[offset[id.Case1_boucle]:(offset[id.Case1_boucle+1]-1)]
      time.measures_i <- unique(time.measures_i)
      Varia <- exp(O_base_i%*%omega + W_base_i%*%random.effects_i$b[(x$control$Objectlsmm$control$nb.e.a+1):(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)])
      cv.Pred <- rbind(cv.Pred, cbind(rep(data.id.Case1$id[id.Case1_boucle], length(CV)),
                                      time.measures_i, CV,  Varia))
    }

  }

  #3. Cas 1bis
  Case1bis <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0); O_GK_T = as.matrix(0); W_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); O_T = as.matrix(0); W_T = as.matrix(0); B_T_01 = as.matrix(0); B_T_02 = as.matrix(0); B_T_12 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L = as.matrix(0); U_GK_L = as.matrix(0);
  Xslope_GK_L = as.matrix(0); Uslope_GK_L = as.matrix(0); O_GK_L = as.matrix(0); W_GK_L = as.matrix(0);
  st_L = as.matrix(0); Bs_L_01 = as.matrix(0); Bs_L_02 = as.matrix(0); Bs_L_12 = as.matrix(0); B_L_01 = as.matrix(0); B_L_02 = as.matrix(0); B_L_12 = as.matrix(0);
  X_L = as.matrix(0); U_L = as.matrix(0); Xslope_L = as.matrix(0); Uslope_L = as.matrix(0); O_L = as.matrix(0); W_L = as.matrix(0); Time_T0 = c(0)
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0); O_GK_T0 = as.matrix(0); W_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)
  if(length(id.Case1bis)>0){
    data.long.Case1bis <- data.long[which(data.long$id %in% id.Case1bis),]
    data.id.Case1bis <- data.long.Case1bis[!duplicated(data.long.Case1bis$id),]
    list.long.Case1bis <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long.Case1bis)
    X_base <- list.long.Case1bis$X; U_base <- list.long.Case1bis$U; y.new <- list.long.Case1bis$y.new
    offset <- list.long.Case1bis$offset

    list.var <- data.manag.sigma(formGroup,formFixedVar, formRandomVar,data.long.Case1bis)
    O_base <- list.var$X
    O_base <- as.matrix(O_base)
    W_base <- list.var$U
    W_base <- as.matrix(W_base)

    list.GK_T <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_T, k = x$control$nb_pointsGK)
    list.GK_L <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_L, k = x$control$nb_pointsGK)
    st_T <- list.GK_T$st
    st_L <- list.GK_L$st
    if(x$control$left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case1bis, a = 0, b = data.id.Case1bis$Time_T0, k = x$control$nb_pointsGK)
      st_T0 <- list.GK_T0$st
      Time_T0 <- data.id.Case1bis$Time_T0
    }
    if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L, x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_L <- list.data_L$Xtime; U_L <- list.data_L$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
      X_GK_L <- list.data.GK_L$Xtime; U_GK_L <- list.data.GK_L$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L, x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_L <- list.data_L$Xtime; Uslope_L <- list.data_L$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
      Xslope_GK_L <- list.data.GK_L$Xtime; Uslope_GK_L <- list.data.GK_L$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case1bis, data.id.Case1bis$Time_T, x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data_L <- data.time(data.id.Case1bis, data.id.Case1bis$Time_L, x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_L <- data.time(list.GK_L$data.id2, c(t(st_L)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      O_T <- list.data_T$Xtime; W_T <- list.data_T$Utime
      O_L <- list.data_L$Xtime; W_L <- list.data_L$Utime
      O_GK_T <- list.data.GK_T$Xtime; W_GK_T <- list.data.GK_T$Utime
      O_GK_L <- list.data.GK_L$Xtime; W_GK_L <- list.data.GK_L$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
        O_GK_T0 <- list.data.GK_T0$Xtime
        W_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long.Case1bis)
    Z_01 <- list.surv$Z
    if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_02, data.long.Case1bis)
    Z_02 <- list.surv$Z
    if(x$control$hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_12, data.long.Case1bis)
    Z_12 <- list.surv$Z
    if(x$control$hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
    if(x$control$hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(x$control$knots_01, data.id.Case1bis$Time_T, ord = 4L)
      B_L_01 <- splineDesign(x$control$knots_01, data.id.Case1bis$Time_L, ord = 4L)
      Bs_T_01 <- splineDesign(x$control$knots_01, c(t(st_T)), ord = 4L)
      Bs_L_01 <- splineDesign(x$control$knots_01, c(t(st_L)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_01 <- splineDesign(x$control$knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(x$control$knots_02, data.id.Case1bis$Time_T, ord = 4L)
      B_L_02 <- splineDesign(x$control$knots_02, data.id.Case1bis$Time_L, ord = 4L)
      Bs_T_02 <- splineDesign(x$control$knots_02, c(t(st_T)), ord = 4L)
      Bs_L_02 <- splineDesign(x$control$knots_02, c(t(st_L)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_02 <- splineDesign(x$control$knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_12 == "Splines"){
      Z_12 <- as.matrix(Z_12[,-1])
      B_T_12 <- splineDesign(x$control$knots_12, data.id.Case1bis$Time_T, ord = 4L)
      B_L_12 <- splineDesign(x$control$knots_12, data.id.Case1bis$Time_L, ord = 4L)
      Bs_T_12 <- splineDesign(x$control$knots_12, c(t(st_T)), ord = 4L)
      Bs_L_12 <- splineDesign(x$control$knots_12, c(t(st_L)), ord = 4L)
    }



    X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);  X_L_i <- c(0); U_L_i <- c(0);
    X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ; X_GK_L_i <- as.matrix(0); U_GK_L_i<- as.matrix(0) ;
    Xslope_T_i <- c(0); Uslope_T_i <- c(0);  Xslope_L_i <- c(0); Uslope_L_i <- c(0);
    Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0); Xslope_GK_L_i<- as.matrix(0) ; Uslope_GK_L_i<- as.matrix(0);
    Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
    O_GK_T0_i <- as.matrix(0); W_GK_T0_i <- as.matrix(0); O_T_i <- c(0); W_T_i <- c(0);  O_L_i <- c(0); W_L_i <- c(0);
    O_GK_T_i <- as.matrix(0); W_GK_T_i<- as.matrix(0) ; O_GK_L_i <- as.matrix(0); W_GK_L_i<- as.matrix(0) ;
    B_T_i_12 <- c(0) ; Bs_T_i_12 <- as.matrix(0) ;
    B_L_i_01 <- c(0);
    Bs_L_i_01 <- as.matrix(0) ;   Bs_L_i_02 <- as.matrix(0) ;   Bs_L_i_12 <- as.matrix(0) ;
    Bs_T0_i_01<-as.matrix(0) ;   Bs_T0_i_02 <- as.matrix(0) ;

    time.measures <- data.id.Case1bis[,x$control$Objectlsmm$control$timeVar]

    for(id.Case1bis_boucle in 1:nbCase1bis){
      delta2_i <- data.id.Case1bis$delta2[id.Case1bis_boucle]
      Z_01_i <- Z_01[id.Case1bis_boucle,]
      Z_02_i <- Z_02[id.Case1bis_boucle,]
      Z_12_i <- Z_12[id.Case1bis_boucle,]
      if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
        X_T_i <- X_T[id.Case1bis_boucle,]
        U_T_i <- U_T[id.Case1bis_boucle,]
        X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])

        X_L_i <- X_L[id.Case1bis_boucle,]
        U_L_i <- U_L[id.Case1bis_boucle,]
        X_GK_L_i <- as.matrix(X_GK_L[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        U_GK_L_i <- as.matrix(U_GK_L[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])

        if(x$control$left_trunc){
          X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
          U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        }
      }
      if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
        Xslope_T_i <- Xslope_T[id.Case1bis_boucle,]
        Uslope_T_i <- Uslope_T[id.Case1bis_boucle,]
        Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])

        Xslope_L_i <- Xslope_L[id.Case1bis_boucle,]
        Uslope_L_i <- Uslope_L[id.Case1bis_boucle,]
        Xslope_GK_L_i <- as.matrix(Xslope_GK_L[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        Uslope_GK_L_i <- as.matrix(Uslope_GK_L[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])

        if(x$control$left_trunc){
          Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
          Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        }
      }
      if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
        O_T_i <- O_T[id.Case1bis_boucle,]
        W_T_i <- W_T[id.Case1bis_boucle,]
        O_GK_T_i <- as.matrix(O_GK_T[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        W_GK_T_i <- as.matrix(W_GK_T[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])

        O_L_i <- O_L[id.Case1bis_boucle,]
        W_L_i <- W_L[id.Case1bis_boucle,]
        O_GK_L_i <- as.matrix(O_GK_L[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        W_GK_L_i <- as.matrix(W_GK_L[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])

        if(x$control$left_trunc){
          O_GK_T0_i <- as.matrix(O_GK_T0[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
          W_GK_T0_i <- as.matrix(W_GK_T0[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),])
        }
      }
      Time_T_i <- data.id.Case1bis$Time_T[id.Case1bis_boucle]
      Time_L_i <- data.id.Case1bis$Time_L[id.Case1bis_boucle]
      if(x$control$left_trunc){
        Time_T0_i <- data.id.Case1bis$Time_T0[id.Case1bis_boucle]
        st_T0_i <- st_T0[id.Case1bis_boucle,]
      }
      st_T_i <- st_T[id.Case1bis_boucle,]
      st_L_i <- st_T[id.Case1bis_boucle,]
      if(x$control$hazard_baseline_01 == "Splines"){
        Bs_L_R_i_01 <- Bs_L_R_01[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),]
        Bs_0_LR_i_01 <- Bs_0_LR_01[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1bis_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1bis_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),]
        }
      }
      if(x$control$hazard_baseline_02 == "Splines"){
        Bs_0_LR_i_02 <- Bs_0_LR_02[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1bis_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1bis_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),]
        }
      }
      if(x$control$hazard_baseline_02 == "Splines"){
        B_T_i_12 <- B_T_12[id.Case1bis_boucle,]
        Bs_T_i_12 <- Bs_T_12[(x$control$nb_pointsGK*(id.Case1bis_boucle-1)+1):(x$control$nb_pointsGK*(id.Case1bis_boucle)),]
        Bs_0_LR_i_12 <- Bs_0_LR_12[((x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case1bis_boucle-1))+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case1bis_boucle),]

      }

      X_base_i <- X_base[offset[id.Case1bis_boucle]:(offset[id.Case1bis_boucle+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[id.Case1bis_boucle+1]-offset[id.Case1bis_boucle])
      X_base_i <- unique(X_base_i)
      U_i <- U_base[offset[id.Case1bis_boucle]:(offset[id.Case1bis_boucle+1]-1),]
      U_i <- matrix(U_i, nrow = offset[id.Case1bis_boucle+1]-offset[id.Case1bis_boucle])
      U_base_i <- unique(U_i)
      y_i <- y.new[offset[id.Case1bis_boucle]:(offset[id.Case1bis_boucle+1]-1)]
      W_base_i <- W_base[offset[id.Case1bis_boucle]:(offset[id.Case1bis_boucle+1]-1),]
      W_base_i <- matrix(W_base_i, nrow = offset[id.Case1bis_boucle+1]-offset[id.Case1bis_boucle])
      W_base_i <- unique(W_base_i)
      O_base_i <- O_base[offset[id.Case1bis_boucle]:(offset[id.Case1bis_boucle+1]-1),]
      O_base_i <- matrix(O_base_i, nrow = offset[id.Case1bis_boucle+1]-offset[id.Case1bis_boucle])
      O_base_i <- unique(O_base_i)





      random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase1Bis, minimize = FALSE,nb.e.a = x$control$Objectlsmm$control$nb.e.a,
                                     nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma, Sigma.re = MatCov,
                                     sharedtype = sharedtype, HB = HB, W_G = W_G, nb_pointsGK = x$control$nb_pointsGK,
                                     alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega= omega,  wk = wk,
                                     delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i, Z_12_i=Z_12_i,
                                     X_T_i=X_T_i,  U_T_i=U_T_i, Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                     X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i, Uslope_GK_T_i=Uslope_GK_T_i,  O_GK_T_i=O_GK_T_i,  W_GK_T_i=W_GK_T_i,
                                     X_L_i=X_L_i,  U_L_i=U_L_i, Xslope_L_i=Xslope_L_i,  Uslope_L_i=Uslope_L_i, O_L_i=O_L_i,  W_L_i=W_L_i,
                                     X_GK_L_i=X_GK_L_i,  U_GK_L_i=U_GK_L_i,  Xslope_GK_L_i=Xslope_GK_L_i,Uslope_GK_L_i=Uslope_GK_L_i, O_GK_L_i=O_GK_L_i,  W_GK_L_i=W_GK_L_i,

                                     X_GK_T0_i=X_GK_T0_i,  U_GK_T0_i=U_GK_T0_i,  Xslope_GK_T0_i=Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i, O_GK_T0_i=O_GK_T0_i,  W_GK_T0_i=W_GK_T0_i,
                                     Time_T_i=Time_T_i,  Time_L_i=Time_L_i,  Time_T0_i=Time_T0_i, st_T_i=st_T_i,  st_L_i=st_L_i,  st_T0_i=st_T0_i,
                                     B_T_i_12=B_T_i_12,
                                     Bs_T_i_12=Bs_T_i_12,
                                     B_L_i_01=B_L_i_01,
                                     Bs_L_i_01=Bs_L_i_01,  Bs_L_i_02=Bs_L_i_02,  Bs_L_i_12=Bs_L_i_12,
                                     Bs_T0_i_01=Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc=x$control$left_trunc,
                                     X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, O_base_i=O_base_i,  W_base_i=W_base_i, index_b_slope = index_b_slope,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

      while(random.effects_i$istop !=1){
        binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
        random.effects_i <-  marqLevAlg(binit, fn = re_lsjm_covDepIDMCase1Bis, minimize = FALSE,nb.e.a = x$control$Objectlsmm$control$nb.e.a,
                                        nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma, Sigma.re = MatCov,
                                        sharedtype = sharedtype, HB = HB, W_G = W_G, nb_pointsGK = x$control$nb_pointsGK,
                                        alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega= omega,  wk = wk,
                                        delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i, Z_12_i=Z_12_i,
                                        X_T_i=X_T_i,  U_T_i=U_T_i, Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                        X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i, Uslope_GK_T_i=Uslope_GK_T_i,  O_GK_T_i=O_GK_T_i,  W_GK_T_i=W_GK_T_i,
                                        X_L_i=X_L_i,  U_L_i=U_L_i, Xslope_L_i=Xslope_L_i,  Uslope_L_i=Uslope_L_i, O_L_i=O_L_i,  W_L_i=W_L_i,
                                        X_GK_L_i=X_GK_L_i,  U_GK_L_i=U_GK_L_i,  Xslope_GK_L_i=Xslope_GK_L_i,Uslope_GK_L_i=Uslope_GK_L_i, O_GK_L_i=O_GK_L_i,  W_GK_L_i=W_GK_L_i,

                                        X_GK_T0_i=X_GK_T0_i,  U_GK_T0_i=U_GK_T0_i,  Xslope_GK_T0_i=Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i, O_GK_T0_i=O_GK_T0_i,  W_GK_T0_i=W_GK_T0_i,
                                        Time_T_i=Time_T_i,  Time_L_i=Time_L_i,  Time_T0_i=Time_T0_i, st_T_i=st_T_i,  st_L_i=st_L_i,  st_T0_i=st_T0_i,
                                        B_T_i_12=B_T_i_12,
                                        Bs_T_i_12=Bs_T_i_12,
                                        B_L_i_01=B_L_i_01,
                                        Bs_L_i_01=Bs_L_i_01,  Bs_L_i_02=Bs_L_i_02,  Bs_L_i_12=Bs_L_i_12,
                                        Bs_T0_i_01=Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc=x$control$left_trunc,
                                        X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, O_base_i=O_base_i,  W_base_i=W_base_i, index_b_slope = index_b_slope,
                                        nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                        file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
        binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a + x$control$Objectlsmm$control$nb.e.a.sigma)
      }

      random.effects.Predictions[id.Case1bis_boucle+nbCase1,] <-c( data.id.Case1bis$id[id.Case1bis_boucle], random.effects_i$b)
      CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
      time.measures_i <- time.measures[offset[id.Case1bis_boucle]:(offset[id.Case1bis_boucle+1]-1)]
      time.measures_i <- unique(time.measures_i)
      Varia <- exp(O_base_i%*%omega + W_base_i%*%random.effects_i$b[(x$control$Objectlsmm$control$nb.e.a+1):(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)])
      cv.Pred <- rbind(cv.Pred, cbind(rep(data.id.Case1bis$id[id.Case1bis_boucle], length(CV)),
                                      time.measures_i, CV, Varia))
    }

  }

  #4. Cas 2

  message("Case2")
  Case2 <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0); O_GK_T = as.matrix(0); W_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); O_T = as.matrix(0); W_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  O_GK_T0 = as.matrix(0); W_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)

  X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
  X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ;
  Xslope_T_i <- c(0); Uslope_T_i <- c(0);
  Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0);
  Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
  O_GK_T0_i <- as.matrix(0); W_GK_T0_i <- as.matrix(0); O_T_i <- c(0); W_T_i <- c(0);
  O_GK_T_i <- as.matrix(0); W_GK_T_i<- as.matrix(0) ;
  B_T_i_02 <- c(0);
  Bs_T_i_01 <- as.matrix(0) ;   Bs_T_i_02 <- as.matrix(0) ;
  Bs_T0_i_01<-as.matrix(0) ;   Bs_T0_i_02 <- as.matrix(0) ;
  st_T0_i <- c(0) ; st_T_i <- c(0)
  if(length(id.Case2)>0){
    data.long.Case2 <- data.long[which(data.long$id %in% id.Case2),]
    data.id.Case2 <- data.long.Case2[!duplicated(data.long.Case2$id),]
    list.long.Case2 <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long.Case2)
    X_base <- list.long.Case2$X; U_base <- list.long.Case2$U; y.new <- list.long.Case2$y.new
    offset <- list.long.Case2$offset

    list.var <- data.manag.sigma(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,data.id.Case2)
    O_base <- list.var$X
    O_base <- as.matrix(O_base)
    W_base <- list.var$U
    W_base <- as.matrix(W_base)

    time.measures <- data.long.Case2[,x$control$Objectlsmm$control$timeVar]

    list.GK_T <- data.GaussKronrod(data.id.Case2, a = 0, b = data.id.Case2$Time_T, k = x$control$nb_pointsGK)
    st_T <- list.GK_T$st
    if(x$control$left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case2, a = 0, b = data.id.Case2$Time_T0, k = x$control$nb_pointsGK)
      st_T0 <- list.GK_T0$st
    }
    if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) ){
      list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) ){
      list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) ){
      list.data_T <- data.time(data.id.Case2, data.id.Case2$Time_T, x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      O_T <- list.data_T$Xtime; W_T <- list.data_T$Utime
      O_GK_T <- list.data.GK_T$Xtime; W_GK_T <- list.data.GK_T$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
        O_GK_T0 <- list.data.GK_T0$Xtime
        W_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long.Case2)
    Z_01 <- list.surv$Z
    if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_02, data.long.Case2)
    Z_02 <- list.surv$Z
    if(x$control$hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    if(x$control$hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(x$control$knots_01, data.id.Case2$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(x$control$knots_01, c(t(st_T)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_01 <- splineDesign(x$control$knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(x$control$knots_02, data.id.Case2$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(x$control$knots_02, c(t(st_T)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_02 <- splineDesign(x$control$knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    for(id.Case2_boucle in 1:nbCase2){

      if("current value" %in% x$control$sharedtype_01 || "current value" %in% x$control$sharedtype_02){
        X_T_i <- X_T[id.Case2_boucle,];U_T_i <- U_T[id.Case2_boucle,]
        X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]);U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),])
        if(x$control$left_trunc){
          X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]);U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),])
        }
      }
      if("slope" %in% x$control$sharedtype_01 || "slope" %in%x$control$ sharedtype_02){
        Xslope_T_i <- Xslope_T[id.Case2_boucle,];Uslope_T_i <- Uslope_T[id.Case2_boucle,]
        Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]);Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),])
        if(x$control$left_trunc){
          Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]);Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),])
        }
      }

      if("variability" %in% x$control$sharedtype_01 || "variability" %in% x$control$sharedtype_02){
        O_T_i <- O_T[id.Case2_boucle,];W_T_i <- W_T[id.Case2_boucle,]
        O_GK_T_i <- as.matrix(O_GK_T[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]);W_GK_T_i <- as.matrix(W_GK_T[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),])
        if(x$control$left_trunc){
          O_GK_T0_i <- as.matrix(O_GK_T0[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]);W_GK_T0_i <- as.matrix(W_GK_T0[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),])
        }
      }
      if("Weibull" %in% c(x$control$hazard_baseline_01,x$control$hazard_baseline_02) ||"Gompertz" %in% c(x$control$hazard_baseline_01,x$control$hazard_baseline_02)){
        st_T_i <- st_T[id.Case2_boucle,]
        if(x$control$left_trunc){
          st_T0_i <- st_T0[id.Case2_boucle,]
        }
      }
      if("Splines" %in% x$control$hazard_baseline_01){
        Bs_T_i_01 <- Bs_T_01[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]
        }
      }
      if("Splines" %in% x$control$hazard_baseline_02){
        B_T_i_02 <- B_T_02[id.Case2_boucle,]
        Bs_T_i_02 <- Bs_T_02[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]
        if(left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(x$control$nb_pointsGK*(id.Case2_boucle-1)+1):(x$control$nb_pointsGK*id.Case2_boucle),]
        }
      }
      Z_01_i <- Z_01[id.Case2_boucle,]
      Z_02_i <- Z_02[id.Case2_boucle,]
      Time_T_i <- data.id.Case2$Time_T[id.Case2_boucle]
      if(x$control$left_trunc){
        Time_T0_i <- data.id.Case2$Time_T0[id.Case2_boucle]
      }
      delta2_i <- data.id.Case2$delta2[id.Case2_boucle]

      X_base_i <- X_base[offset[id.Case2_boucle]:(offset[id.Case2_boucle+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[id.Case2_boucle+1]-offset[id.Case2_boucle])
      X_base_i <- unique(X_base_i)
      U_i <- U_base[offset[id.Case2_boucle]:(offset[id.Case2_boucle+1]-1),]
      U_i <- matrix(U_i, nrow = offset[id.Case2_boucle+1]-offset[id.Case2_boucle])
      U_base_i <- unique(U_i)
      y_i <- y.new[offset[id.Case2_boucle]:(offset[id.Case2_boucle+1]-1)]

      W_base_i <- W_base[offset[id.Case2_boucle]:(offset[id.Case2_boucle+1]-1),]
      W_base_i <- matrix(W_base_i, nrow = offset[id.Case2_boucle+1]-offset[id.Case2_boucle])
      W_base_i <- unique(W_base_i)
      O_base_i <- O_base[offset[id.Case2_boucle]:(offset[id.Case2_boucle+1]-1),]
      O_base_i <- matrix(O_base_i, nrow = offset[id.Case2_boucle+1]-offset[id.Case2_boucle])
      O_base_i <- unique(O_base_i)


      random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase2, minimize = FALSE,

                                     nb.e.a = x$control$Objectlsmm$control$nb.e.a, nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma, Sigma.re = MatCov,
                                     sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                     alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,
                                     wk = wk,
                                     delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i,  X_T_i=X_T_i,  U_T_i=U_T_i,
                                     Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  O_T_i=O_T_i,  W_T_i=W_T_i,
                                     X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                     Uslope_GK_T_i=Uslope_GK_T_i, O_GK_T_i=O_GK_T_i,  W_GK_T_i=W_GK_T_i,
                                     X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                     O_GK_T0_i=O_GK_T0_i, W_GK_T0_i = W_GK_T0_i,
                                     Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                     B_T_i_02=B_T_i_02,
                                     Bs_T_i_01=Bs_T_i_01, Bs_T_i_02= Bs_T_i_02,
                                     Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc = x$control$left_trunc,
                                     X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, O_base_i=O_base_i,  W_base_i=W_base_i, index_b_slope = index_b_slope,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

      while(random.effects_i$istop !=1){
        binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
        random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase2, minimize = FALSE,

                                       nb.e.a = x$control$Objectlsmm$control$nb.e.a, nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma, Sigma.re = MatCov,
                                       sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                       alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,
                                       wk = wk,
                                       delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i,  X_T_i=X_T_i,  U_T_i=U_T_i,
                                       Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  O_T_i=O_T_i,  W_T_i=W_T_i,
                                       X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                       Uslope_GK_T_i=Uslope_GK_T_i, O_GK_T_i=O_GK_T_i,  W_GK_T_i=W_GK_T_i,
                                       X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                       O_GK_T0_i=O_GK_T0_i, W_GK_T0_i = W_GK_T0_i,
                                       Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                       B_T_i_02=B_T_i_02,
                                       Bs_T_i_01=Bs_T_i_01, Bs_T_i_02= Bs_T_i_02,
                                       Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc = x$control$left_trunc,
                                       X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, O_base_i=O_base_i,  W_base_i=W_base_i, index_b_slope = index_b_slope,
                                       nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                       file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
        binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)
      }

      random.effects.Predictions[id.Case2_boucle+nbCase1+nbCase1bis,] <- c(data.id.Case2$id[id.Case2_boucle],random.effects_i$b)
      time.measures_i <- time.measures[offset[id.Case2_boucle]:(offset[id.Case2_boucle+1]-1)]
      time.measures_i <- unique(time.measures_i)
      CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
      Varia <- exp(O_base_i%*%omega + W_base_i%*%random.effects_i$b[(x$control$Objectlsmm$control$nb.e.a+1):(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)])
      cv.Pred <- rbind(cv.Pred, cbind(rep(data.id.Case2$id[id.Case2_boucle], length(CV)),
                                      time.measures_i, CV, Varia))
    }

  }

  #5. Cas 3

  message("Case3")
  Case3 <- NULL
  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0); O_GK_T = as.matrix(0); W_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); O_T = as.matrix(0); W_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0); B_T_12 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Bs_T_12 = as.matrix(0); X_GK_L_T = as.matrix(0); U_GK_L_T = as.matrix(0);
  Xslope_GK_L_T = as.matrix(0); Uslope_GK_L_T = as.matrix(0); O_GK_L_T = as.matrix(0); W_GK_L_T = as.matrix(0);
  Time_L_T = c(0);
  st_L_T = as.matrix(0); B_GK_L_T_01=as.matrix(0); B_GK_L_T_02 = as.matrix(0); B_GK_L_T_12 = as.matrix(0); Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  O_GK_T0 = as.matrix(0); W_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0); X_0_LT = as.matrix(0); U_0_LT = as.matrix(0); Xslope_0_LT = as.matrix(0);
  Uslope_0_LT = as.matrix(0); Bs_0_LT_01 <- as.matrix(0);Bs_0_LT_02 <- as.matrix(0); Bs_0_LT_12 <- as.matrix(0)
  if(length(id.Case3)>0){
    data.long.Case3 <- data.long[which(data.long$id %in% id.Case3),]
    data.id.Case3 <- data.long.Case3[!duplicated(data.long.Case3$id),]
    list.long.Case3 <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long.Case3)
    X_base <- list.long.Case3$X; U_base <- list.long.Case3$U; y.new <- list.long.Case3$y.new
    offset <- list.long.Case3$offset

    list.var <- data.manag.sigma(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,data.long.Case3)
    O_base <- list.var$X
    O_base <- as.matrix(O_base)
    W_base <- list.var$U
    W_base <- as.matrix(W_base)

    time.measures <- data.long.Case3[,x$control$Objectlsmm$control$timeVar]

    list.GK_T <- data.GaussKronrod(data.id.Case3, a = 0, b = data.id.Case3$Time_T, k = x$control$nb_pointsGK)
    list.GK_L_T <- data.GaussKronrod(data.id.Case3, a = data.id.Case3$Time_L, b = data.id.Case3$Time_T, k = x$control$nb_pointsGK)
    st_T <- list.GK_T$st
    st_L_T <- list.GK_L_T$st
    if(x$control$left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id.Case3, a = 0, b = data.id.Case3$Time_T0, k = x$control$nb_pointsGK)
      st_T0 <- list.GK_T0$st
    }
    if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case3, data.id.Case3$Time_T, x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      list.data.GK_L_T <- data.time(list.GK_L_T$data.id2, c(t(st_L_T)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
      X_T <- list.data_T$Xtime; U_T <- list.data_T$Utime
      X_GK_T <- list.data.GK_T$Xtime; U_GK_T <- list.data.GK_T$Utime
      X_GK_L_T <- list.data.GK_L_T$Xtime; U_GK_L_T <- list.data.GK_L_T$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
        X_GK_T0 <- list.data.GK_T0$Xtime
        U_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case3, data.id.Case3$Time_T, x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      list.data.GK_L_T <- data.time(list.GK_L_T$data.id2, c(t(st_L_T)), x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
      Xslope_T <- list.data_T$Xtime; Uslope_T <- list.data_T$Utime
      Xslope_GK_T <- list.data.GK_T$Xtime; Uslope_GK_T <- list.data.GK_T$Utime
      Xslope_GK_L_T <- list.data.GK_L_T$Xtime; Uslope_GK_L_T <- list.data.GK_L_T$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
        Xslope_GK_T0 <- list.data.GK_T0$Xtime
        Uslope_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
      list.data_T <- data.time(data.id.Case3, data.id.Case3$Time_T, x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_T <- data.time(list.GK_T$data.id2, c(t(st_T)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      list.data.GK_L_T <- data.time(list.GK_L_T$data.id2, c(t(st_L_T)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
      O_T <- list.data_T$Xtime; W_T <- list.data_T$Utime
      O_GK_T <- list.data.GK_T$Xtime; W_GK_T <- list.data.GK_T$Utime
      O_GK_L_T <- list.data.GK_L_T$Xtime; W_GK_L_T <- list.data.GK_L_T$Utime

      if(x$control$left_trunc){
        list.data.GK_T0 <- data.time(list.GK_T0$data.id2, c(t(st_T0)),
                                     x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
        O_GK_T0 <- list.data.GK_T0$Xtime
        W_GK_T0 <- list.data.GK_T0$Utime
      }
    }

    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long.Case3)
    Z_01 <- list.surv$Z
    if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_02, data.long.Case3)
    Z_02 <- list.surv$Z
    if(x$control$hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_12, data.long.Case3)
    Z_12 <- list.surv$Z
    if(x$control$hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
    if(x$control$hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(x$control$knots_01, data.id.Case3$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(x$control$knots_01, c(t(st_T)), ord = 4L)
      B_GK_L_T_01 <- splineDesign(x$control$knots_01, c(t(st_L_T)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_01 <- splineDesign(x$control$knots_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(x$control$knots_02, data.id.Case3$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(x$control$knots_02, c(t(st_T)), ord = 4L)
      B_GK_L_T_02 <- splineDesign(x$control$knots_02, c(t(st_L_T)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_02 <- splineDesign(x$control$knots_02, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_12 == "Splines"){
      Z_12 <- as.matrix(Z_12[,-1])
      B_T_12 <- splineDesign(x$control$knots_12, data.id.Case3$Time_T, ord = 4L)
      Bs_T_12 <- splineDesign(x$control$knots_12, c(t(st_T)), ord = 4L)
      B_GK_L_T_12 <- splineDesign(x$control$knots_12, c(t(st_L_T)), ord = 4L)
    }

    ## Pour l'intégrale (à optmiser plus tard)
    st_0_LT <- c()
    if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
      X_0_LT <- c()
      U_0_LT <- c()
    }
    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
      Xslope_0_LT <- c()
      Uslope_0_LT <- c()
    }
    if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
      O_0_LT <- c()
      W_0_LT <- c()
    }

    if(x$control$hazard_baseline_01 == "Splines"){
      Bs_0_LT_01 <- c();
    }
    if(x$control$hazard_baseline_02== "Splines"){
      Bs_0_LT_02 <- c();
    }
    if(x$control$hazard_baseline_12== "Splines"){
      Bs_0_LT_12 <- c();
    }

    for(id.integrale in 1:nbCase3){
      # print(id.integrale)
      data.id.integrale <- data.id.Case3[id.integrale,]
      st_L_T_i <- st_L_T[id.integrale,]
      for(st.integrale in st_L_T_i){
        list.GK_0_stLT <- data.GaussKronrod(data.id.integrale, a = 0, b = st.integrale, k = x$control$nb_pointsGK)
        st_0_stLT_i <- list.GK_0_stLT$st
        st_0_LT <- rbind(st_0_LT, st_0_stLT_i)
        if(("current value" %in% x$control$sharedtype_01) || ("current value" %in% x$control$sharedtype_02) || ("current value" %in% x$control$sharedtype_12)){
          list.data.GK_0_stLT <- data.time(list.GK_0_stLT$data.id2, c(t(st_0_stLT_i)),x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
          X_0_stLT_i <- list.data.GK_0_stLT$Xtime; U_0_stLT_i <- list.data.GK_0_stLT$Utime
          X_0_LT <- rbind(X_0_LT,X_0_stLT_i); U_0_LT <- rbind(U_0_LT,U_0_stLT_i)
        }
        if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) || ("slope" %in% x$control$sharedtype_12)){
          list.data.GK_0_stLT <- data.time(list.GK_0_stLT$data.id2, c(t(st_0_stLT_i)),x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
          Xslope_0_stLT_i <- list.data.GK_0_stLT$Xtime; Uslope_0_stLT_i <- list.data.GK_0_stLT$Utime
          Xslope_0_LT <- rbind(Xslope_0_LT,Xslope_0_stLT_i); Uslope_0_LT <- rbind(Uslope_0_LT,Uslope_0_stLT_i)
        }
        if(("variability" %in% x$control$sharedtype_01) || ("variability" %in% x$control$sharedtype_02) || ("variability" %in% x$control$sharedtype_12)){
          list.data.GK_0_stLT <- data.time(list.GK_0_stLT$data.id2, c(t(st_0_stLT_i)),x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
          O_0_stLT_i <- list.data.GK_0_stLT$Xtime; W_0_stLT_i <- list.data.GK_0_stLT$Utime
          O_0_LT <- rbind(X_0_LT,X_0_stLT_i); W_0_LT <- rbind(U_0_LT,U_0_stLT_i)
        }
        if(x$control$hazard_baseline_01 == "Splines"){
          Bs_0_LT_01 <- rbind(Bs_0_LT_01,splineDesign(x$control$knots_01, c(t(st_0_stLT_i)), ord = 4L))
        }
        if(x$control$hazard_baseline_02 == "Splines"){
          Bs_0_LT_02 <- rbind(Bs_0_LT_02,splineDesign(x$control$knots_02, c(t(st_0_stLT_i)), ord = 4L))
        }
        if(x$control$hazard_baseline_12 == "Splines"){
          Bs_0_LT_12 <- rbind(Bs_0_LT_12,splineDesign(x$control$knots_12, c(t(st_0_stLT_i)), ord = 4L))
        }
      }



    }


    X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
    X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ; X_GK_L_T_i<- as.matrix(0) ; U_GK_L_T_i<- as.matrix(0) ; X_GK_0_LT_i<- as.matrix(0); U_GK_0_LT_i<- as.matrix(0);
    Xslope_T_i <- c(0); Uslope_T_i <- c(0); Xslope_GK_0_LT_i <- as.matrix(0); Uslope_GK_0_LT_i <- as.matrix(0);
    Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0); Xslope_GK_L_T_i <- as.matrix(0); Uslope_GK_L_T_i <- as.matrix(0)
    Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0)
    O_GK_T0_i <- as.matrix(0); W_GK_T0_i <- as.matrix(0); O_T_i <- c(0); W_T_i <- c(0);
    O_GK_T_i <- as.matrix(0); W_GK_T_i<- as.matrix(0) ; O_GK_L_T_i<- as.matrix(0) ; W_GK_L_T_i<- as.matrix(0) ; O_GK_0_LT_i<- as.matrix(0); W_GK_0_LT_i<- as.matrix(0);
    B_T_i_12 <- c(0); B_T_i_02 <- c(0)
    Bs_T_i_12 <- as.matrix(0); Bs_T_i_01 <- as.matrix(0);  Bs_T_i_02 <- as.matrix(0);
    Bs_0_LT_i_01<- as.matrix(0);  Bs_0_LT_i_02<- as.matrix(0);  Bs_0_LT_i_12<- as.matrix(0);
    Bs_L_T_i_01<- as.matrix(0);
    Bs_T0_i_01<- as.matrix(0);  Bs_T0_i_02<- as.matrix(0);
    st_T_i <- c(0); st_T0_i <- c(0); st_0_LT_i <- as.matrix(0);

    for(id.Case3_boucle in 1:nbCase3){
      #binit <- mvtnorm::rmvnorm(1, mean = rep(0, x$control$Objectlsmm$control$nb.e.a+2), MatCov)
      if("current value" %in% x$control$sharedtype_01 || "current value" %in% x$control$sharedtype_02 || "current value" %in% x$control$sharedtype_12){
        X_T_i <- X_T[id.Case3_boucle,];U_T_i <- U_T[id.Case3_boucle,]
        X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        X_GK_L_T_i <- as.matrix(X_GK_L_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);U_GK_L_T_i <- as.matrix(U_GK_L_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        #####()
        X_GK_0_LT_i <- as.matrix(X_0_LT[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),])
        U_GK_0_LT_i <- as.matrix(U_0_LT[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),])
        if(x$control$left_trunc){
          X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        }
      }
      if("slope" %in% x$control$sharedtype_01 || "slope" %in% x$control$sharedtype_02 || "slope" %in% x$control$sharedtype_12){
        Xslope_T_i <- Xslope_T[id.Case3_boucle,];Uslope_T_i <- Uslope_T[id.Case3_boucle,]
        Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        Xslope_GK_L_T_i <- as.matrix(Xslope_GK_L_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);Uslope_GK_L_T_i <- as.matrix(Uslope_GK_L_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        #####()
        Xslope_GK_0_LT_i <- as.matrix(Xslope_0_LT[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),])
        Uslope_GK_0_LT_i <- as.matrix(Uslope_0_LT[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),])
        if(x$control$left_trunc){
          Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        }
      }
      if("variability" %in% x$control$sharedtype_01 || "variability" %in% x$control$sharedtype_02 || "variability" %in% x$control$sharedtype_12){
        O_T_i <- O_T[id.Case3_boucle,];W_T_i <- W_T[id.Case3_boucle,]
        O_GK_T_i <- as.matrix(O_GK_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);W_GK_T_i <- as.matrix(W_GK_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        O_GK_L_T_i <- as.matrix(O_GK_L_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);W_GK_L_T_i <- as.matrix(W_GK_L_T[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        #####()
        O_GK_0_LT_i <- as.matrix(O_0_LT[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),])
        W_GK_0_LT_i <- as.matrix(W_0_LT[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),])
        if(x$control$left_trunc){
          O_GK_T0_i <- as.matrix(O_GK_T0[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]);W_GK_T0_i <- as.matrix(W_GK_T0[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),])
        }
      }
      if("Weibull" %in% c(x$control$hazard_baseline_01,x$control$hazard_baseline_02, x$control$hazard_baseline_12) ||"Gompertz" %in% c(x$control$hazard_baseline_01,x$control$hazard_baseline_02, x$control$hazard_baseline_12)){
        st_T_i <- st_T[id.Case3_boucle,]
        st_0_LT_i <- st_0_LT[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        st_L_T_i <- st_L_T[id.Case3_boucle,]
        if(x$control$left_trunc){
          st_T0_i <- st_T0[id.Case3_boucle,]
        }
      }
      if("Splines" %in% x$control$hazard_baseline_01){
        Bs_T_i_01 <- Bs_T_01[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        Bs_L_T_i_01 <- B_GK_L_T_01[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        Bs_0_LT_i_01 <- Bs_0_LT_01[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        }
      }
      if("Splines" %in% x$control$hazard_baseline_02){
        B_T_i_02 <- B_T_02[id.Case3_boucle,]
        Bs_T_i_02 <- Bs_T_02[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        Bs_0_LT_i_02 <- Bs_0_LT_02[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),]
        if(x$control$left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        }
      }
      if("Splines" %in% x$control$hazard_baseline_12){
        B_T_i_12 <- B_T_12[id.Case3_boucle,]
        Bs_T_i_12 <- Bs_T_12[(x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*id.Case3_boucle),]
        Bs_0_LT_i_12 <- Bs_0_LT_12[(x$control$nb_pointsGK*x$control$nb_pointsGK*(id.Case3_boucle-1)+1):(x$control$nb_pointsGK*x$control$nb_pointsGK*id.Case3_boucle),]
      }
      Z_01_i <- Z_01[id.Case3_boucle,]
      Z_02_i <- Z_02[id.Case3_boucle,]
      Z_12_i <- Z_12[id.Case3_boucle,]
      Time_T_i <- data.id.Case3$Time_T[id.Case3_boucle]
      Time_L_i <- data.id.Case3$Time_L[id.Case3_boucle]
      Time_L_T_i <- Time_T_i-Time_L_i
      ck <- ((sk_GK+1)/4)*Time_L_T_i+Time_L_i/2
      if(x$control$left_trunc){
        Time_T0_i <- data.id.Case3$Time_T0[id.Case3_boucle]
      }
      delta2_i <- data.id.Case3$delta2[id.Case3_boucle]

      X_base_i <- X_base[offset[id.Case3_boucle]:(offset[id.Case3_boucle+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[id.Case3_boucle+1]-offset[id.Case3_boucle])
      X_base_i <- unique(X_base_i)
      U_i <- U_base[offset[id.Case3_boucle]:(offset[id.Case3_boucle+1]-1),]
      U_i <- matrix(U_i, nrow = offset[id.Case3_boucle+1]-offset[id.Case3_boucle])
      U_base_i <- unique(U_i)
      y_i <- y.new[offset[id.Case3_boucle]:(offset[id.Case3_boucle+1]-1)]

      W_base_i <- W_base[offset[id.Case3_boucle]:(offset[id.Case3_boucle+1]-1),]
      W_base_i <- matrix(W_base_i, nrow = offset[id.Case3_boucle+1]-offset[id.Case3_boucle])
      W_base_i <- unique(W_base_i)
      O_base_i <- O_base[offset[id.Case3_boucle]:(offset[id.Case3_boucle+1]-1),]
      O_base_i <- matrix(O_base_i, nrow = offset[id.Case3_boucle+1]-offset[id.Case3_boucle])
      O_base_i <- unique(O_base_i)

      random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase3, minimize = FALSE,nb.e.a = x$control$Objectlsmm$control$nb.e.a, nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma,
                                     Sigma.re = MatCov,
                                     sharedtype = sharedtype, HB = HB, W_G = W_G, nb_pointsGK = x$control$nb_pointsGK,
                                     alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,
                                     wk = wk, rep_wk = rep_wk,
                                     delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i, Z_12_i=Z_12_i,
                                     X_T_i=X_T_i,  U_T_i=U_T_i, Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                     X_GK_T_i = X_GK_T_i,  U_GK_T_i = U_GK_T_i,  Xslope_GK_T_i = Xslope_GK_T_i,Uslope_GK_T_i = Uslope_GK_T_i, O_GK_T_i = O_GK_T_i,  W_GK_T_i = W_GK_T_i,
                                     X_GK_L_T_i = X_GK_L_T_i,  U_GK_L_T_i = U_GK_L_T_i,  Xslope_GK_L_T_i = Xslope_GK_L_T_i,  Uslope_GK_L_T_i = Uslope_GK_L_T_i, O_GK_L_T_i = O_GK_L_T_i,  W_GK_L_T_i = W_GK_L_T_i,
                                     X_GK_0_LT_i = X_GK_0_LT_i,  U_GK_0_LT_i = U_GK_0_LT_i,  Xslope_GK_0_LT_i = Xslope_GK_0_LT_i,  Uslope_GK_0_LT_i = Uslope_GK_0_LT_i, O_GK_0_LT_i = O_GK_0_LT_i,  W_GK_0_LT_i = W_GK_0_LT_i,
                                     X_GK_T0_i = X_GK_T0_i,  U_GK_T0_i = U_GK_T0_i,  Xslope_GK_T0_i = Xslope_GK_T0_i,  Uslope_GK_T0_i = Uslope_GK_T0_i, O_GK_T0_i = O_GK_T0_i,  W_GK_T0_i = W_GK_T0_i,
                                     Time_T_i = Time_T_i,  Time_L_T_i = Time_L_T_i,  Time_T0_i = Time_T0_i, st_T_i = st_T_i,  st_0_LT_i = st_0_LT_i,  st_L_T_i = st_L_T_i,  st_T0_i = st_T0_i,
                                     ck = ck,
                                     B_T_i_02 = B_T_i_02,  B_T_i_12 = B_T_i_12,
                                     Bs_T_i_01 = Bs_T_i_01,  Bs_T_i_02 = Bs_T_i_02,  Bs_T_i_12 = Bs_T_i_12,
                                     Bs_0_LT_i_01 = Bs_0_LT_i_01,  Bs_0_LT_i_02 = Bs_0_LT_i_02,  Bs_0_LT_i_12 = Bs_0_LT_i_12,
                                     Bs_L_T_i_01 = Bs_L_T_i_01,
                                     Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02 = Bs_T0_i_02,

                                     left_trunc= x$control$left_trunc,
                                     X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i, O_base_i = O_base_i, W_base_i = W_base_i,  index_b_slope = index_b_slope,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
      while(random.effects_i$istop !=1){
        binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
        random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepIDMCase3, minimize = FALSE,nb.e.a = x$control$Objectlsmm$control$nb.e.a, nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma,
                                       Sigma.re = MatCov,
                                       sharedtype = sharedtype, HB = HB, G_W = G_W, nb_pointsGK = x$control$nb_pointsGK,
                                       alpha_y_slope_var = alpha_y_slope_var, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,
                                       wk = wk, rep_wk = rep_wk,
                                       delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i, Z_12_i=Z_12_i,
                                       X_T_i=X_T_i,  U_T_i=U_T_i, Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                       X_GK_T_i = X_GK_T_i,  U_GK_T_i = U_GK_T_i,  Xslope_GK_T_i = Xslope_GK_T_i,Uslope_GK_T_i = Uslope_GK_T_i, O_GK_T_i = O_GK_T_i,  W_GK_T_i = W_GK_T_i,
                                       X_GK_L_T_i = X_GK_L_T_i,  U_GK_L_T_i = U_GK_L_T_i,  Xslope_GK_L_T_i = Xslope_GK_L_T_i,  Uslope_GK_L_T_i = Uslope_GK_L_T_i, O_GK_L_T_i = O_GK_L_T_i,  W_GK_L_T_i = W_GK_L_T_i,
                                       X_GK_0_LT_i = X_GK_0_LT_i,  U_GK_0_LT_i = U_GK_0_LT_i,  Xslope_GK_0_LT_i = Xslope_GK_0_LT_i,  Uslope_GK_0_LT_i = Uslope_GK_0_LT_i, O_GK_0_LT_i = O_GK_0_LT_i,  W_GK_0_LT_i = W_GK_0_LT_i,
                                       X_GK_T0_i = X_GK_T0_i,  U_GK_T0_i = U_GK_T0_i,  Xslope_GK_T0_i = Xslope_GK_T0_i,  Uslope_GK_T0_i = Uslope_GK_T0_i, O_GK_T0_i = O_GK_T0_i,  W_GK_T0_i = W_GK_T0_i,
                                       Time_T_i = Time_T_i,  Time_L_T_i = Time_L_T_i,  Time_T0_i = Time_T0_i, st_T_i = st_T_i,  st_0_LT_i = st_0_LT_i,  st_L_T_i = st_L_T_i,  st_T0_i = st_T0_i,
                                       ck = ck,
                                       B_T_i_02 = B_T_i_02,  B_T_i_12 = B_T_i_12,
                                       Bs_T_i_01 = Bs_T_i_01,  Bs_T_i_02 = Bs_T_i_02,  Bs_T_i_12 = Bs_T_i_12,
                                       Bs_0_LT_i_01 = Bs_0_LT_i_01,  Bs_0_LT_i_02 = Bs_0_LT_i_02,  Bs_0_LT_i_12 = Bs_0_LT_i_12,
                                       Bs_L_T_i_01 = Bs_L_T_i_01,
                                       Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02 = Bs_T0_i_02,

                                       left_trunc= x$control$left_trunc,
                                       X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i, O_base_i = O_base_i, W_base_i = W_base_i,  index_b_slope = index_b_slope,
                                       nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                       file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
        binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)
      }



      random.effects.Predictions[id.Case3_boucle+nbCase1+nbCase1bis+nbCase2,] <- c(data.id.Case3$id[id.Case3_boucle] ,random.effects_i$b)
      time.measures_i <- time.measures[offset[id.Case3_boucle]:(offset[id.Case3_boucle+1]-1)]
      time.measures_i <- unique(time.measures_i)
      CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
      Varia <- exp(O_base_i%*%omega + W_base_i%*%random.effects_i$b[(x$control$Objectlsmm$control$nb.e.a+1):(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)])
      cv.Pred <- rbind(cv.Pred, cbind(rep(data.id.Case3$id[id.Case3_boucle], length(CV)),
                                      time.measures_i, CV, Varia))
    }


  }



  pred_haz_01 <- 0
  pred_haz_02 <- 0
  pred_haz_12 <- 0
  data.id <- x$control$Objectlsmm$control$data.long[!duplicated(x$control$Objectlsmm$control$data.long$id),]
  gread.time <- seq(min(data.id[,x$control$Objectlsmm$control$timeVar]), max(data.id[,x$control$Objectlsmm$control$timeVar]), by = (max(data.id[,x$control$Objectlsmm$control$timeVar])-min(data.id[,x$control$Objectlsmm$control$timeVar]))/100)
  data.GaussKronrod.sort.unique <- data.GaussKronrod(data.id = data.id, a = 0,b = gread.time, k = x$control$nb_pointsGK)
  st_calc.sort.unique <- data.GaussKronrod.sort.unique$st
  P.sort.unique <- data.GaussKronrod.sort.unique$P

  Cum_risk_01 <- c()
  Cum_risk_02 <- c()
  Cum_risk_12 <- c()
  wk <- gaussKronrod()$wk

  list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long)
  Z_01 <- list.surv$Z
  if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
  list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_02, data.long)
  Z_02 <- list.surv$Z
  if(x$control$hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
  list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_12, data.long)
  Z_12 <- list.surv$Z
  if(x$control$hazard_baseline_12 == "Gompertz"){Z_12 <- as.matrix(Z_12[,-1])}
  if(x$control$hazard_baseline_01 == "Splines"){
    Z_01 <- as.matrix(Z_01[,-1])
  #  B_T_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, data.id$Time_T, ord = 4L)
  #  Bs_T_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, c(t(st_T)), ord = 4L)
  #  if(x$control$left_trunc){
  #    Bs_T0_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, c(t(st_T0)), ord = 4L)
  #  }
  }
  if(x$control$hazard_baseline_02 == "Splines"){
    Z_02 <- as.matrix(Z_02[,-1])
  #  B_T_02 <- splineDesign(x$control$knots.hazard_baseline.splines_02, data.id$Time_T, ord = 4L)
  #  Bs_T_02 <- splineDesign(x$control$knots.hazard_baseline.splines_02, c(t(st_T)), ord = 4L)
  #  if(x$control$left_trunc){
  #    Bs_T0_02 <- splineDesign(x$control$knots.hazard_baseline.splines_02, c(t(st_T0)), ord = 4L)
  #  }
  }
  if(x$control$hazard_baseline_12 == "Splines"){
    Z_12 <- as.matrix(Z_12[,-1])
  #  B_T_12 <- splineDesign(x$control$knots.hazard_baseline.splines_12, data.id$Time_T, ord = 4L)
  #  Bs_T_12 <- splineDesign(x$control$knots.hazard_baseline.splines_12, c(t(st_T)), ord = 4L)
  #  if(x$control$left_trunc){
  #    Bs_T0_12 <- splineDesign(x$control$knots.hazard_baseline.splines_12, c(t(st_T0)), ord = 4L)
  #  }
  }

  cat("Cumulative risks")
  for(id_boucle in 1:nrow(data.id)){
    Cum_risk_01i <- c()
    Cum_risk_02i  <- c()
    Cum_risk_12i <- c()
    for(j in 1:length((gread.time))){
      pred_haz_01 <- 0
      pred_haz_02 <- 0
      pred_haz_12 <- 0
      if("variability" %in% x$control$sharedtype_01 || "variability" %in% x$control$sharedtype_02){
        list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),], st_calc.sort.unique[j,],
                                                            x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
        Os.j <- list.data.GK.current.sigma.sort.unique$Xtime
        Ws.j <- list.data.GK.current.sigma.sort.unique$Utime
        Sigma.current.GK <- exp(omega%*%t(Os.j) + random.effects_i$b[(x$control$Objectlsmm$control$nb.e.a+1):(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)]%*%t(Ws.j))
        if("variability" %in% x$control$sharedtype_01){
          pred_haz_01 <- pred_haz_01 + alpha.var_01*Sigma.current.GK
        }
        if("variability" %in% x$control$sharedtype_02){
          pred_haz_02 <- pred_haz_02 + alpha.var_02*Sigma.current.GK
        }
        if("variability" %in% x$control$sharedtype_12){
          pred_haz_12 <- pred_haz_12 + alpha.var_12*Sigma.current.GK
        }
      }

      if("current value" %in% x$control$sharedtype_01 || "current value" %in% x$control$sharedtype_02){
        list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),], st_calc.sort.unique[j,],
                                                            x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
        Xs.j <- list.data.GK.current.sigma.sort.unique$Xtime
        Us.j <- list.data.GK.current.sigma.sort.unique$Utime
        current.GK <- beta%*%t(Xs.j) + random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]%*%t(Us.j)
        if("current value" %in% x$control$sharedtype_01){
          pred_haz_01 <- pred_haz_01 + alpha.current_01*current.GK
        }
        if("current value" %in% x$control$sharedtype_02){
          pred_haz_02 <- pred_haz_02 + alpha.current_02*current.GK
        }
        if("current value" %in% x$control$sharedtype_12){
          pred_haz_12 <- pred_haz_12 + alpha.current_12*current.GK
        }
      }

      if("slope" %in% x$control$sharedtype_01 || "slope" %in% x$control$sharedtype_02){
        list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),], st_calc.sort.unique[j,],
                                                            x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
        Xs.slope.j <- list.data.GK.current.sigma.sort.unique$Xtime
        Us.slope.j <- list.data.GK.current.sigma.sort.unique$Utime
        bslope <- random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
        bslope <- bslope[x$control$index_b_slope]
        bslope <- as.matrix(bslope, nrow = 1)
        slope.GK <- beta_slope%*%t(Xs.slope.j) + bslope%*%t(Us.slope.j)
        if("slope" %in% x$control$sharedtype_01){
          pred_haz_01 <- pred_haz_01 + alpha.slope_01*slope.GK
        }
        if("slope" %in% x$control$sharedtype_02){
          pred_haz_02 <- pred_haz_02 + alpha.slope_02*slope.GK
        }
        if("slope" %in% x$control$sharedtype_12){
          pred_haz_12 <- pred_haz_02 + alpha.slope_12*slope.GK
        }
      }

      if(x$control$hazard_baseline_01 == "Exponential"){
        h_0_01 <- 1
        h_0.GK_01 <- wk
      }
      if(x$control$hazard_baseline_01 == "Weibull"){
        st_j <- st_calc.sort.unique[j,]
        h_0.GK_01 <- shape_01*(st_j**(shape_01-1))*wk    #### AJOUTER GOMPERTZ
      }
      if(x$control$hazard_baseline_01 == "Gompertz"){
        stop("Not implemented.")    #### AJOUTER GOMPERTZ
      }
      if(x$control$hazard_baseline_01 == "Splines"){
        st_j <- st_calc.sort.unique[j,]
        Bs_j <- splines::splineDesign(x$control$knots.hazard_baseline.splines_01, st_j, ord = 4L)
        #Bs_j <- Bs[(x$control$nb_pointsGK*(j-1)+1):(x$control$nb_pointsGK*j),]
        mat_h0s <- matrix(gamma_01,ncol=1)
        h_0.GK_01 <- (wk*exp(Bs_j%*%mat_h0s))
      }

      Z_01_i <- Z_01[id_boucle,]
      if(length(Z_01_i)==0){
        pred_surv_01 <- 0
      }
      else{
        pred_surv_01 <- (alpha_01%*%Z_01_i)[1,1]
      }

      pred_haz_01 <- pred_haz_01 + pred_surv_01

      Cum_risk_01i <- c(Cum_risk_01i, P.sort.unique[j]*sum(exp(pred_haz_01)%*%h_0.GK_01))

      if(x$control$hazard_baseline_02 == "Exponential"){
        h_0_02 <- 1
        h_0.GK_02 <- wk
      }
      if(x$control$hazard_baseline_02 == "Weibull"){
        st_j <- st_calc.sort.unique[j,]
        h_0.GK_02 <- shape_02*(st_j**(shape_02-1))*wk    #### AJOUTER GOMPERTZ
      }
      if(x$control$hazard_baseline_02 == "Gompertz"){
        stop("Not implemented.")    #### AJOUTER GOMPERTZ
      }
      if(x$control$hazard_baseline_02 == "Splines"){
        st_j <- st_calc.sort.unique[j,]
        Bs_j <- splines::splineDesign(x$control$knots.hazard_baseline.splines_02, st_j, ord = 4L)
        #Bs_j <- Bs[(x$control$nb_pointsGK*(j-1)+1):(x$control$nb_pointsGK*j),]
        mat_h0s <- matrix(gamma_02,ncol=1)
        h_0.GK_02 <- (wk*exp(Bs_j%*%mat_h0s))
      }

      Z_02_i <- Z_02[id_boucle,]
      if(length(Z_02_i)==0){
        pred_surv_02 <- 0
      }
      else{
        pred_surv_02 <- (alpha_02%*%Z_02_i)[1,1]
      }

      pred_haz_02 <- pred_haz_02 + pred_surv_02

      Cum_risk_02i <- c(Cum_risk_02i, P.sort.unique[j]*sum(exp(pred_haz_02)%*%h_0.GK_02))

      if(x$control$hazard_baseline_12 == "Exponential"){
        h_0_12 <- 1
        h_0.GK_12 <- wk
      }
      if(x$control$hazard_baseline_12 == "Weibull"){
        st_j <- st_calc.sort.unique[j,]
        h_0.GK_12 <- shape_12*(st_j**(shape_12-1))*wk    #### AJOUTER GOMPERTZ
      }
      if(x$control$hazard_baseline_12 == "Gompertz"){
        stop("Not implemented.")    #### AJOUTER GOMPERTZ
      }
      if(x$control$hazard_baseline_12 == "Splines"){
        st_j <- st_calc.sort.unique[j,]
        Bs_j <- splines::splineDesign(x$control$knots.hazard_baseline.splines_12, st_j, ord = 4L)
        #Bs_j <- Bs[(x$control$nb_pointsGK*(j-1)+1):(x$control$nb_pointsGK*j),]
        mat_h0s <- matrix(gamma_12,ncol=1)
        h_0.GK_12 <- (wk*exp(Bs_j%*%mat_h0s))
      }

      Z_12_i <- Z_12[id_boucle,]
      if(length(Z_12_i)==0){
        pred_surv_12 <- 0
      }
      else{
        pred_surv_12 <- (alpha_12%*%Z_12_i)[1,1]
      }

      pred_haz_12 <- pred_haz_12 + pred_surv_12

      Cum_risk_12i <- c(Cum_risk_12i, P.sort.unique[j]*sum(exp(pred_haz_12)%*%h_0.GK_12))

    }
    Cum_risk_01 <- rbind(Cum_risk_01,Cum_risk_01i)
    Cum_risk_02 <- rbind(Cum_risk_02,Cum_risk_02i)
    Cum_risk_12 <- rbind(Cum_risk_12,Cum_risk_12i)
  }

  Cum_01 <- apply(Cum_risk_01,2,mean)
  Cum_02 <- apply(Cum_risk_02,2,mean)
  Cum_12 <- apply(Cum_risk_12,2,mean)




  cv.Pred <- as.data.frame(cv.Pred)
  colnames(cv.Pred) <- c("id", "time", "CV", "Residual_SD")
  list(random.effects.Predictions = random.effects.Predictions, cv.Pred = cv.Pred, Cum_01 = Cum_01, Cum_02 = Cum_02, Cum_12 = Cum_12)

}
