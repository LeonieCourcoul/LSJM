#' @importFrom splines splineDesign
#' @importFrom stats model.frame model.matrix qnorm dnorm
#' @importFrom spacefillr generate_sobol_owen_set
predyn_ponct_lsjm_covDepIDM <- function(Objectlsjm, data.long.until.time.s, s, window, event){


  if(is.null(Objectlsjm$result_step2)){
    param <- Objectlsjm$result_step1$b
    nbQMC <- Objectlsjm$control$S1
  }
  else{
    param <- Objectlsjm$result_step2$b
    nbQMC <- Objectlsjm$control$S2
  }

  ## Param
  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(Objectlsjm$control$hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]**2
    curseur <- curseur + 1
  }

  if(Objectlsjm$control$hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]**2
    Gompertz.2_01 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(Objectlsjm$control$hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+Objectlsjm$control$nb.knots.splines[1]+2+1)]
    curseur <- curseur + Objectlsjm$control$nb.knots.splines[1]+2 + 2
  }
  ### Covariables :
  nb.alpha_01 <- Objectlsjm$control$nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("random effects" %in% Objectlsjm$control$sharedtype_01){
    alpha_b_01 <- param[curseur:(curseur+Objectlsjm$control$Objectlsmm$control$nb.e.a-1)]
    curseur <- curseur + Objectlsjm$control$Objectlsmm$control$nb.e.a
  }
  if("value" %in% Objectlsjm$control$sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% Objectlsjm$control$sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% Objectlsjm$control$sharedtype_01){
    alpha.var_01 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Risque 02
  if(Objectlsjm$control$hazard_baseline_02 == "Weibull"){
    shape_02 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(Objectlsjm$control$hazard_baseline_02 == "Gompertz"){
    Gompertz.1_02 <- param[curseur]**2
    Gompertz.2_02 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(Objectlsjm$control$hazard_baseline_02 == "Splines"){
    gamma_02 <- param[(curseur):(curseur+Objectlsjm$control$nb.knots.splines[2]+2+1)]
    curseur <- curseur + Objectlsjm$control$nb.knots.splines[2]+2+ 2
  }
  ### Covariables :
  nb.alpha_02 <- Objectlsjm$control$nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("random effects" %in% Objectlsjm$control$sharedtype_02){
    alpha_b_02 <- param[curseur:(curseur+Objectlsjm$control$Objectlsmm$control$nb.e.a-1)]
    curseur <- curseur + Objectlsjm$control$Objectlsmm$control$nb.e.a
  }
  if("value" %in% Objectlsjm$control$sharedtype_02){
    alpha.current_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% Objectlsjm$control$sharedtype_02){
    alpha.slope_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% Objectlsjm$control$sharedtype_02){
    alpha.var_02 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Risque 12
  ### Hazard baseline
  if(Objectlsjm$control$hazard_baseline_12 == "Weibull"){
    shape_12 <- param[curseur]**2
    curseur <- curseur + 1
  }

  if(Objectlsjm$control$hazard_baseline_12 == "Gompertz"){
    Gompertz.1_12 <- param[curseur]**2
    Gompertz.2_12 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(Objectlsjm$control$hazard_baseline_12 == "Splines"){
    gamma_12 <- param[(curseur):(curseur+Objectlsjm$control$nb.knots.splines[3]+2+1)]
    curseur <- curseur + Objectlsjm$control$nb.knots.splines[3]+2 + 2
  }
  ### Covariables :
  nb.alpha_12 <- Objectlsjm$control$nb.alpha[3]
  if(nb.alpha_12 >=1){
    alpha_12 <-  param[(curseur):(curseur+nb.alpha_12-1)]
    curseur <- curseur+nb.alpha_12
  }
  ### Association
  if("value" %in% Objectlsjm$control$sharedtype_12){
    alpha.current_12 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% Objectlsjm$control$sharedtype_12){
    alpha.slope_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% Objectlsjm$control$sharedtype_12){
    alpha.var_12 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ Objectlsjm$control$Objectlsmm$control$nb.beta-1)]
  curseur <- curseur+Objectlsjm$control$Objectlsmm$control$nb.beta
  omega <- param[(curseur):(curseur+Objectlsjm$control$Objectlsmm$control$nb.omega-1)]
  curseur <- curseur+Objectlsjm$control$Objectlsmm$control$nb.omega


  Zq1 <- generate_sobol_owen_set(nbQMC,  Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)
  Zq <- apply(Zq1, 2, qnorm)

  if(Objectlsjm$control$Objectlsmm$control$correlated_re){
    C1 <- matrix(rep(0,(Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)**2),nrow=Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma,ncol=Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)
    C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
    Cholesky <- C1
    Cholesky <- as.matrix(Cholesky)
    random.effects <- Zq%*%t(Cholesky)
    b_al <- random.effects[,1:Objectlsjm$control$Objectlsmm$control$nb.e.a]
    b_al <- matrix(b_al, ncol = Objectlsjm$control$Objectlsmm$control$nb.e.a)
    b_om <- random.effects[,(Objectlsjm$control$Objectlsmm$control$nb.e.a+1):(Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)]
    b_om <- matrix(b_om, ncol = Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)
  }
  else{
    borne1 <- curseur + choose(n = Objectlsjm$control$Objectlsmm$control$nb.e.a, k = 2) + Objectlsjm$control$Objectlsmm$control$nb.e.a - 1
    C1 <- matrix(rep(0,(Objectlsjm$control$Objectlsmm$control$nb.e.a)**2),nrow=Objectlsjm$control$Objectlsmm$control$nb.e.a,ncol=Objectlsjm$control$Objectlsmm$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    borne3 <- borne1 + choose(n = Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma, k = 2) + Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma
    C3 <- matrix(rep(0,(Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)**2),nrow=Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma,ncol=Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
    MatCovb <- as.matrix(C1)
    MatCovSig <- as.matrix(C3)
    b_al <- Zq[,1:Objectlsjm$control$Objectlsmm$control$nb.e.a]%*%t(MatCovb)
    b_al <- matrix(b_al, ncol = Objectlsjm$control$Objectlsmm$control$nb.e.a)
    b_om <- Zq[,(Objectlsjm$control$Objectlsmm$control$nb.e.a+1):(Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)]%*%t(MatCovSig)
    b_om <- matrix(b_om, ncol = Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)
  }

  if( "slope" %in%  Objectlsjm$control$sharedtype_01 || "slope" %in%  Objectlsjm$control$sharedtype_02){
    beta_slope <- beta[Objectlsjm$control$index_beta_slope]
    b_al_slope <- matrix(b_al[,Objectlsjm$control$index_b_slope], ncol = length(Objectlsjm$control$index_b_slope))
  }

  data.long.until.time.s.id <- data.long.until.time.s[1,]

  ### Longitudinal part
  list.long <- data.manag.long(Objectlsjm$control$Objectlsmm$control$formGroup,Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,data.long.until.time.s)
  X_base <- list.long$X
  U <- list.long$U
  y.new.prog <- list.long$y.new.prog

  list.var <- data.manag.sigma(Objectlsjm$control$Objectlsmm$control$formGroup,Objectlsjm$control$Objectlsmm$control$formFixedVar, Objectlsjm$control$Objectlsmm$control$formRandomVar,data.long.until.time.s)
  O_base <- list.var$X
  W_base <- list.var$U
  if(is.null(nrow(X_base))){
    sigma.long <- exp((omega%*%O_base)[1,1] + b_om%*%W_base)
    CV <- (beta%*%X_base)[1,1] + b_al%*%U
    f_Y_b_sigma <- dnorm(x=y.new.prog, mean = CV, sd = sigma.long)
  }else{
    f_Y_b_sigma <- rep(1,nbQMC)
    for(k in 1:nrow(X_base)){
      sigma.long <- exp((omega%*%O_base[k,])[1,1] + b_om%*%W_base[k,])
      CV <- (beta%*%X_base[k,])[1,1] + b_al%*%U[k,]
      f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y.new.prog[k], mean = CV, sd = sigma.long)
    }
  }

  # Survival part

  data.GaussKronrod.1 <- data.GaussKronrod(data.long.until.time.s.id,a=s,b=s+window,k = Objectlsjm$control$nb_pointsGK)
  P.1 <- data.GaussKronrod.1$P
  st.1 <- data.GaussKronrod.1$st
  wk.1 <- data.GaussKronrod.1$wk
  data.id.1 <- data.GaussKronrod.1$data.id2
  data.GaussKronrod.den <- data.GaussKronrod(data.long.until.time.s.id,a = 0,b=s,k = Objectlsjm$control$nb_pointsGK)
  P.den <- data.GaussKronrod.den$P
  st.den <- data.GaussKronrod.den$st
  wk.den <- data.GaussKronrod.den$wk
  data.id.den <- data.GaussKronrod.den$data.id2

  if("value" %in% Objectlsjm$control$sharedtype_01 ||"value" %in% Objectlsjm$control$sharedtype_02){
    list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                       Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs <- list.data.GK.current$Xtime
    Us <- list.data.GK.current$Utime

    list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                           Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs.den <- list.data.GK.current.den$Xtime
    Us.den <- list.data.GK.current.den$Utime

  }

  if("slope" %in% Objectlsjm$control$sharedtype_01 ||"slope" %in% Objectlsjm$control$sharedtype_02){
    list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                       Objectlsjm$control$formSlopeFixed, Objectlsjm$control$formSlopeRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs.slope <- list.data.GK.current$Xtime
    Us.slope <- list.data.GK.current$Utime

    list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                           Objectlsjm$control$formSlopeFixed, Objectlsjm$control$formSlopeRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs.slope.den <- list.data.GK.current.den$Xtime
    Us.slope.den <- list.data.GK.current.den$Utime

  }

  if("variability" %in% Objectlsjm$control$sharedtype_01 ||"variability" %in% Objectlsjm$control$sharedtype_02){
    list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                       Objectlsjm$control$Objectlsmm$control$formFixedVar, Objectlsjm$control$Objectlsmm$control$formRandomVar,Objectlsjm$control$Objectlsmm$control$timeVar)
    Os <- list.data.GK.current$Xtime
    Ws <- list.data.GK.current$Utime

    list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                           Objectlsjm$control$Objectlsmm$control$formFixedVar, Objectlsjm$control$Objectlsmm$control$formRandomVar,Objectlsjm$control$Objectlsmm$control$timeVar)
    Os.den <- list.data.GK.current.den$Xtime
    Ws.den <- list.data.GK.current.den$Utime

  }

  if(Objectlsjm$control$hazard_baseline_01 == "Exponential"){
    mfZ <- model.frame(Objectlsjm$control$formSurv_01, data = data.long.until.time.s.id)
    Z_01 <- model.matrix(Objectlsjm$control$formSurv_01, mfZ)
  }else{
    if(Objectlsjm$control$hazard_baseline_01 == "Weibull" || Objectlsjm$control$hazard_baseline_01 == "Gompertz"){
      mfZ <- model.frame(Objectlsjm$control$formSurv_01, data = data.long.until.time.s.id)
      Z_01 <- model.matrix(Objectlsjm$control$formSurv_01, mfZ)
    }else{
      if(Objectlsjm$control$hazard_baseline_01 == "Splines"){
        mfZ <- model.frame(Objectlsjm$control$formSurv_01, data = data.long.until.time.s.id)
        Z_01 <- model.matrix(Objectlsjm$control$formSurv_01, mfZ)
        Z_01 <- Z_01[,-1]
        Bs_01 <- splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_01, c(t(st.1)), ord = 4L)
        Bs.den_01 <- splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_01, c(t(st.den)), ord = 4L)
      }else{
        stop("This type of base survival function is not implemented.")
      }
    }
  }

  if(Objectlsjm$control$hazard_baseline_02 == "Exponential"){
    mfZ <- model.frame(Objectlsjm$control$formSurv_02, data = data.long.until.time.s.id)
    Z_02 <- model.matrix(Objectlsjm$control$formSurv_02, mfZ)
  }else{
    if(Objectlsjm$control$hazard_baseline_02 == "Weibull" || Objectlsjm$control$hazard_baseline_02 == "Gompertz"){
      mfZ <- model.frame(Objectlsjm$control$formSurv_02, data = data.long.until.time.s.id)
      Z_02 <- model.matrix(Objectlsjm$control$formSurv_02, mfZ)
    }else{
      if(Objectlsjm$control$hazard_baseline_02 == "Splines"){
        mfZ <- model.frame(Objectlsjm$control$formSurv_02, data = data.long.until.time.s.id)
        Z_02 <- model.matrix(Objectlsjm$control$formSurv_02, mfZ)
        Z_02 <- Z_02[,-1]
        Bs_02 <- splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_02, c(t(st.1)), ord = 4L)
        Bs.den_02 <- splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_02, c(t(st.den)), ord = 4L)
      }else{
        stop("This type of base survival function is not implemented.")
      }
    }
  }

  ### Between 0 and u (double integral)
  st_0_u <- c(); X_0_u <- c(); U_0_u <- c(); Xslope_0_u <- c(); Uslope_0_u <- c(); Bs_0_u_01 <- c(); Bs_0_u_02 <- c(); O_0_u <- c(); W_0_u <- c()
  data.id.integrale <- data.long.until.time.s.id
  for(st.integrale in st.1){
    list.GK_0_st.2 <- data.GaussKronrod(data.id.integrale, a= 0, b = st.integrale, k = Objectlsjm$control$nb_pointsGK)
    st.2 <- list.GK_0_st.2$st
    st_0_u <- rbind(st_0_u, st.2)
    if(("variability" %in% Objectlsjm$control$sharedtype_01) || ("variability" %in% Objectlsjm$control$sharedtype_02)){
      list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),Objectlsjm$control$Objectlsmm$control$formFixedVar, Objectlsjm$control$Objectlsmm$control$formRandomVar,Objectlsjm$control$Objectlsmm$control$timeVar)
      O_0_st_u <- list.data.GK_0_u$Xtime; W_0_st_u <- list.data.GK_0_u$Utime
      O_0_u <- rbind(O_0_u,O_0_st_u); W_0_u <- rbind(W_0_u,W_0_st_u)
    }
    if(("value" %in% Objectlsjm$control$sharedtype_01) || ( "value" %in% Objectlsjm$control$sharedtype_02)){
      list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
      X_0_st_u <- list.data.GK_0_u$Xtime; U_0_st_u <- list.data.GK_0_u$Utime
      X_0_u <- rbind(X_0_u,X_0_st_u); U_0_u <- rbind(U_0_u,U_0_st_u)
    }
    if(("slope" %in% Objectlsjm$control$sharedtype_01) || ("slope" %in% Objectlsjm$control$sharedtype_02)){
      list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),Objectlsjm$control$formSlopeFixed, Objectlsjm$control$formSlopeRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
      Xslope_0_st_u <- list.data.GK_0_u$Xtime; Uslope_0_st_u <- list.data.GK_0_u$Utime
      Xslope_0_u <- rbind(Xslope_0_u,Xslope_0_st_u); Uslope_0_u <- rbind(Uslope_0_u,Uslope_0_st_u)
    }
    if(Objectlsjm$control$hazard_baseline_01 == "Splines"){
      Bs_0_u_01 <- rbind(Bs_0_u_01,splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_01, c(t(st.2)), ord = 4L))
    }
    if(Objectlsjm$control$hazard_baseline_02 == "Splines"){
      Bs_0_u_02 <- rbind(Bs_0_u_02,splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_02, c(t(st.2)), ord = 4L))
    }
  }

  # Computation
  etaBaseline_s_t_0k <- 0; survLong_s_t_0k <- 0; etaBaseline_0_s_01 <- 0; survLong_0_s_01 <- 0;
  etaBaseline_0_s_02 <- 0; survLong_0_s_02 <- 0; etaBaseline_0_u_01 <- 0; survLong_0_u_01 <- 0;
  etaBaseline_0_u_02 <- 0; survLong_0_u_02 <- 0


  if(c("random effects") %in% Objectlsjm$control$sharedtype_01){
    survLong_0_s_01 <- survLong_0_s_01 + c(b_al%*%alpha_b_01)
    survLong_0_u_01 <- survLong_0_u_01 + c(b_al%*%alpha_b_01)
    if(event == 1){
      survLong_s_t_0k <- survLong_s_t_0k + c(b_al%*%alpha_b_01)
    }
  }
  if(c("random effects") %in% Objectlsjm$control$sharedtype_02){
    survLong_0_s_02 <- survLong_0_s_02 + c(b_al%*%alpha_b_02)
    survLong_0_u_02 <- survLong_0_u_02 + c(b_al%*%alpha_b_02)
    if(event == 2){
      survLong_s_t_0k <- survLong_s_t_0k + c(b_al%*%alpha_b_02)
    }
  }


  if((c("value") %in% Objectlsjm$control$sharedtype_01 )|| (c("value") %in% Objectlsjm$control$sharedtype_02)){
    current.GK <- matrix(rep(beta%*%t(Xs),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(Us)
    current.GK.den <- matrix(rep(beta%*%t(Xs.den),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(Us.den)
    current.GK.0_u <- matrix(rep(beta%*%t(X_0_u),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(U_0_u)
    if(c("value") %in% Objectlsjm$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + alpha.current_01*current.GK.den
      survLong_0_u_01 <- survLong_0_u_01 + alpha.current_01*current.GK.0_u
      if(event == 1){
        survLong_s_t_0k <- survLong_s_t_0k + alpha.current_01*current.GK
      }
    }
    if(c("value") %in% Objectlsjm$control$sharedtype_02){
      survLong_0_s_02 <- survLong_0_s_02 + alpha.current_02*current.GK.den
      survLong_0_u_02 <- survLong_0_u_02 + alpha.current_02*current.GK.0_u
      if(event == 2){
        survLong_s_t_0k <- survLong_s_t_0k + alpha.current_02*current.GK
      }
    }
  }

  if((c("slope") %in% Objectlsjm$control$sharedtype_01 )|| (c("slope") %in% Objectlsjm$control$sharedtype_02)){
    slope.GK <- matrix(rep(beta_slope%*%t(Xs.slope),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Us.slope)
    slope.GK.den <- matrix(rep(beta_slope%*%t(Xs.slope.den),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Us.slope.den)
    slope.GK.0_u <- matrix(rep(beta_slope%*%t(Xslope_0_u),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Uslope_0_u)
    if(c("slope") %in% Objectlsjm$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + alpha.slope_01*slope.GK.den
      survLong_0_u_01 <- survLong_0_u_01 + alpha.slope_01*slope.GK.0_u
      if(event == 1){
        survLong_s_t_0k <- survLong_s_t_0k + alpha.slope_01*slope.GK
      }
    }
    if(c("slope") %in% Objectlsjm$control$sharedtype_02){
      survLong_0_s_02 <- survLong_0_s_02 + alpha.slope_02*slope.GK.den
      survLong_0_u_02 <- survLong_0_u_02 + alpha.slope_02*slope.GK.0_u
      if(event == 2){
        survLong_s_t_0k <- survLong_s_t_0k + alpha.slope_02*slope.GK
      }
    }
  }

  if((c("variability") %in% Objectlsjm$control$sharedtype_01 )|| (c("variability") %in% Objectlsjm$control$sharedtype_02)){
    var.GK <- matrix(rep(omega%*%t(Os),nbQMC),nrow=nbQMC,byrow = T) + b_om%*%t(Ws)
    var.GK.den <- matrix(rep(omega%*%t(Os.den),nbQMC),nrow=nbQMC,byrow = T) + b_om%*%t(Ws.den)
    var.GK.0_u <- matrix(rep(omega%*%t(O_0_u),nbQMC),nrow=nbQMC,byrow = T) + b_om%*%t(W_0_u)
    if(c("variability") %in% Objectlsjm$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + alpha.var_01*exp(var.GK.den)
      survLong_0_u_01 <- survLong_0_u_01 + alpha.var_01*exp(var.GK.0_u)
      if(event == 1){
        survLong_s_t_0k <- survLong_s_t_0k + alpha.var_01*exp(var.GK)
      }
    }
    if(c("variability") %in% Objectlsjm$control$sharedtype_02){
      survLong_0_s_02 <- survLong_0_s_02 + alpha.var_02*exp(var.GK.den)
      survLong_0_u_02 <- survLong_0_u_02 + alpha.var_02*exp(var.GK.0_u)
      if(event == 2){
        survLong_s_t_0k <- survLong_s_t_0k + alpha.var_02*exp(var.GK)
      }
    }
  }

  wk <- wk.1
  if(Objectlsjm$control$hazard_baseline_01 == "Exponential"){
    h_0.GK_0_s_01 <- matrix(wk.1, nrow = 1)
    h_0.GK_0_u_01 <- rep(wk.1, length(wk))
    if(event == 1){
      h_0.GK_s_t_0k <- wk
    }
  }
  else{
    if(Objectlsjm$control$hazard_baseline_01 == "Weibull"){
      h_0.GK_0_s_01 <- shape_01*(st.den**(shape_01-1))*wk
      h_0.GK_0_u_01 <- shape_01*(st_0_u**(shape_01-1))*matrix(rep(wk.1, length(wk)), ncol = 15, byrow = T)
      if(event == 1){
        h_0.GK_s_t_0k <- shape_01*(st.1**(shape_01-1))*wk
      }
    }
    else{
      if(Objectlsjm$control$hazard_baseline_01 == "Gompertz"){
        h_0.GK_0_s_01 <- Gompertz.1_01*exp(Gompertz.2_01*st.den)*wk
        h_0.GK_0_u_01 <- Gompertz.1_01*exp(Gompertz.2_01*st_0_u)*matrix(rep(wk.1, length(wk)), ncol = 15, byrow = T)
        if(event == 1){
          h_0.GK_s_t_0k <- Gompertz.1_01*exp(Gompertz.2_01*st.1)*wk
        }
      }
      else{
        if(Objectlsjm$control$hazard_baseline_01 == "Splines"){
          mat_h0s <- matrix(gamma_01,ncol=1)
          h_0.GK_0_s_01 <- t((wk*exp(Bs.den_01%*%mat_h0s)))
          h_0.GK_0_u_01 <- exp(Bs_0_u_01%*%mat_h0s)*rep(wk, length(wk))

          if(event == 1){
            h_0.GK_s_t_0k <- (wk*exp(Bs_01%*%mat_h0s))
          }
        }
      }
    }
  }

  if(Objectlsjm$control$hazard_baseline_02 == "Exponential"){
    h_0.GK_0_s_02 <- matrix(wk.1, nrow = 1)
    h_0.GK_0_u_02 <- rep(wk.1, length(wk))
    if(event == 2){
      h_0.GK_s_t_0k <- wk
    }
  }
  else{
    if(Objectlsjm$control$hazard_baseline_02 == "Weibull"){
      h_0.GK_0_s_02 <- shape_02*(st.den**(shape_02-1))*wk
      h_0.GK_0_u_02 <- shape_02*(st_0_u**(shape_02-1))*matrix(rep(wk.1, length(wk)), ncol = 15, byrow = T)
      if(event == 2){
        h_0.GK_s_t_0k <- shape_02*(st.1**(shape_02-1))*wk
      }
    }
    else{
      if(Objectlsjm$control$hazard_baseline_02 == "Gompertz"){
        h_0.GK_0_s_02 <- Gompertz.1_02*exp(Gompertz.2_02*st.den)*wk
        h_0.GK_0_u_02 <- Gompertz.1_02*exp(Gompertz.2_02*st_0_u)*matrix(rep(wk.1, length(wk)), ncol = 15, byrow = T)
        if(event == 2){
          h_0.GK_s_t_0k <- Gompertz.1_02*exp(Gompertz.2_02*st.1)*wk
        }
      }
      else{
        if(Objectlsjm$control$hazard_baseline_02 == "Splines"){
          mat_h0s <- matrix(gamma_01,ncol=1)
          h_0.GK_0_s_02 <- t((wk*exp(Bs.den_02%*%mat_h0s)))
          h_0.GK_0_u_02 <- exp(Bs_0_u_02%*%mat_h0s)*rep(wk, length(wk))
          if(event == 2){
            h_0.GK_s_t_0k <- (wk*exp(Bs_02%*%mat_h0s))
          }
        }
      }
    }
  }

  if(length(Z_01)==0){
    pred_surv_01 <- 0
  }
  else{
    pred_surv_01 <- (alpha_01%*%Z_01)[1,1]
  }
  if(length(Z_02)==0){
    pred_surv_02 <- 0
  }
  else{
    pred_surv_02 <- (alpha_02%*%Z_02)[1,1]
  }

  if(event == 1){
    etaBaseline_s_t_0k <- etaBaseline_s_t_0k + pred_surv_01
  }
  else{
    etaBaseline_s_t_0k <- etaBaseline_s_t_0k + pred_surv_02
  }

  etaBaseline_0_s_01 <- etaBaseline_0_s_01 + pred_surv_01
  etaBaseline_0_u_01 <- etaBaseline_0_u_01 + pred_surv_01
  etaBaseline_0_s_02 <- etaBaseline_0_s_02 + pred_surv_02
  etaBaseline_0_u_02 <- etaBaseline_0_u_02 + pred_surv_02





  survLong_s_t_0k <- exp(survLong_s_t_0k)*matrix(rep(t(h_0.GK_s_t_0k), nbQMC), nrow = nbQMC, byrow = TRUE)
  h_0k <- exp(etaBaseline_s_t_0k)*survLong_s_t_0k

  survLong_0_s_01 <- exp(survLong_0_s_01)%*%t(h_0.GK_0_s_01)
  A_0_s_01 <- exp(etaBaseline_0_s_01)*P.den*survLong_0_s_01

  survLong_0_s_02 <- exp(survLong_0_s_02)%*%t(h_0.GK_0_s_02)
  A_0_s_02 <- exp(etaBaseline_0_s_02)*P.den*survLong_0_s_02

  survLong_0_u_01 <- exp(survLong_0_u_01)*matrix(rep(c(t(h_0.GK_0_u_01)),each = nbQMC), nrow = nbQMC, byrow = F)
  survLong_0_u_02 <- exp(survLong_0_u_02)*matrix(rep(c(t(h_0.GK_0_u_02)),each = nbQMC), nrow = nbQMC, byrow = F)


  survLong_red1 <- c()
  survLong_red2 <- c()
  for(nb.col in 1:Objectlsjm$control$nb_pointsGK){
    survLong_red1 <- cbind(survLong_red1, rowSums(survLong_0_u_01[,(Objectlsjm$control$nb_pointsGK*(nb.col-1)+1):(Objectlsjm$control$nb_pointsGK*nb.col)]))
    survLong_red2 <- cbind(survLong_red2, rowSums(survLong_0_u_02[,(Objectlsjm$control$nb_pointsGK*(nb.col-1)+1):(Objectlsjm$control$nb_pointsGK*nb.col)]))
  }
  A1_comp <- 0.5*matrix(rep(st.1, nbQMC), nrow = nbQMC, byrow = T)*exp(etaBaseline_0_u_01)*survLong_red1
  A2_comp <- 0.5*matrix(rep(st.1, nbQMC), nrow = nbQMC, byrow = T)*exp(etaBaseline_0_u_02)*survLong_red2

  Surv.num <- P.1*rowSums(h_0k*exp(-A1_comp - A2_comp))
  Surv.den <- exp(-A_0_s_01-A_0_s_02)


  #browser()
  #f_Y_b_sigma_exp <- exp(f_Y_b_sigma)
  numerateur <- Surv.num*f_Y_b_sigma
  denominateur <- Surv.den*f_Y_b_sigma
  pred.current <- mean(numerateur)/mean(denominateur)

  pred.current


}
