#' @rdname surv_marg
#' @importFrom splines splineDesign
#' @importFrom spacefillr generate_sobol_owen_set
#' @importFrom stats model.frame model.matrix
#' @export

surv_marg.lsjm_covDepCR <- function(object, individual, time){

  if(is.null(object$result_step2)){
    param <- object$result_step1$b
    nbQMC <- object$control$S1
  }
  else{
    param <- object$result_step2$b
    nbQMC <- object$control$S2
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


  if( "slope" %in%  object$control$sharedtype_01 || "slope" %in%  object$control$sharedtype_02){
    beta_slope <- beta[object$control$index_beta_slope]
    b_al_slope <- matrix(b_al[,object$control$index_b_slope], ncol = length(object$control$index_b_slope))
  }


  # Survival part
  if(!(object$control$Objectlsmm$control$timeVar %in% colnames(individual))){
    individual[[object$control$Objectlsmm$control$timeVar]] <- NA
  }
  if(!("id" %in% colnames(individual))){
    individual[["id"]] <- 1
  }
  if(!(all.vars(object$control$Objectlsmm$control$formFixed)[1] %in% colnames(individual))){
    individual[[all.vars(object$control$Objectlsmm$control$formFixed)[1]]] <- 100
  }
  data.GaussKronrod.den <-  data.GaussKronrod(individual,a = 0,b=time,k = object$control$nb_pointsGK)
  P.den <- data.GaussKronrod.den$P
  st.den <- data.GaussKronrod.den$st
  wk.den <- data.GaussKronrod.den$wk
  data.id.den <- data.GaussKronrod.den$data.id2
  #
  if("value" %in% object$control$sharedtype_01 ||"value" %in% object$control$sharedtype_02){
    list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                                 object$control$Objectlsmm$control$formFixed, object$control$Objectlsmm$control$formRandom,object$control$Objectlsmm$control$timeVar)
    Xs.den <- list.data.GK.current.den$Xtime
    Us.den <- list.data.GK.current.den$Utime

  }

  if("slope" %in% object$control$sharedtype_01 ||"slope" %in% object$control$sharedtype_02){
    list.data.GK.current.den <-   data.time(data.id.den, c(t(st.den)),
                                                  object$control$formSlopeFixed, object$control$formSlopeRandom,object$control$Objectlsmm$control$timeVar)
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

  if(object$control$hazard_baseline_01 == "Exponential"){
    mfZ <- model.frame(object$control$formSurv_01, data = individual)
    mfZ2 <- model.frame(object$control$formSurv_01, data = object$control$Objectlsmm$control$data.long[!duplicated(object$control$Objectlsmm$control$data.long$id),])
    mfZ <- rbind(mfZ, mfZ2)
    Z_01 <- model.matrix(object$control$formSurv_01, mfZ)[1,]
  }else{
    if(object$control$hazard_baseline_01 == "Weibull" || object$control$hazard_baseline_01 == "Gompertz"){

      mfZ <- model.frame(object$control$formSurv_01, data = individual)
      mfZ2 <- model.frame(object$control$formSurv_01, data = object$control$Objectlsmm$control$data.long[!duplicated(object$control$Objectlsmm$control$data.long$id),])
      mfZ <- rbind(mfZ, mfZ2)
      Z_01 <- model.matrix(object$control$formSurv_01, mfZ)[1,]
    }else{
      if(object$control$hazard_baseline_01 == "Splines"){
        mfZ <- model.frame(object$control$formSurv_01, data = individual)
        mfZ2 <- model.frame(object$control$formSurv_01, data = object$control$Objectlsmm$control$data.long[!duplicated(object$control$Objectlsmm$control$data.long$id),])
        mfZ <- rbind(mfZ, mfZ2)
        Z_01 <- model.matrix(object$control$formSurv_01, mfZ)[1,]
        Z_01 <- Z_01[,-1]
        Bs_01 <- splineDesign(object$control$knots.hazard_baseline.splines_01, c(t(st.1)), ord = 4L)
        Bs.den_01 <- splineDesign(object$control$knots.hazard_baseline.splines_01, c(t(st.den)), ord = 4L)
      }else{
        stop("This type of base survival function is not implemented.")
      }
    }
  }

  if(object$control$hazard_baseline_02 == "Exponential"){
    mfZ <- model.frame(object$control$formSurv_02, data = individual)
    mfZ2 <- model.frame(object$control$formSurv_02, data = object$control$Objectlsmm$control$data.long[!duplicated(object$control$Objectlsmm$control$data.long$id),])
    mfZ <- rbind(mfZ, mfZ2)
    Z_02 <- model.matrix(object$control$formSurv_02, mfZ)[1,]
  }else{
    if(object$control$hazard_baseline_02 == "Weibull" || object$control$hazard_baseline_02 == "Gompertz"){
      mfZ <- model.frame(object$control$formSurv_02, data = individual)
      mfZ2 <- model.frame(object$control$formSurv_02, data = object$control$Objectlsmm$control$data.long[!duplicated(object$control$Objectlsmm$control$data.long$id),])
      mfZ <- rbind(mfZ, mfZ2)
      Z_02 <- model.matrix(object$control$formSurv_02, mfZ)[1,]
    }else{
      if(object$control$hazard_baseline_02 == "Splines"){
        mfZ <- model.frame(object$control$formSurv_02, data = individual)
        mfZ2 <- model.frame(object$control$formSurv_02, data = object$control$Objectlsmm$control$data.long[!duplicated(object$control$Objectlsmm$control$data.long$id),])
        mfZ <- rbind(mfZ, mfZ2)
        Z_02 <- model.matrix(object$control$formSurv_02, mfZ)[1,]
        Z_02 <- Z_02[,-1]
        Bs_02 <- splineDesign(object$control$knots.hazard_baseline.splines_02, c(t(st.1)), ord = 4L)
        Bs.den_02 <- splineDesign(object$control$knots.hazard_baseline.splines_02, c(t(st.den)), ord = 4L)
      }else{
        stop("This type of base survival function is not implemented.")
      }
    }
  }

  # Computation
  etaBaseline_0_s_01 <- 0; survLong_0_s_01 <- 0;
  etaBaseline_0_s_02 <- 0; survLong_0_s_02 <- 0;


  if(c("random effects") %in% object$control$sharedtype_01){
    survLong_0_s_01 <- survLong_0_s_01 + c(b_al%*%alpha_b_01)
  }
  if(c("random effects") %in% object$control$sharedtype_02){
    survLong_0_s_02 <- survLong_0_s_02 + c(b_al%*%alpha_b_02)
  }


  if((c("value") %in% object$control$sharedtype_01 )|| (c("value") %in% object$control$sharedtype_02)){
    current.GK.den <- matrix(rep(beta%*%t(Xs.den),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(Us.den)
    if(c("value") %in% object$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + alpha.current_01*current.GK.den
    }
    if(c("value") %in% object$control$sharedtype_02){
      survLong_0_s_02 <- survLong_0_s_02 + alpha.current_02*current.GK.den
    }
  }

  if((c("slope") %in% object$control$sharedtype_01 )|| (c("slope") %in% object$control$sharedtype_02)){
    slope.GK.den <- matrix(rep(beta_slope%*%t(Xs.slope.den),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Us.slope.den)
    if(c("slope") %in% object$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + alpha.slope_01*slope.GK.den
    }
    if(c("slope") %in% object$control$sharedtype_02){
      survLong_0_s_02 <- survLong_0_s_02 + alpha.slope_02*slope.GK.den
    }
  }


  if((c("variability") %in% Objectlsjm$control$sharedtype_01 )|| (c("variability") %in% Objectlsjm$control$sharedtype_02)){
    var.GK.den <- matrix(rep(omega%*%t(Os.den),nbQMC),nrow=nbQMC,byrow = T) + b_om%*%t(Ws.den)
    if(c("variability") %in% Objectlsjm$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + alpha.var_01*exp(var.GK.den)
    }
    if(c("variability") %in% Objectlsjm$control$sharedtype_02){
      survLong_0_s_02 <- survLong_0_s_02 + alpha.var_02*exp(var.GK.den)
    }
  }



  wk <- wk.den
  if(object$control$hazard_baseline_01 == "Exponential"){
    h_0.GK_0_s_01 <- matrix(wk.1, nrow = 1)

  }
  else{
    if(object$control$hazard_baseline_01 == "Weibull"){
      h_0.GK_0_s_01 <- shape_01*(st.den**(shape_01-1))*wk

    }
    else{
      if(object$control$hazard_baseline_01 == "Gompertz"){
        h_0.GK_0_s_01 <- Gompertz.1_01*exp(Gompertz.2_01*st.den)*wk

      }
      else{
        if(object$control$hazard_baseline_01 == "Splines"){
          mat_h0s <- matrix(gamma_01,ncol=1)
          h_0.GK_0_s_01 <- t((wk*exp(Bs.den_01%*%mat_h0s)))

        }
      }
    }
  }

  if(object$control$hazard_baseline_02 == "Exponential"){
    h_0.GK_0_s_02 <- matrix(wk.1, nrow = 1)

  }
  else{
    if(object$control$hazard_baseline_02 == "Weibull"){
      h_0.GK_0_s_02 <- shape_02*(st.den**(shape_02-1))*wk

    }
    else{
      if(object$control$hazard_baseline_02 == "Gompertz"){
        h_0.GK_0_s_02 <- Gompertz.1_02*exp(Gompertz.2_02*st.den)*wk

      }
      else{
        if(object$control$hazard_baseline_02 == "Splines"){
          mat_h0s <- matrix(gamma_01,ncol=1)
          h_0.GK_0_s_02 <- t((wk*exp(Bs.den_02%*%mat_h0s)))

        }
      }
    }
  }

  #

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



  etaBaseline_0_s_01 <- etaBaseline_0_s_01 + pred_surv_01
  etaBaseline_0_s_02 <- etaBaseline_0_s_02 + pred_surv_02





  survLong_0_s_01 <- exp(survLong_0_s_01)%*%t(h_0.GK_0_s_01)
  A_0_s_01 <- exp(etaBaseline_0_s_01)*P.den*survLong_0_s_01

  survLong_0_s_02 <- exp(survLong_0_s_02)%*%t(h_0.GK_0_s_02)
  A_0_s_02 <- exp(etaBaseline_0_s_02)*P.den*survLong_0_s_02

  Surv.den <- exp(-A_0_s_01-A_0_s_02)


  #
  #f_Y_b_sigma_exp <- exp(f_Y_b_sigma)
  proba <- mean(Surv.den)
  proba



}
