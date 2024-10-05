predyn_boot_lsjm_classicSingle <- function(Objectlsjm, data.long.until.time.s, s, window,  nb.draws){


  if(is.null(Objectlsjm$result_step2)){
    grad <- Objectlsjm$result_step1$grad
    v_mat <- Objectlsjm$result_step1$v
    nbQMC <- Objectlsjm$control$S1
  }
  else{
    grad <- Objectlsjm$result_step2$grad
    v_mat <- Objectlsjm$result_step2$v
    nbQMC <- Objectlsjm$control$S2
  }

  result <- c()
  Hess <- matrix(rep(0,length(grad)**2),nrow=length(grad),ncol=length(grad))
  Hess[upper.tri(Hess, diag=T)] <- v_mat
  Hess2 = Hess + t(Hess)
  diag(Hess2) <- diag(Hess2) - diag(Hess)

  data.long.until.time.s.id <- data.long.until.time.s[1,]

  ### Longitudinal part
  list.long <- data.manag.long(Objectlsjm$control$Objectlsmm$control$formGroup,Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,data.long.until.time.s)
  X_base <- list.long$X
  U <- list.long$U
  y.new.prog <- list.long$y.new.prog



  # Survival part

  data.GaussKronrod.1 <- data.GaussKronrod(data.long.until.time.s.id,a=s,b=s+window,k = Objectlsjm$control$nb_pointsGK)
  P.1 <- data.GaussKronrod.1$P
  st.1 <- data.GaussKronrod.1$st
  wk.1 <- data.GaussKronrod.1$wk
  data.id.1 <- data.GaussKronrod.1$data.id2
  data.GaussKronrod.den <- data.GaussKronrod(data.long.until.time.s.id,a=0,b=s,k = Objectlsjm$control$nb_pointsGK)
  P.den <- data.GaussKronrod.den$P
  st.den <- data.GaussKronrod.den$st
  wk.den <- data.GaussKronrod.den$wk
  data.id.den <- data.GaussKronrod.den$data.id2

  if("value" %in% Objectlsjm$control$sharedtype_01 ){
    list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                       Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs <- list.data.GK.current$Xtime
    Us <- list.data.GK.current$Utime

    list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                           Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs.den <- list.data.GK.current.den$Xtime
    Us.den <- list.data.GK.current.den$Utime

  }

  if("slope" %in% Objectlsjm$control$sharedtype_01){
    list.data.GK.current <-  data.time(data.id.1, c(t(st.1)),
                                       Objectlsjm$control$formSlopeFixed, Objectlsjm$control$formSlopeRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs.slope <- list.data.GK.current$Xtime
    Us.slope <- list.data.GK.current$Utime

    list.data.GK.current.den <-  data.time(data.id.den, c(t(st.den)),
                                           Objectlsjm$control$formSlopeFixed, Objectlsjm$control$formSlopeRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
    Xs.slope.den <- list.data.GK.current.den$Xtime
    Us.slope.den <- list.data.GK.current.den$Utime

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
        Bs_01 <- splines::splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_01, c(t(st.1)), ord = 4L)
        Bs.den_01 <- splines::splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_01, c(t(st.den)), ord = 4L)
      }else{
        stop("This type of base survival function is not implemented.")
      }
    }
  }



  ### Between 0 and u (double integral)
  st_0_u <- c(); X_0_u <- c(); U_0_u <- c(); Xslope_0_u <- c(); Uslope_0_u <- c(); Bs_0_u_01 <- c(); O_0_u <- c(); W_0_u <- c()
  data.id.integrale <- data.long.until.time.s.id
  for(st.integrale in st.1){
    list.GK_0_st.2 <- data.GaussKronrod(data.id.integrale, a = 0, b = st.integrale, k = Objectlsjm$control$nb_pointsGK)
    st.2 <- list.GK_0_st.2$st
    st_0_u <- rbind(st_0_u, st.2)
    if(("value" %in% Objectlsjm$control$sharedtype_01) ){
      list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),Objectlsjm$control$Objectlsmm$control$formFixed, Objectlsjm$control$Objectlsmm$control$formRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
      X_0_st_u <- list.data.GK_0_u$Xtime; U_0_st_u <- list.data.GK_0_u$Utime
      X_0_u <- rbind(X_0_u,X_0_st_u); U_0_u <- rbind(U_0_u,U_0_st_u)
    }
    if(("slope" %in% Objectlsjm$control$sharedtype_01) ){
      list.data.GK_0_u <- data.time(list.GK_0_st.2$data.id2, c(t(st.2)),Objectlsjm$control$formSlopeFixed, Objectlsjm$control$formSlopeRandom,Objectlsjm$control$Objectlsmm$control$timeVar)
      Xslope_0_st_u <- list.data.GK_0_u$Xtime; Uslope_0_st_u <- list.data.GK_0_u$Utime
      Xslope_0_u <- rbind(Xslope_0_u,Xslope_0_st_u); Uslope_0_u <- rbind(Uslope_0_u,Uslope_0_st_u)
    }
    if(Objectlsjm$control$hazard_baseline_01 == "Splines"){
      Bs_0_u_01 <- rbind(Bs_0_u_01,splineDesign(Objectlsjm$control$knots.hazard_baseline.splines_01, c(t(st.2)), ord = 4L))
    }
  }


  for(l in 1:nb.draws){
    if(is.null(Objectlsjm$result_step2)){
      param <- Objectlsjm$result_step1$b
    }
    else{
      param <- Objectlsjm$result_step2$b
    }
    param <- mvtnorm::rmvnorm(1, mean = param_mean, sigma = Hess2)
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
    if("value" %in% Objectlsjm$control$sharedtype_01){
      alpha.current_01 <-  param[curseur]
      curseur <- curseur + 1
    }
    if("slope" %in% Objectlsjm$control$sharedtype_01){
      alpha.slope_01 <- param[curseur]
      curseur <- curseur + 1
    }



    ## Marker
    ### Fixed effects
    beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
    curseur <- curseur+x$control$Objectlsmm$control$nb.beta
    sigma_epsilon <- param[curseur]
    curseur <- curseur +1

    Zq1 <- spacefillr::generate_sobol_owen_set(nbQMC,  Objectlsjm$control$Objectlsmm$control$nb.e.a+Objectlsjm$control$Objectlsmm$control$nb.e.a.sigma)
    Zq <- apply(Zq1, 2, qnorm)

    borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
    C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    Cholesky <- C1
    Cholesky <- as.matrix(Cholesky)
    random.effects <- Zq%*%t(Cholesky)
    b_al <- random.effects[,1:Objectlsjm$control$Objectlsmm$control$nb.e.a]
    b_al <- matrix(b_al, ncol = Objectlsjm$control$Objectlsmm$control$nb.e.a)



    if( "slope" %in%  Objectlsjm$control$sharedtype_01 ){
      beta_slope <- beta[Objectlsjm$control$index_beta_slope]
      b_al_slope <- matrix(b_al[,Objectlsjm$control$index_b_slope], ncol = length(Objectlsjm$control$index_b_slope))
    }


    if(is.null(nrow(X_base))){
      CV <- (beta%*%X_base)[1,1] + b_al%*%U
      f_Y_b_sigma <- dnorm(x=y.new.prog, mean = CV, sd = sigma_epsilon)
    }else{
      f_Y_b_sigma <- rep(1,nbQMC)
      for(k in 1:nrow(X_base)){
        CV <- (beta%*%X_base[k,])[1,1] + b_al%*%U[k,]
        f_Y_b_sigma <- f_Y_b_sigma*dnorm(x = y.new.prog[k], mean = CV, sd = sigma_epsilon)
      }
    }

    # Computation
    etaBaseline_s_t_0k <- 0; survLong_s_t_0k <- 0; etaBaseline_0_s_01 <- 0; survLong_0_s_01 <- 0;
     etaBaseline_0_u_01 <- 0; survLong_0_u_01 <- 0;


    if(c("random effects") %in% Objectlsjm$control$sharedtype_01){
      survLong_0_s_01 <- survLong_0_s_01 + c(b_al%*%alpha_b_01)
      survLong_0_u_01 <- survLong_0_u_01 + c(b_al%*%alpha_b_01)

        survLong_s_t_0k <- survLong_s_t_0k + c(b_al%*%alpha_b_01)

    }

    if((c("value") %in% Objectlsjm$control$sharedtype_01 )){
      current.GK <- matrix(rep(beta%*%t(Xs),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(Us)
      current.GK.den <- matrix(rep(beta%*%t(Xs.den),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(Us.den)
      current.GK.0_u <- matrix(rep(beta%*%t(X_0_u),nbQMC),nrow=nbQMC,byrow = T) + b_al%*%t(U_0_u)
        survLong_0_s_01 <- survLong_0_s_01 + alpha.current_01*current.GK.den
        survLong_0_u_01 <- survLong_0_u_01 + alpha.current_01*current.GK.0_u
          survLong_s_t_0k <- survLong_s_t_0k + alpha.current_01*current.GK

      }


    if((c("slope") %in% Objectlsjm$control$sharedtype_01 )){
      slope.GK <- matrix(rep(beta_slope%*%t(Xs.slope),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Us.slope)
      slope.GK.den <- matrix(rep(beta_slope%*%t(Xs.slope.den),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Us.slope.den)
      slope.GK.0_u <- matrix(rep(beta_slope%*%t(Xslope_0_u),nbQMC),nrow=nbQMC,byrow = T) + b_al_slope%*%t(Uslope_0_u)
        survLong_0_s_01 <- survLong_0_s_01 + alpha.slope_01*slope.GK.den
        survLong_0_u_01 <- survLong_0_u_01 + alpha.slope_01*slope.GK.0_u
          survLong_s_t_0k <- survLong_s_t_0k + alpha.slope_01*slope.GK

    }


    wk <- wk.1
    if(Objectlsjm$control$hazard_baseline_01 == "Exponential"){
      h_0.GK_0_s_01 <- matrix(wk.1, nrow = 1)
      h_0.GK_0_u_01 <- rep(wk.1, length(wk))
        h_0.GK_s_t_0k <- wk

    }
    else{
      if(Objectlsjm$control$hazard_baseline_01 == "Weibull"){
        h_0.GK_0_s_01 <- shape_01*(st.den**(shape_01-1))*wk
        h_0.GK_0_u_01 <- shape_01*(st_0_u**(shape_01-1))*matrix(rep(wk.1, length(wk)), ncol = 15, byrow = T)
          h_0.GK_s_t_0k <- shape_01*(st.1**(shape_01-1))*wk

      }
      else{
        if(Objectlsjm$control$hazard_baseline_01 == "Gompertz"){
          h_0.GK_0_s_01 <- Gompertz.1_01*exp(Gompertz.2_01*st.den)*wk
          h_0.GK_0_u_01 <- Gompertz.1_01*exp(Gompertz.2_01*st_0_u)*matrix(rep(wk.1, length(wk)), ncol = 15, byrow = T)
            h_0.GK_s_t_0k <- Gompertz.1_01*exp(Gompertz.2_01*st.1)*wk

        }
        else{
          if(Objectlsjm$control$hazard_baseline_01 == "Splines"){
            mat_h0s <- matrix(gamma_01,ncol=1)
            h_0.GK_0_s_01 <- t((wk*exp(Bs.den_01%*%mat_h0s)))
            h_0.GK_0_u_01 <- exp(Bs_0_u_01%*%mat_h0s)*rep(wk, length(wk))

              h_0.GK_s_t_0k <- (wk*exp(Bs_01%*%mat_h0s))

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


    etaBaseline_s_t_0k <- etaBaseline_s_t_0k + pred_surv_01


    etaBaseline_0_s_01 <- etaBaseline_0_s_01 + pred_surv_01
    etaBaseline_0_u_01 <- etaBaseline_0_u_01 + pred_surv_01





    survLong_s_t_0k <- exp(survLong_s_t_0k)*matrix(rep(t(h_0.GK_s_t_0k), nbQMC), nrow = nbQMC, byrow = TRUE)
    h_0k <- exp(etaBaseline_s_t_0k)*survLong_s_t_0k

    survLong_0_s_01 <- exp(survLong_0_s_01)%*%t(h_0.GK_0_s_01)
    A_0_s_01 <- exp(etaBaseline_0_s_01)*P.den*survLong_0_s_01



    survLong_0_u_01 <- exp(survLong_0_u_01)*matrix(rep(c(t(h_0.GK_0_u_01)),each = nbQMC), nrow = nbQMC, byrow = F)



    survLong_red1 <- c()
    for(nb.col in 1:Objectlsjm$control$nb_pointsGK){
      survLong_red1 <- cbind(survLong_red1, rowSums(survLong_0_u_01[,(Objectlsjm$control$nb_pointsGK*(nb.col-1)+1):(Objectlsjm$control$nb_pointsGK*nb.col)]))
    }
    A1_comp <- 0.5*matrix(rep(st.1, nbQMC), nrow = nbQMC, byrow = T)*exp(etaBaseline_0_u_01)*survLong_red1

    Surv.num <- P.1*rowSums(h_0k*exp(-A1_comp))
    Surv.den <- exp(-A_0_s_01)


    #browser()
    #f_Y_b_sigma_exp <- exp(f_Y_b_sigma)
    numerateur <- Surv.num*f_Y_b_sigma
    denominateur <- Surv.den*f_Y_b_sigma
    pred.current <- mean(numerateur)/mean(denominateur)
    result <- c(result, pred.current)

  }

  result


}
