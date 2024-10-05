#' ranef : Compute the random effects of the longitudinal submodel
#'
#' @param object A lsmm or lsjm object
#'
#' @name ranef
#' @rdname ranef
#' @export
#'

ranef.lsjm_covDepSingle <- function(object,...){

  x <- object
  if(!inherits(x, "lsjm_covDepSingle")) stop("use only \"lsjm_covDepSingle\" objects")
  if(x$result_step1$istop != 1|| (!is.null(x$result_step2) && x$result_step2$istop !=1)){
    stop("The model didn't reach convergence.")
  }
  if(is.null(x$result_step2)){
    param <- x$result_step1$b
  }
  else{
    param <- x$result_step2$b
  }


  shape_01 <- 0;
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0;
  alpha.current_01 <- 0;   alpha.slope_01 <- 0;
  alpha.var_01 <- 0;
  alpha_01 <- c(0);
  gamma_01 <- c(0);
  beta_slope <- c(0);alpha_b_01 <- c(0);

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
    gamma_01 <- param[(curseur):(curseur+x$control$nb.knots.splines[1]-2+1)]
    curseur <- curseur + x$control$nb.knots.splines[1]-2 + 2
  }
  ### Covariables :
  nb.alpha_01 <- x$control$nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("random effects" %in% x$control$sharedtype_01){
    alpha_b_01 <- param[curseur:(curseur+x$control$Objectlsmm$control$nb.e.a-1)]
    curseur <- curseur + x$control$Objectlsmm$control$nb.e.a
  }
  if("value" %in% x$control$sharedtype_01){
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



  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
  if( "slope" %in%  x$control$sharedtype_01){
    beta_slope <- beta[x$control$index_beta_slope]
  }
  omega <- param[(curseur):(curseur+x$control$Objectlsmm$control$nb.omega-1)]
  curseur <- curseur+x$control$Objectlsmm$control$nb.omega
  if(x$control$correlated_re){
    C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)**2),nrow=x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma,ncol=x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)
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

  sharedtype <- c("value" %in% x$control$sharedtype_01, "slope" %in% x$control$sharedtype_01, "variability" %in% x$control$sharedtype_01, "random effects" %in% x$control$sharedtype_01)
  HB <- list(x$control$hazard_baseline_01)
  Weibull <- c(shape_01)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01)
  alpha_y_slope_var <- c(alpha.current_01,alpha.slope_01, alpha.var_01)

  alpha_z <- list(alpha_01)
  gamma_z0 <- list(gamma_01)
  wk = gaussKronrod()$wk
  rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk))
  sk_GK <- gaussKronrod()$sk

  data.long <- x$control$data.long

  random.effects.Predictions <- matrix(NA, nrow = lenght(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)

  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  O_GK_T = as.matrix(0); W_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0); O_T = as.matrix(0); W_T = as.matrix(0);
  B_T_01 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0); O_GK_T0 = as.matrix(0); W_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0);

  X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
  X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ;
  Xslope_T_i <- c(0); Uslope_T_i <- c(0);
  Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0);
  Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
  O_GK_T0_i <- as.matrix(0); W_GK_T0_i <- as.matrix(0); O_T_i <- c(0); W_T_i <- c(0);
  O_GK_T_i <- as.matrix(0); W_GK_T_i<- as.matrix(0) ;
  B_T_i_01 <- c(0);
  Bs_T_i_01 <- as.matrix(0) ;
  Bs_T0_i_01<-as.matrix(0) ;
  st_T0_i <- c(0) ; st_T_i <- c(0)

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long)
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  offset <- list.long$offset

  list.var <- data.manag.sigma(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,data.long)
  O_base <- list.var$X
  O_base <- as.matrix(O_base)
  W_base <- list.var$U
  W_base <- as.matrix(W_base)

  time.measures <- data.long[,x$control$Objectlsmm$control$timeVar]

  list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = x$control$nb_pointsGK)
  st_T <- list.GK_T$st
  if(x$control$left_trunc){
    list.GK_T0 <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T0, k = x$control$nb_pointsGK)
    st_T0 <- list.GK_T0$st
  }
  if(("value" %in% x$control$sharedtype_01)){
    list.data_T <- data.time(data.id, data.id$Time_T, x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
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

  if(("slope" %in% x$control$sharedtype_01) ){
    list.data_T <- data.time(data.id, data.id$Time_T, x$control$formSlopeFixed, x$control$formSlopeRandom, x$control$Objectlsmm$control$timeVar)
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

  if(("variability" %in% x$control$sharedtype_01) ){
    list.data_T <- data.time(data.id, data.id$Time_T, x$control$Objectlsmm$control$formFixedVar, x$control$Objectlsmm$control$formRandomVar,x$control$Objectlsmm$control$timeVar)
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

  list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long)
  Z_01 <- list.surv$Z
  if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
  if(x$control$hazard_baseline_01 == "Splines"){
    Z_01 <- as.matrix(Z_01[,-1])
    B_T_01 <- splineDesign(x$control$knots_01, data.id$Time_T, ord = 4L)
    Bs_T_01 <- splineDesign(x$control$knots_01, c(t(st_T)), ord = 4L)
    if(x$control$left_trunc){
      Bs_T0_01 <- splineDesign(x$control$knots_01, c(t(st_T0)), ord = 4L)
    }
  }

  for(id_boucle in 1:length(unique(data.long$id))){

    if("value" %in% x$control$sharedtype_01){
      X_T_i <- X_T[id_boucle,];U_T_i <- U_T[id_boucle,]
      X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]);U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),])
      if(x$control$left_trunc){
        X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]);U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),])
      }
    }
    if("slope" %in% x$control$sharedtype_01){
      Xslope_T_i <- Xslope_T[id_boucle,];Uslope_T_i <- Uslope_T[id_boucle,]
      Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]);Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),])
      if(x$control$left_trunc){
        Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]);Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),])
      }
    }
    if("variability" %in% x$control$sharedtype_01){
      O_T_i <- O_T[id_boucle,];W_T_i <- W_T[id_boucle,]
      O_GK_T_i <- as.matrix(O_GK_T[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]);W_GK_T_i <- as.matrix(W_GK_T[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),])
      if(x$control$left_trunc){
        O_GK_T0_i <- as.matrix(O_GK_T0[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]);W_GK_T0_i <- as.matrix(W_GK_T0[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),])
      }
    }

    if("Weibull" %in% c(x$control$hazard_baseline_01) ||"Gompertz" %in% c(x$control$hazard_baseline_01)){
      st_T_i <- st_T[id_boucle,]
      if(x$control$left_trunc){
        st_T0_i <- st_T0[id_boucle,]
      }
    }
    if("Splines" %in% x$control$hazard_baseline_01){
      B_T_i_01 <- B_T_01[id_boucle,]
      Bs_T_i_01 <- Bs_T_01[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]
      if(x$control$left_trunc){
        Bs_T0_i_01 <- Bs_T0_01[(x$control$nb_pointsGK*(id_boucle-1)+1):(x$control$nb_pointsGK*id_boucle),]
      }
    }

    Z_01_i <- Z_01[id_boucle,]
    Time_T_i <- data.id$Time_T[id_boucle]
    if(x$control$left_trunc){
      Time_T0_i <- data.id$Time_T0[id_boucle]
    }
    delta1_i <- data.id$delta1[id_boucle]

    X_base_i <- X_base[offset[id_boucle]:(offset[id_boucle+1]-1),]
    X_base_i <- matrix(X_base_i, nrow = offset[id_boucle+1]-offset[id_boucle])
    X_base_i <- unique(X_base_i)
    U_i <- U_base[offset[id_boucle]:(offset[id_boucle+1]-1),]
    U_i <- matrix(U_i, nrow = offset[id_boucle+1]-offset[id_boucle])
    U_base_i <- unique(U_i)
    y_i <- y.new[offset[id_boucle]:(offset[id_boucle+1]-1)]

    W_base_i <- W_base[offset[id_boucle]:(offset[id_boucle+1]-1),]
    W_base_i <- matrix(W_base_i, nrow = offset[id_boucle+1]-offset[id_boucle])
    W_base_i <- unique(W_base_i)
    O_base_i <- O_base[offset[id_boucle]:(offset[id_boucle+1]-1),]
    O_base_i <- matrix(O_base_i, nrow = offset[id_boucle+1]-offset[id_boucle])
    O_base_i <- unique(O_base_i)



    random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepSingle, minimize = FALSE,

                                   nb.e.a = x$control$Objectlsmm$control$nb.e.a, nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma,
                                   Sigma.re = MatCov,
                                   sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                   alpha_y_slope_var = alpha_y_slope_var, alpha_b_01 = alpha_b_01,  alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,  wk = wk,
                                   delta1_i = delta1_i, Z_01_i=Z_01_i,    X_T_i=X_T_i,  U_T_i=U_T_i,
                                   Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                   X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                   Uslope_GK_T_i=Uslope_GK_T_i, O_GK_T_i=O_GK_T_i,  W_GK_T_i=W_GK_T_i,
                                   X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i, O_GK_T0_i=O_GK_T0_i, W_GK_T0_i = W_GK_T0_i,
                                   Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                   B_T_i_01=B_T_i_01,
                                   Bs_T_i_01=Bs_T_i_01,
                                   Bs_T0_i_01 = Bs_T0_i_01,  left_trunc = x$control$left_trunc,
                                   X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, O_base_i=O_base_i,  W_base_i=W_base_i, index_b_slope = index_b_slope,
                                   nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                   file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

    while(random.effects_i$istop !=1){
      binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
      random.effects_i <- marqLevAlg(binit, fn = re_lsjm_covDepSingle, minimize = FALSE,

                                     nb.e.a = x$control$Objectlsmm$control$nb.e.a, nb.e.a.sigma = x$control$Objectlsmm$control$nb.e.a.sigma,
                                     Sigma.re = MatCov,
                                     sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                     alpha_y_slope_var = alpha_y_slope_var, alpha_b_01 = alpha_b_01,  alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope, omega = omega,  wk = wk,
                                     delta1_i = delta1_i, Z_01_i=Z_01_i,    X_T_i=X_T_i,  U_T_i=U_T_i,
                                     Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i, O_T_i=O_T_i,  W_T_i=W_T_i,
                                     X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                     Uslope_GK_T_i=Uslope_GK_T_i, O_GK_T_i=O_GK_T_i,  W_GK_T_i=W_GK_T_i,
                                     X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i, O_GK_T0_i=O_GK_T0_i, W_GK_T0_i = W_GK_T0_i,
                                     Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                     B_T_i_01=B_T_i_01,
                                     Bs_T_i_01=Bs_T_i_01,
                                     Bs_T0_i_01 = Bs_T0_i_01,  left_trunc = x$control$left_trunc,
                                     X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, O_base_i=O_base_i,  W_base_i=W_base_i, index_b_slope = index_b_slope,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
      binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)
    }

    random.effects.Predictions[id_boucle,] <- c(data.id$id[id_boucle],random.effects_i$b)
    time.measures_i <- time.measures[offset[id_boucle]:(offset[id_boucle+1]-1)]
    time.measures_i <- unique(time.measures_i)
    CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
    Varia <- exp(O_base_i%*%omega + W_base_i%*%random.effects_i$b[(x$control$Objectlsmm$control$nb.e.a+1):(x$control$Objectlsmm$control$nb.e.a+x$control$Objectlsmm$control$nb.e.a.sigma)])
    cv.Pred <- rbind(cv.Pred, cbind(rep(data.id$id[id_boucle], length(CV)),
                                    time.measures_i, CV, Varia))
  }



  cv.Pred <- as.data.frame(cv.Pred)
  colnames(cv.Pred) <- c("id", "time", "CV", "Residual_SD")
  list(random.effects.Predictions = random.effects.Predictions, cv.Pred = cv.Pred)





}
