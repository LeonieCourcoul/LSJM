ranef.lsjm_interintraSingle <- function(object,...){

  x <- object
  param <- x$result_step2$b
  cv.Pred <- c()
  x$control$nproc <- 1


  shape_01 <- 0;
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0;
  alpha.inter_01 <- 0 ;   alpha.intra_01 <- 0;
  alpha.current_01 <- 0;   alpha.slope_01 <- 0;
  alpha_01 <- c(0);
  gamma_01 <- c(0);
  beta_slope <- c(0); mu.inter <- 0; sigma.epsilon.inter <-0; mu.intra <- 0;sigma.epsilon.intra <- 0

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
  if("current value" %in% x$control$sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% x$control$sharedtype_01){
    alpha.inter_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% x$control$sharedtype_01){
    alpha.intra_01 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
  if( "slope" %in%  x$control$sharedtype_01){
    beta_slope <- beta[x$control$index_beta_slope]
  }
  curseur <- curseur+x$control$Objectlsmm$control$nb.beta
  if(x$control$Objectlsmm$control$var_inter){
    mu.inter <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.inter <- param[curseur]
    curseur <- curseur +1
  }
  if(x$control$Objectlsmm$control$var_intra){
    mu.intra <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.intra <- param[curseur]
    curseur <- curseur +1
  }
  ### Cholesky matrix for random effects
  if(x$control$Objectlsmm$control$var_inter && x$control$Objectlsmm$control$var_intra){
    if(x$control$Objectlsmm$control$correlated_re){

      C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+2)**2),nrow=x$control$Objectlsmm$control$nb.e.a+2,ncol=x$control$Objectlsmm$control$nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }
    else{
      borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      C3 <- matrix(rep(0,2*x$control$Objectlsmm$control$nb.e.a), ncol = x$control$Objectlsmm$control$nb.e.a)
      C4 <- matrix(rep(0,2*x$control$Objectlsmm$control$nb.e.a), nrow = x$control$Objectlsmm$control$nb.e.a)
      Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
      Cholesky <- as.matrix(Cholesky)
    }
  }
  else{
    if(x$control$Objectlsmm$control$var_inter){
      if(x$control$Objectlsmm$control$correlated_re){
        C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+1)**2),nrow=x$control$Objectlsmm$control$nb.e.a+1,ncol=x$control$Objectlsmm$control$nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
        Cholesky <- C1
      }
      else{
        borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
        C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        C3 <- matrix(rep(0,x$control$Objectlsmm$control$nb.e.a), ncol = 1)
        C4 <- matrix(rep(0,x$control$Objectlsmm$control$nb.e.a), nrow = 1)
        Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
        Cholesky <- as.matrix(Cholesky)
      }
    }
    else{
      if(x$control$Objectlsmm$control$var_intra){
        if(x$control$Objectlsmm$control$correlated_re){
          C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+1)**2),nrow=x$control$Objectlsmm$control$nb.e.a+1,ncol=x$control$Objectlsmm$control$nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
          Cholesky <- C1
        }
        else{
          borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
          C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          C3 <- matrix(rep(0,x$control$Objectlsmm$control$nb.e.a), ncol = 1)
          C4 <- matrix(rep(0,x$control$Objectlsmm$control$nb.e.a), nrow = 1)
          Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
          Cholesky <- as.matrix(Cholesky)
        }
      }
    }
    if(!x$control$Objectlsmm$control$var_inter && !x$control$Objectlsmm$control$var_intra){
      C1 <- matrix(rep(0,(length(param)-curseur)**2),nrow=length(param)-curseur,ncol=length(param)-curseur)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }

  }

  MatCov <- Cholesky%*%t(Cholesky)

  sharedtype <- c("current value" %in% x$control$sharedtype_01, "slope" %in% x$control$sharedtype_01, "inter visit variability" %in% x$control$sharedtype_01, "intra visit variability" %in% x$control$sharedtype_01)
  HB <- list(x$control$hazard_baseline_01)
  Weibull <- c(shape_01)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01)
  alpha_y_slope <- c(alpha.current_01,alpha.slope_01)
  alpha_inter_intra <- c(alpha.inter_01 , alpha.intra_01)
  alpha_z <- list(alpha_01)
  gamma_z0 <- list(gamma_01)
  wk = gaussKronrod()$wk
  rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk))
  sk_GK <- gaussKronrod()$sk

  data.long <- x$control$data.long

  random.effects.Predictions <- matrix(NA, nrow = lenght(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)

  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  B_T_01 = as.matrix(0);
  Bs_T_01 = as.matrix(0);   Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0);

  X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
  X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ;
  Xslope_T_i <- c(0); Uslope_T_i <- c(0);
  Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0);
  Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
  B_T_i_01 <- c(0);
  Bs_T_i_01 <- as.matrix(0) ;
  Bs_T0_i_01<-as.matrix(0) ;
  st_T0_i <- c(0) ; st_T_i <- c(0)

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long)
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  offset <- list.long$offset
  ID.visit <- data.long.Case1[all.vars(x$control$Objectlsmm$control$formGroupVisit)][,1];
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

  time.measures <- data.long[,x$control$Objectlsmm$control$timeVar]

  list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = x$control$nb_pointsGK)
  st_T <- list.GK_T$st
  if(x$control$left_trunc){
    list.GK_T0 <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T0, k = x$control$nb_pointsGK)
    st_T0 <- list.GK_T0$st
  }
  if(("current value" %in% x$control$sharedtype_01)  ){
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

  if(("slope" %in% x$control$sharedtype_01)  ){
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

    if("current value" %in% x$control$sharedtype_01 ){
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
    ID.visit_i <- ID.visit[offset[id_boucle]:(offset[id_boucle+1]-1)]
    offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
    len_visit_i <- length(unique(ID.visit_i))



    random.effects_i <- marqLevAlg(binit, fn = re_lsjm_interintraSingle, minimize = FALSE,

                                   nb.e.a = x$control$Objectlsmm$control$nb.e.a, variability_inter_visit = x$control$Objectlsmm$control$var_inter,
                                   variability_intra_visit = x$control$Objectlsmm$control$var_intra,
                                   Sigma.re = MatCov,
                                   sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                   alpha_y_slope = alpha_y_slope, alpha_inter_intra = alpha_inter_intra, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                   mu.inter = mu.inter , sigma.epsilon.inter = sigma.epsilon.inter, mu.intra = mu.intra,sigma.epsilon.intra = sigma.epsilon.intra,
                                   delta1_i = delta1_i,  Z_01_i=Z_01_i,   X_T_i=X_T_i,  U_T_i=U_T_i,
                                   Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                   Uslope_GK_T_i=Uslope_GK_T_i,
                                   X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                   Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                   B_T_i_01=B_T_i_01,
                                   Bs_T_i_01=Bs_T_i_01,
                                   Bs_T0_i_01 = Bs_T0_i_01,   left_trunc = x$control$left_trunc,
                                   len_visit_i = len_visit_i,
                                   X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, offset_ID_i = offset_ID_i, index_b_slope = index_b_slope,
                                   nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                   file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

    while(random.effects_i$istop !=1){
      binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
      random.effects_i <- marqLevAlg(binit, fn = re_lsjm_interintraSingle, minimize = FALSE,

                                     nb.e.a = x$control$Objectlsmm$control$nb.e.a, variability_inter_visit = x$control$Objectlsmm$control$var_inter,
                                     variability_intra_visit = x$control$Objectlsmm$control$var_intra,
                                     Sigma.re = MatCov,
                                     sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                     alpha_y_slope = alpha_y_slope, alpha_inter_intra = alpha_inter_intra, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                     mu.inter = mu.inter , sigma.epsilon.inter = sigma.epsilon.inter, mu.intra = mu.intra,sigma.epsilon.intra = sigma.epsilon.intra,
                                     delta1_i = delta1_i,  Z_01_i=Z_01_i,   X_T_i=X_T_i,  U_T_i=U_T_i,
                                     Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                     Uslope_GK_T_i=Uslope_GK_T_i,
                                     X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                     Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                     B_T_i_01=B_T_i_01,
                                     Bs_T_i_01=Bs_T_i_01,
                                     Bs_T0_i_01 = Bs_T0_i_01,   left_trunc = x$control$left_trunc,
                                     len_visit_i = len_visit_i,
                                     X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, offset_ID_i = offset_ID_i, index_b_slope = index_b_slope,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
      binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)
    }

    random.effects.Predictions[id_boucle,] <- c(data.id$id[id_boucle],random.effects_i$b)
    time.measures_i <- time.measures[offset[id_boucle]:(offset[id_boucle+1]-1)]
    time.measures_i <- unique(time.measures_i)
    CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
    if(x$control$Objectlsmm$control$var_inter && x$control$Objectlsmm$control$var_intra){
      Varia.inter <- exp(mu.inter + random.effects_i$b[x$control$Objectlsmm$control$nb.e.a+1])
      Varia.intra <- exp(mu.intra + random.effects_i$b[x$control$Objectlsmm$control$nb.e.a+2])
    }
    else{
      if(x$control$Objectlsmm$control$var_inter){
        Varia.inter <- exp(mu.inter + random.effects_i$b[x$control$Objectlsmm$control$nb.e.a+1])
        Varia.intra <- sigma.epsilon.intra
      }
      else{
        if(x$control$Objectlsmm$control$var_intra){
          Varia.intra <- exp(mu.intra + random.effects_i$b[x$control$Objectlsmm$control$nb.e.a+1])
          Varia.inter <- sigma.epsilon.inter
        }
        else{
          Varia.intra <- sigma.epsilon.intra
          Varia.inter <- sigma.epsilon.inter
        }
      }
    }
    cv.Pred <- rbind(cv.Pred, cbind(rep(data.id$id[id_boucle], length(CV)),
                                    time.measures_i, CV,  Varia.inter, Varia.intra))
  }



  cv.Pred <- as.data.frame(cv.Pred)
  colnames(cv.Pred) <- c("id", "time", "CV", "Residual_SD_inter", "Residual_SD_intra")
  list(random.effects.Predictions = random.effects.Predictions, cv.Pred = cv.Pred)





}
