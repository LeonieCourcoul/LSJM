#' ranef : Compute the random effects of the longitudinal submodel
#'
#' @param object A lsmm or lsjm object
#'
#' @name ranef
#' @rdname ranef
#' @export
#'
ranef.lsjm_classicCR <- function(object,...){

  x <- object
  if(!inherits(x, "lsjm_classicCR")) stop("use only \"lsjm_classicCR\" objects")
  if(x$result_step1$istop != 1|| (!is.null(x$result_step2) && x$result_step2$istop !=1)){
    stop("The model didn't reach convergence.")
  }
  if(is.null(x$result_step2)){
    param <- x$result_step1$b
  }
  else{
    param <- x$result_step2$b
  }


  shape_01 <- 0; shape_02 <- 0;
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0; Gompertz.1_02 <- 0; Gompertz.2_02 <- 0;
  alpha.current_01 <- 0; alpha.current_02 <- 0;  alpha.slope_01 <- 0; alpha.slope_02 <- 0;
  alpha_01 <- c(0); alpha_02 <- c(0);
  gamma_01 <- c(0); gamma_02 <- c(0);
  beta_slope <- c(0);  sigma_epsilon <- 0; alpha_b_01 <- c(0); alpha_b_02 <- c(0)

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
  if("random effects" %in%x$control$sharedtype_01){
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
  if("random effects" %in%x$control$sharedtype_02){
    alpha_b_02 <- param[curseur:(curseur+x$control$Objectlsmm$control$nb.e.a-1)]
    curseur <- curseur + x$control$Objectlsmm$control$nb.e.a
  }
  if("value" %in% x$control$sharedtype_02){
    alpha.current_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_02){
    alpha.slope_02 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
  if( "slope" %in%  x$control$sharedtype_01 || "slope" %in%  x$control$sharedtype_02){
    beta_slope <- beta[x$control$index_beta_slope]
  }
  curseur <- curseur+x$control$Objectlsmm$control$nb.beta
  sigma_epsilon <- param[curseur]
  curseur <- curseur +1

  borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
  C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
  Cholesky <- C1
  Cholesky <- as.matrix(Cholesky)

  MatCov <- Cholesky%*%t(Cholesky)

  sharedtype <- c("value" %in% x$control$sharedtype_01, "slope" %in% x$control$sharedtype_01,
                  "value" %in% x$control$sharedtype_02, "slope" %in% x$control$sharedtype_02,
                  "random effects" %in% x$control$sharedtype_01, "random effects" %in% x$control$sharedtype_02)
  HB <- list(x$control$hazard_baseline_01, x$control$hazard_baseline_02)
  Weibull <- c(shape_01, shape_02)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01, Gompertz.1_02, Gompertz.2_02)
  alpha_y_slope <- c(alpha.current_01,alpha.current_02,alpha.slope_01,alpha.slope_02)
  alpha_z <- list(alpha_01, alpha_02)
  gamma_z0 <- list(gamma_01, gamma_02)
  wk = gaussKronrod()$wk
  rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk))
  sk_GK <- gaussKronrod()$sk
  # Différencier chaque cas, pour chaque cas créer les matrices puis faire une fonction à maximiser (dedans un appel à C++ très proche des fonctions déjà existantes)
  # => il faut appeler marqlevalg individu par individu !
  #1. On associe chaque individu à un cas
  data.long <- x$control$Objectlsmm$control$data.long

  Time_T <- x$control$Time[["Time_T"]]
  data.long$Time_T <- data.long[all.vars(Time_T)][,1]
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
  knots_01 <- NULL
  knots_02 <- NULL


  random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)

  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Time_T0 = c(0); Time_T0_i <- c(0)
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)

  X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
  X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ;
  Xslope_T_i <- c(0); Uslope_T_i <- c(0);
  Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0);
  Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
  B_T_i_02 <- c(0); B_T_i_01 <- c(0);Time_T0_i <- c(0)
  Bs_T_i_01 <- as.matrix(0) ;   Bs_T_i_02 <- as.matrix(0) ;
  Bs_T0_i_01<-as.matrix(0) ;   Bs_T0_i_02 <- as.matrix(0) ;
  st_T0_i <- c(0) ; st_T_i <- c(0)

    data.id <- data.long[!duplicated(data.long$id),]
    list.long <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long)
    X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
    offset <- list.long$offset

    time.measures <- data.long[,x$control$Objectlsmm$control$timeVar]

    list.GK_T <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T, k = x$control$nb_pointsGK)
    st_T <- list.GK_T$st
    if(x$control$left_trunc){
      list.GK_T0 <- data.GaussKronrod(data.id, a = 0, b = data.id$Time_T0, k = x$control$nb_pointsGK)
      st_T0 <- list.GK_T0$st
    }
    if(("value" %in% x$control$sharedtype_01) || ("value" %in% x$control$sharedtype_02) ){
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

    if(("slope" %in% x$control$sharedtype_01) || ("slope" %in% x$control$sharedtype_02) ){
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
    list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_02, data.long)
    Z_02 <- list.surv$Z
    if(x$control$hazard_baseline_02 == "Gompertz"){Z_02 <- as.matrix(Z_02[,-1])}
    if(x$control$hazard_baseline_01 == "Splines"){
      Z_01 <- as.matrix(Z_01[,-1])
      B_T_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, data.id$Time_T, ord = 4L)
      Bs_T_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, c(t(st_T)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, c(t(st_T0)), ord = 4L)
      }
    }
    if(x$control$hazard_baseline_02 == "Splines"){
      Z_02 <- as.matrix(Z_02[,-1])
      B_T_02 <- splineDesign(x$control$knots.hazard_baseline.splines_02, data.id$Time_T, ord = 4L)
      Bs_T_02 <- splineDesign(x$control$knots.hazard_baseline.splines_02, c(t(st_T)), ord = 4L)
      if(x$control$left_trunc){
        Bs_T0_02 <- splineDesign(x$control$knots.hazard_baseline.splines_02, c(t(st_T0)), ord = 4L)
      }
    }


    n_cores <- min(x$control$nproc,detectCores() - 1)  # Utiliser tous les cœurs sauf 1 pour éviter de surcharger
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)



     #Parallélisation avec foreach
    random.effects.Predictions <- foreach(id.boucle = 1:length(unique(data.long$id)),
                                          .combine = 'rbind', .packages = c("mvtnorm", "marqLevAlg")) %dopar% {

    #for(id.boucle in 1:length(unique(data.long$id))){


     # print(id.boucle)


      if("value" %in% x$control$sharedtype_01 || "value" %in% x$control$sharedtype_02){
        X_T_i <- X_T[id.boucle,];U_T_i <- U_T[id.boucle,]
        X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
        if(x$control$left_trunc){
          X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
        }
      }
      if("slope" %in% x$control$sharedtype_01 || "slope" %in%x$control$sharedtype_02){
        Xslope_T_i <- Xslope_T[id.boucle,];Uslope_T_i <- Uslope_T[id.boucle,]
        Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
        if(x$control$left_trunc){
          Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
        }
      }

                                            if("Weibull" %in% c(x$control$hazard_baseline_01,x$control$hazard_baseline_02) ||"Gompertz" %in% c(x$control$hazard_baseline_01,x$control$hazard_baseline_02)){
                                              st_T_i <- st_T[id.boucle,]
                                              if(x$control$left_trunc){
                                                st_T0_i <- st_T0[id.boucle,]
                                              }
                                            }
                                            if("Splines" %in% x$control$hazard_baseline_01){
                                              B_T_i_01 <- B_T_01[id.boucle,]
                                              Bs_T_i_01 <- Bs_T_01[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]
                                              if(x$control$left_trunc){
                                                Bs_T0_i_01 <- Bs_T0_01[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]
                                              }
                                            }
                                            if("Splines" %in% x$control$hazard_baseline_02){
                                              B_T_i_02 <- B_T_02[id.boucle,]
                                              Bs_T_i_02 <- Bs_T_02[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]
                                              if(left_trunc){
                                                Bs_T0_i_02 <- Bs_T0_02[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]
                                              }
                                            }
                                            Z_01_i <- Z_01[id.boucle,]
                                            Z_02_i <- Z_02[id.boucle,]
                                            Time_T_i <- data.id$Time_T[id.boucle]
                                            if(x$control$left_trunc){
                                              Time_T0_i <- data.id$Time_T0[id.boucle]
                                            }
                                            delta2_i <- data.id$delta2[id.boucle]
                                            delta1_i <- data.id$delta1[id.boucle]

                                            X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                                            X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                                            X_base_i <- unique(X_base_i)
                                            U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                                            U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                                            U_base_i <- unique(U_i)
                                            y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]




                                            random.effects_i <- marqLevAlg(binit, fn = re_lsjm_classicCR, minimize = FALSE,

                                                                           nb.e.a = x$control$Objectlsmm$control$nb.e.a, Sigma.re = MatCov,
                                                                           sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                                                           alpha_y_slope = alpha_y_slope, alpha_b_01 = alpha_b_01, alpha_b_02 = alpha_b_02, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                                                           sigma_epsilon = sigma_epsilon,
                                                                           delta1_i = delta1_i,delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i,  X_T_i=X_T_i,  U_T_i=U_T_i,
                                                                           Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                                                           Uslope_GK_T_i=Uslope_GK_T_i,
                                                                           X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                                                           Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                                                           B_T_i_01=B_T_i_01,B_T_i_02=B_T_i_02,
                                                                           Bs_T_i_01=Bs_T_i_01, Bs_T_i_02= Bs_T_i_02,
                                                                           Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc = x$control$left_trunc,
                                                                           X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i,  index_b_slope = x$control$index_b_slope,
                                                                           nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                                           file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

                                            while(random.effects_i$istop !=1){
                                              binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
                                              random.effects_i <- marqLevAlg(binit, fn = re_lsjm_classicCR, minimize = FALSE,

                                                                             nb.e.a = x$control$Objectlsmm$control$nb.e.a, Sigma.re = MatCov,
                                                                             sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                                                             alpha_y_slope = alpha_y_slope, alpha_b_01 = alpha_b_01, alpha_b_02 = alpha_b_02,alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                                                             sigma_epsilon = sigma_epsilon,
                                                                             delta1_i = delta1_i,delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i,  X_T_i=X_T_i,  U_T_i=U_T_i,
                                                                             Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                                                             Uslope_GK_T_i=Uslope_GK_T_i,
                                                                             X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                                                             Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                                                             B_T_i_01=B_T_i_01,B_T_i_02=B_T_i_02,
                                                                             Bs_T_i_01=Bs_T_i_01, Bs_T_i_02= Bs_T_i_02,
                                                                             Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc = x$control$left_trunc,
                                                                             X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i,  index_b_slope = index_b_slope,
                                                                             nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                                             file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
                                              binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)
                                            }
                                            return(c(data.id$id[id.boucle], random.effects_i$b))
                                          }

    stopCluster(cl)




    random.effects.Predictions <- as.data.frame(random.effects.Predictions)
    name_b <- grep("*cov*", rownames(x$table.res), value = TRUE)
    colnames(random.effects.Predictions) <- c("id",unique(unique(gsub("\\*.*", "", gsub("__", "_", name_b)))))

    random.effects.Predictions


}
