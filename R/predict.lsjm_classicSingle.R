#' @rdname predict
#' @export
#'

predict.lsjm_classicSingle <- function(Objectlsjm, which = "RE", Objectranef = NULL, data.long = NULL){

  if(missing(Objectlsjm)) stop("The argument Objectlsjm must be specified")
  if(!inherits((Objectlsjm),"lsjm_classicSingle")) stop("use only \"lsjm_classicSingle\" objects")
  #if(missing(data.long)) stop("The argument data.long must be specified")
  #if(!inherits((data.long),"data.frame")) stop("use only \"data.frame\" objects")
  if(missing(which)) stop("The argument which must be specified")
  if(!inherits((which),"character")) stop("The argument which must be a character object")
  #if(!is.null(Objectranef)&&!inherits((Objectranef),"ranef.lsmm_classic")) stop("The argument Objectranef must be a \"ranef.lsmm_classic\" object")

  x <- Objectlsjm
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
  alpha.current_01 <- 0;  alpha.slope_01 <- 0;
  alpha_01 <- c(0);
  gamma_01 <- c(0);
  beta_slope <- c(0);  sigma_epsilon <- 0; alpha_b_01 <- c(0);

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



  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
  if( "slope" %in%  x$control$sharedtype_01 ){
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
                  "random effects" %in% x$control$sharedtype_01)
  HB <- list(x$control$hazard_baseline_01)
  Weibull <- c(shape_01)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01)
  alpha_y_slope <- c(alpha.current_01,alpha.slope_01)
  alpha_z <- list(alpha_01)
  gamma_z0 <- list(gamma_01)
  wk = gaussKronrod()$wk
  rep_wk = rep(gaussKronrod()$wk, length(gaussKronrod()$wk))
  sk_GK <- gaussKronrod()$sk
  if(is.null(data.long)){
    data.long <- x$control$Objectlsmm$control$data.long
  }

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
  knots_01 <- NULL
  random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)

  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  B_T_01 = as.matrix(0);
  Bs_T_01 = as.matrix(0);  Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0);

  X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
  X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ;
  Xslope_T_i <- c(0); Uslope_T_i <- c(0);
  Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0);
  Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
  B_T_i_01 <- c(0); Time_T0_i <- c(0);
  Bs_T_i_01 <- as.matrix(0) ;
  Bs_T0_i_01<-as.matrix(0) ;
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
  if(("value" %in% x$control$sharedtype_01) ){
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

  list.surv <- data.manag.surv(x$control$Objectlsmm$control$formGroup, x$control$formSurv_01, data.long)
  Z_01 <- list.surv$Z
  if(x$control$hazard_baseline_01 == "Gompertz"){Z_01 <- as.matrix(Z_01[,-1])}
  if(x$control$hazard_baseline_01 == "Splines"){
    Z_01 <- as.matrix(Z_01[,-1])
    B_T_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, data.id$Time_T, ord = 4L)
    Bs_T_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, c(t(st_T)), ord = 4L)
    if(x$control$left_trunc){
      Bs_T0_01 <- splineDesign(x$control$knots.hazard_baseline.splines_01, c(t(st_T0)), ord = 4L)
    }
  }

  Cum_risk1 <- c()
  Time.sort.unique <- unique(sort(data.id$Time_T))
  data.GaussKronrod.sort.unique <- data.GaussKronrod(data.id = data.id, a=0,  b = Time.sort.unique, k = x$control$nb_pointsGK)
  st_calc.sort.unique <- data.GaussKronrod.sort.unique$st
  P.sort.unique <- data.GaussKronrod.sort.unique$P

  if(is.null(Objectranef) || ('RE' %in% which)){

    n_cores <- min(x$control$nproc,detectCores() - 1)  # Utiliser tous les cœurs sauf 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    results <- foreach(id.boucle = 1:length(unique(data.long$id)),
                       .combine = function(...) {
                         lst <- list(...)
                         list(random.effects.Predictions = do.call(rbind, lapply(lst, `[[`, 1)),
                              cv.Pred = do.call(rbind, lapply(lst, `[[`, 2)),
                              Cum_risk1 = do.call(rbind, lapply(lst, `[[`, 3)))
                       },
                       .multicombine = TRUE,
                       .packages = c("mvtnorm", "marqLevAlg")) %dopar% {
#
                        # for(id.boucle in 1:length(unique(data.long$id))){
                        #   print(id.boucle)
                         #  browser()

                         if("value" %in% x$control$sharedtype_01 ){
                           X_T_i <- X_T[id.boucle,];U_T_i <- U_T[id.boucle,]
                           X_GK_T_i <- as.matrix(X_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);U_GK_T_i <- as.matrix(U_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
                           if(x$control$left_trunc){
                             X_GK_T0_i <- as.matrix(X_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);U_GK_T0_i <- as.matrix(U_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
                           }
                         }
                         if("slope" %in% x$control$sharedtype_01){
                           Xslope_T_i <- Xslope_T[id.boucle,];Uslope_T_i <- Uslope_T[id.boucle,]
                           Xslope_GK_T_i <- as.matrix(Xslope_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);Uslope_GK_T_i <- as.matrix(Uslope_GK_T[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
                           if(x$control$left_trunc){
                             Xslope_GK_T0_i <- as.matrix(Xslope_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),]);Uslope_GK_T0_i <- as.matrix(Uslope_GK_T0[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),])
                           }
                         }

                         if("Weibull" %in% c(x$control$hazard_baseline_01) ||"Gompertz" %in% c(x$control$hazard_baseline_01)){
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
                         Z_01_i <- Z_01[id.boucle,]
                         Time_T_i <- data.id$Time_T[id.boucle]
                         if(x$control$left_trunc){
                           Time_T0_i <- data.id$Time_T0[id.boucle]
                         }
                         delta1_i <- data.id$delta1[id.boucle]

                         X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                         X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                         X_base_i <- unique(X_base_i)
                         U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                         U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                         U_base_i <- unique(U_i)
                         y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]


                         random.effects_i <- marqLevAlg(binit, fn = re_lsjm_classicSingle, minimize = FALSE,

                                                        nb.e.a = x$control$Objectlsmm$control$nb.e.a, Sigma.re = MatCov,
                                                        sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                                        alpha_y_slope = alpha_y_slope, alpha_b_01 = alpha_b_01,  alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                                        sigma_epsilon = sigma_epsilon,
                                                        delta1_i = delta1_i, Z_01_i=Z_01_i, X_T_i=X_T_i,  U_T_i=U_T_i,
                                                        Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                                        Uslope_GK_T_i=Uslope_GK_T_i,
                                                        X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                                        Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                                        B_T_i_01=B_T_i_01,
                                                        Bs_T_i_01=Bs_T_i_01,
                                                        Bs_T0_i_01 = Bs_T0_i_01,   left_trunc = x$control$left_trunc,
                                                        X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i,  index_b_slope = x$control$index_b_slope,
                                                        nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                        file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

                         while(random.effects_i$istop !=1){
                           binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
                           random.effects_i <- marqLevAlg(binit, fn = re_lsjm_classicCR, minimize = FALSE,

                                                          nb.e.a = x$control$Objectlsmm$control$nb.e.a, Sigma.re = MatCov,
                                                          sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                                          alpha_y_slope = alpha_y_slope, alpha_b_01 = alpha_b_01, alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                                          sigma_epsilon = sigma_epsilon,
                                                          delta1_i = delta1_i, Z_01_i=Z_01_i,   X_T_i=X_T_i,  U_T_i=U_T_i,
                                                          Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                                          Uslope_GK_T_i=Uslope_GK_T_i,
                                                          X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                                          Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                                          B_T_i_01=B_T_i_01,
                                                          Bs_T_i_01=Bs_T_i_01,
                                                          Bs_T0_i_01 = Bs_T0_i_01,    left_trunc = x$control$left_trunc,
                                                          X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i,  index_b_slope = x$control$index_b_slope,
                                                          nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                          file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
                           binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)
                         }

                         random_effects_pred <- c(data.id$id[id.boucle], random.effects_i$b)

                         if ('Y' %in% which) {
                           CV <- X_base_i %*% beta + U_base_i %*% random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
                           time.measures_i <- time.measures[offset[id.boucle]:(offset[id.boucle + 1] - 1)]
                           cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)), time.measures_i, CV, rep(sigma_epsilon,length(CV)))
                         }
                         else{
                           cv_pred <- NULL
                         }

                         if("Cum" %in% which){
                           Cum_risk_1i <- c()
                           for(j in 1:nrow(st_calc.sort.unique)){
                             pred_haz_01 <- 0

                             if("random effects" %in% x$control$sharedtype_01){
                               pred_haz_01 <- pred_haz_01 + alpha_b_01*random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
                             }



                             if("value" %in% x$control$sharedtype_01 ){
                               list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                                   x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
                               Xs.j <- list.data.GK.current.sigma.sort.unique$Xtime
                               Us.j <- list.data.GK.current.sigma.sort.unique$Utime
                               current.GK <- beta%*%t(Xs.j) + random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]%*%t(Us.j)
                               if("value" %in% x$control$sharedtype_01){
                                 pred_haz_01 <- pred_haz_01 + alpha.current_01*current.GK
                               }
                             }

                             if("slope" %in% x$control$sharedtype_01 ){
                               list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                                   x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
                               Xs.slope.j <- list.data.GK.current.sigma.sort.unique$Xtime
                               Us.slope.j <- list.data.GK.current.sigma.sort.unique$Utime
                               bslope <- random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
                               bslope <- bslope[x$control$index_b_slope]
                               bslope <- matrix(bslope, nrow = 1)
                               slope.GK <- beta_slope%*%t(Xs.slope.j) + bslope%*%t(Us.slope.j)
                               if("slope" %in% x$control$sharedtype_01){
                                 pred_haz_01 <- pred_haz_01 + alpha.slope_01*slope.GK
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

                             Z_01_i <- Z_01[id.boucle,]
                             if(length(Z_01_i)==0){
                               pred_surv_01 <- 0
                             }
                             else{
                               pred_surv_01 <- (alpha_01%*%Z_01_i)[1,1]
                             }

                             pred_haz_01 <- pred_haz_01 + pred_surv_01

                             Cum_risk_1i <- c(Cum_risk_1i, P.sort.unique[j]*sum(exp(pred_haz_01)%*%h_0.GK_01))



                           }

                           Cum_risk1 <- Cum_risk_1i
                         }
                         else{
                           Cum_risk1 <- NULL
                         }

                         list(random_effects_pred, cv_pred, Cum_risk1)

                       }

    # Extraire les quatres tables finales
    random.effects.Predictions <- results$random.effects.Predictions
    cv.Pred <- results$cv.Pred
    Cum_risk01 <- results$Cum_risk1

    # Fermer le cluster
    stopCluster(cl)





  }
  else{
    if('Y' %in% which){
      for(id.boucle in 1:length(unique(data.long$id))){
        X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle + 1] - 1),]
        X_base_i <- matrix(X_base_i, nrow = offset[id.boucle + 1] - offset[id.boucle])
        X_base_i <- unique(X_base_i)

        U_i <- U_base[offset[id.boucle]:(offset[id.boucle + 1] - 1),]
        U_i <- matrix(U_i, nrow = offset[id.boucle + 1] - offset[id.boucle])
        U_base_i <- unique(U_i)

        random.effects_i <- as.matrix(Objectranef[id.boucle,-1], nrow = 1)
        time.measures_i <- time.measures[offset[id.boucle]:(offset[id.boucle+1]-1)]
        CV <- X_base_i %*% beta + U_base_i %*% random.effects_i[1,1:(x$control$Objectlsmm$control$Snb.e.a)]
        cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)),
                         time.measures_i, CV, rep(sigma_epsilon, length(CV)))
        cv.Pred <- rbind(cv.Pred,cv_pred)
      }
    }
    if("Cum" %in% which){
      for(id.boucle in 1:length(unique(data.long$id))){
        Cum_risk_1i <- c()
        for(j in 1:nrow(st_calc.sort.unique)){
          pred_haz_01 <- 0

          if("random effects" %in% x$control$sharedtype_01){
            pred_haz_01 <- pred_haz_01 + alpha_b_01*random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
          }


          if("value" %in% x$control$sharedtype_01 ){
            list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
            Xs.j <- list.data.GK.current.sigma.sort.unique$Xtime
            Us.j <- list.data.GK.current.sigma.sort.unique$Utime
            current.GK <- beta%*%t(Xs.j) + random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]%*%t(Us.j)
            if("value" %in% x$control$sharedtype_01){
              pred_haz_01 <- pred_haz_01 + alpha.current_01*current.GK
            }
          }

          if("slope" %in% x$control$sharedtype_01 ){
            list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
            Xs.slope.j <- list.data.GK.current.sigma.sort.unique$Xtime
            Us.slope.j <- list.data.GK.current.sigma.sort.unique$Utime
            bslope <- random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
            bslope <- bslope[x$control$index_b_slope]
            bslope <- matrix(bslope, nrow = 1)
            slope.GK <- beta_slope%*%t(Xs.slope.j) + bslope%*%t(Us.slope.j)
            if("slope" %in% x$control$sharedtype_01){
              pred_haz_01 <- pred_haz_01 + alpha.slope_01*slope.GK
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

          Z_01_i <- Z_01[id.boucle,]
          if(length(Z_01_i)==0){
            pred_surv_01 <- 0
          }
          else{
            pred_surv_01 <- (alpha_01%*%Z_01_i)[1,1]
          }

          pred_haz_01 <- pred_haz_01 + pred_surv_01

          Cum_risk_1i <- c(Cum_risk_1i, P.sort.unique[j]*sum(exp(pred_haz_01)%*%h_0.GK_01))


        }
        Cum_risk01 <- rbind(Cum_risk01, Cum_risk_1i)

      }



    }
  }

  if('RE' %in% which){
    random.effects.Predictions <- as.data.frame(random.effects.Predictions)
    name_b <- grep("*cov*", rownames(x$table.res), value = TRUE)
    colnames(random.effects.Predictions) <- c("id",unique(unique(gsub("\\*.*", "", gsub("__", "_", name_b)))))
  }
  else{
    random.effects.Predictions <- Objectranef
  }
  if("Y" %in% which){
    cv.Pred <- as.data.frame(cv.Pred)
    colnames(cv.Pred) <- c("id", "time", "predY", "predSD")
  }
  else{
    cv.Pred <- NULL
  }



  resultat <- list(predictRE = random.effects.Predictions, predictY = cv.Pred, predictCum_01 = Cum_risk01)
  class(resultat) <- c("predict.lsjm_classicSingle")
  resultat






}
