#' @rdname predict
#' @export
#'

predict.lsjm_interintraCR <- function(Objectlsmm, which = "RE", Objectranef = NULL, data.long = NULL){

  if(missing(Objectlsmm)) stop("The argument Objectlsmm must be specified")
  if(!inherits((Objectlsmm),"lsjm_interintraCR")) stop("use only \"lsjm_interintraCR\" objects")
 # if(missing(data.long)) stop("The argument data.long must be specified")
  #if(!inherits((data.long),"data.frame")) stop("use only \"data.frame\" objects")
  if(missing(which)) stop("The argument which must be specified")
  if(!inherits((which),"character")) stop("The argument which must be a character object")

  x <- Objectlsmm
  if(x$result_step1$istop != 1|| (!is.null(x$result_step2) && x$result_step2$istop !=1)){
    stop("The model didn't reach convergence.")
  }
  if(is.null(x$result_step2)){
    param <- x$result_step1$b
  }
  else{
    param <- x$result_step2$b
  }

  Cum_risk01 <- NULL
  Cum_risk02 <- NULL

  shape_01 <- 0; shape_02 <- 0;
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0; Gompertz.1_02 <- 0; Gompertz.2_02 <- 0;
  alpha.inter_01 <- 0 ; alpha.inter_02 <- 0;  alpha.intra_01 <- 0;alpha.intra_02 <- 0;
  alpha.current_01 <- 0; alpha.current_02 <- 0;  alpha.slope_01 <- 0; alpha.slope_02 <- 0;
  alpha_01 <- c(0); alpha_02 <- c(0);
  gamma_01 <- c(0); gamma_02 <- c(0);
  beta_slope <- c(0); mu.inter <- 0; sigma.epsilon.inter <-0; mu.intra <- 0;sigma.epsilon.intra <- 0
  alpha_b_01 <- c(0); alpha_b_02 <- c(0)

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
  if("variability inter" %in% x$control$sharedtype_01){
    alpha.inter_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability intra" %in% x$control$sharedtype_01){
    alpha.intra_01 <- param[curseur]
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
    gamma_02 <- param[(curseur):(curseur+x$control$nb.knots.splines[2]-2+1)]
    curseur <- curseur + x$control$nb.knots.splines[2]-2+ 2
  }
  ### Covariables :
  nb.alpha_02 <- x$control$nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("random effects" %in% x$control$sharedtype_02){
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
  if("variability inter" %in% x$control$sharedtype_02){
    alpha.inter_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability intra" %in% x$control$sharedtype_02){
    alpha.intra_02 <- param[curseur]
    curseur <- curseur + 1
  }

  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+ x$control$Objectlsmm$control$nb.beta-1)]
  if( "slope" %in%  x$control$sharedtype_01 || "slope" %in%  x$control$sharedtype_02){
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

  sharedtype <- c("value" %in% x$control$sharedtype_01, "slope" %in% x$control$sharedtype_01, "variability inter" %in% x$control$sharedtype_01, "variability intra" %in% x$control$sharedtype_01,
                  "value" %in% x$control$sharedtype_02, "slope" %in% x$control$sharedtype_02, "variability inter" %in% x$control$sharedtype_02, "variability intra" %in% x$control$sharedtype_02,
                  "random effects" %in% x$control$sharedtype_01, "random effects" %in% x$control$sharedtype_02)
  HB <- list(x$control$hazard_baseline_01, x$control$hazard_baseline_02)
  Weibull <- c(shape_01, shape_02)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01, Gompertz.1_02, Gompertz.2_02)
  alpha_y_slope <- c(alpha.current_01,alpha.current_02,alpha.slope_01,alpha.slope_02)
  alpha_inter_intra <- c(alpha.inter_01,alpha.inter_02,alpha.intra_01,alpha.intra_02)
  alpha_z <- list(alpha_01, alpha_02)
  gamma_z0 <- list(gamma_01, gamma_02)
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
  delta2 <- x$control$deltas[["delta2"]]
  data.long$delta2 <- data.long[all.vars(delta2)][,1]
  knots_01 <- NULL
  knots_02 <- NULL
  random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1)
  if(x$control$Objectlsmm$control$var_inter && x$control$Objectlsmm$control$var_intra){
    random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+2+1)
    binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+2)
  }
  else{
    if(x$control$Objectlsmm$control$var_inter || x$control$Objectlsmm$control$var_intra){
      random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1+1)
      binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+1)
    }
    else{
      random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$Objectlsmm$control$nb.e.a+1)
      binit <- matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)
    }
  }

  st_T = as.matrix(0); X_GK_T = as.matrix(0); U_GK_T = as.matrix(0); Xslope_GK_T = as.matrix(0); Uslope_GK_T = as.matrix(0);
  X_T = as.matrix(0); U_T = as.matrix(0); Xslope_T = as.matrix(0); Uslope_T = as.matrix(0);
  B_T_01 = as.matrix(0); B_T_02 = as.matrix(0);
  Bs_T_01 = as.matrix(0); Bs_T_02 = as.matrix(0); Time_T0 = c(0);
  st_T0 = as.matrix(0); X_GK_T0 = as.matrix(0); U_GK_T0 = as.matrix(0); Xslope_GK_T0 = as.matrix(0); Uslope_GK_T0 = as.matrix(0);
  Bs_T0_01 = as.matrix(0); Bs_T0_02 = as.matrix(0)

  X_GK_T0_i <- as.matrix(0); U_GK_T0_i <- as.matrix(0); X_T_i <- c(0); U_T_i <- c(0);
  X_GK_T_i <- as.matrix(0); U_GK_T_i<- as.matrix(0) ;
  Xslope_T_i <- c(0); Uslope_T_i <- c(0); Time_T0_i <- c(0);
  Xslope_GK_T_i<- as.matrix(0) ; Uslope_GK_T_i<- as.matrix(0);
  Xslope_GK_T0_i<- as.matrix(0); Uslope_GK_T0_i<- as.matrix(0);
  B_T_i_02 <- c(0); B_T_i_01 <- c(0);
  Bs_T_i_01 <- as.matrix(0) ;   Bs_T_i_02 <- as.matrix(0) ;
  Bs_T0_i_01<-as.matrix(0) ;   Bs_T0_i_02 <- as.matrix(0) ;
  st_T0_i <- c(0) ; st_T_i <- c(0); Time_T0_i <- c(0)

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(x$control$Objectlsmm$control$formGroup,x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,data.long)
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  offset <- list.long$offset
  ID.visit <- data.long[all.vars(x$control$Objectlsmm$control$formGroupVisit)][,1];
  offset_ID <- c()
  len_visit <- c(0)
  Ind <- nrow(data.id)
  for(oo in 1:length(Ind)){
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
    B_T_01 <- splineDesign(x$control$knots_01, data.id$Time_T, ord = 4L)
    Bs_T_01 <- splineDesign(x$control$knots_01, c(t(st_T)), ord = 4L)
    if(x$control$left_trunc){
      Bs_T0_01 <- splineDesign(x$control$knots_01, c(t(st_T0)), ord = 4L)
    }
  }
  if(x$control$hazard_baseline_02 == "Splines"){
    Z_02 <- as.matrix(Z_02[,-1])
    B_T_02 <- splineDesign(x$control$knots_02, data.id$Time_T, ord = 4L)
    Bs_T_02 <- splineDesign(x$control$knots_02, c(t(st_T)), ord = 4L)
    if(x$control$left_trunc){
      Bs_T0_02 <- splineDesign(x$control$knots_02, c(t(st_T0)), ord = 4L)
    }
  }

  Cum_risk2 <- c()
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
                              Cum_risk1 = do.call(rbind, lapply(lst, `[[`, 3)),
                              Cum_risk2 = do.call(rbind, lapply(lst, `[[`, 4)))
                       },
                       .multicombine = TRUE,
                       .packages = c("mvtnorm", "marqLevAlg")) %dopar% {

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
                         ID.visit_i <- ID.visit[offset[id.boucle]:(offset[id.boucle+1]-1)]
                         offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
                         len_visit_i <- length(unique(ID.visit_i))

                         random.effects_i <- marqLevAlg(binit, fn = re_lsjm_interintraCR, minimize = FALSE,

                                                        nb.e.a = x$control$Objectlsmm$control$nb.e.a, variability_inter_visit = x$control$Objectlsmm$control$var_inter,
                                                        variability_intra_visit = x$control$Objectlsmm$control$var_intra,
                                                        Sigma.re = MatCov,
                                                        sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                                        alpha_y_slope = alpha_y_slope, alpha_inter_intra = alpha_inter_intra, alpha_b_01 = alpha_b_01, alpha_b_02 = alpha_b_02,
                                                        alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                                        mu.inter = mu.inter , sigma.epsilon.inter = sigma.epsilon.inter, mu.intra = mu.intra,sigma.epsilon.intra = sigma.epsilon.intra,
                                                        delta1_i = delta1_i,delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i,  X_T_i=X_T_i,  U_T_i=U_T_i,
                                                        Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                                        Uslope_GK_T_i=Uslope_GK_T_i,
                                                        X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                                        Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                                        B_T_i_01=B_T_i_01,B_T_i_02=B_T_i_02,
                                                        Bs_T_i_01=Bs_T_i_01, Bs_T_i_02= Bs_T_i_02,
                                                        Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc = x$control$left_trunc,
                                                        len_visit_i = len_visit_i,
                                                        X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, offset_ID_i = offset_ID_i, index_b_slope = x$control$index_b_slope,
                                                        nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                        file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)

                         while(random.effects_i$istop !=1){
                           binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
                           random.effects_i <- marqLevAlg(binit, fn = re_lsjm_interintraCR, minimize = FALSE,

                                                          nb.e.a = x$control$Objectlsmm$control$nb.e.a, variability_inter_visit = x$control$Objectlsmm$control$var_inter,
                                                          variability_intra_visit = x$control$Objectlsmm$control$var_intra,
                                                          Sigma.re = MatCov,
                                                          sharedtype = sharedtype, HB = HB, Gompertz = Gompertz, Weibull = Weibull, nb_pointsGK = x$control$nb_pointsGK,
                                                          alpha_y_slope = alpha_y_slope, alpha_inter_intra = alpha_inter_intra, alpha_b_01 = alpha_b_01, alpha_b_02 = alpha_b_02,
                                                          alpha_z = alpha_z,  gamma_z0 = gamma_z0,  beta = beta,  beta_slope = beta_slope,  wk = wk,
                                                          mu.inter = mu.inter , sigma.epsilon.inter = sigma.epsilon.inter, mu.intra = mu.intra,sigma.epsilon.intra = sigma.epsilon.intra,
                                                          delta1_i = delta1_i,delta2_i = delta2_i, Z_01_i=Z_01_i, Z_02_i=Z_02_i,  X_T_i=X_T_i,  U_T_i=U_T_i,
                                                          Xslope_T_i = Xslope_T_i,  Uslope_T_i = Uslope_T_i,  X_GK_T_i=X_GK_T_i,  U_GK_T_i=U_GK_T_i,  Xslope_GK_T_i=Xslope_GK_T_i,
                                                          Uslope_GK_T_i=Uslope_GK_T_i,
                                                          X_GK_T0_i=X_GK_T0_i, U_GK_T0_i = U_GK_T0_i,Xslope_GK_T0_i=  Xslope_GK_T0_i,  Uslope_GK_T0_i=Uslope_GK_T0_i,
                                                          Time_T_i=Time_T_i,   Time_T0_i= Time_T0_i, st_T_i=st_T_i,    st_T0_i=st_T0_i,
                                                          B_T_i_01=B_T_i_01,B_T_i_02=B_T_i_02,
                                                          Bs_T_i_01=Bs_T_i_01, Bs_T_i_02= Bs_T_i_02,
                                                          Bs_T0_i_01 = Bs_T0_i_01,  Bs_T0_i_02=Bs_T0_i_02,  left_trunc = x$control$left_trunc,
                                                          len_visit_i = len_visit_i,
                                                          X_base_i=X_base_i,  U_base_i=U_base_i,   y_i=y_i, offset_ID_i = offset_ID_i, index_b_slope = x$control$index_b_slope,
                                                          nproc = 1, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                                          file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
                           if(x$control$var_inter && x$control$var_intra){
                             binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+2)
                           }
                           else{
                             if(x$control$var_inter || x$control$var_intra){
                               binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a+1)
                             }
                             else{
                               binit <-matrix(0, nrow = 1, ncol = x$control$Objectlsmm$control$nb.e.a)
                             }
                           }
                         }

                         random_effects_pred <- c(data.id$id[id.boucle], random.effects_i$b)

                         if ('Y' %in% which) {
                           CV <- X_base_i %*% beta + U_base_i %*% random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]
                           time.measures_i <- unique(time.measures[offset[id.boucle]:(offset[id.boucle+1]-1)])
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

                           cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)),
                                            time.measures_i, CV, rep(Varia.inter,length(time.measures_i)), rep(Varia.intra,length(time.measures_i)))
                         } else {
                           cv_pred <- NULL
                         }

                         if("Cum" %in% which){
                           Cum_risk_02i <- c()
                           Cum_risk_01i <- c()
                           for(j in 1:nrow(st_calc.sort.unique)){
                             pred_haz_01 <- 0
                             pred_haz_02 <- 0

                             if("variability inter" %in% x$control$sharedtype_02){
                               pred_haz_02 <- pred_haz_02 + alpha.inter_02*exp(mu.inter+random.effects_i$b[3])
                             }
                             if("variability inter" %in% x$control$sharedtype_01){
                               pred_haz_01 <- pred_haz_01 + alpha.inter_01*exp(mu.inter+random.effects_i$b[3])
                             }


                             if("variability intra" %in% x$control$sharedtype_02){
                               pred_haz_02 <- pred_haz_02 + alpha.intra_02*exp(mu.intra+random.effects_i$b[4])

                             }
                             if("variability intra" %in% x$control$sharedtype_01){
                               pred_haz_01 <- pred_haz_01 + alpha.intra_01*exp(mu.intra+random.effects_i$b[4])

                             }


                             if("value" %in% x$control$sharedtype_01 || "value" %in% x$control$sharedtype_02){
                               list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                                   x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
                               Xs.j <- list.data.GK.current.sigma.sort.unique$Xtime
                               Us.j <- list.data.GK.current.sigma.sort.unique$Utime
                               current.GK <- beta%*%t(Xs.j) + random.effects_i$b[1:(x$control$Objectlsmm$control$nb.e.a)]%*%t(Us.j)
                               if("value" %in% x$control$sharedtype_01){
                                 pred_haz_01 <- pred_haz_01 + alpha.current_01*current.GK
                               }
                               if("value" %in% x$control$sharedtype_02){
                                 pred_haz_02 <- pred_haz_02 + alpha.current_02*current.GK
                               }
                             }

                             if("slope" %in% x$control$sharedtype_01 || "slope" %in% x$control$sharedtype_02){
                               list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
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

                             Z_02_i <- Z_02[id.boucle,]
                             if(length(Z_02_i)==0){
                               pred_surv_02 <- 0
                             }
                             else{
                               pred_surv_02 <- (alpha_02%*%Z_02_i)[1,1]
                             }

                             pred_haz_02 <- pred_haz_02 + pred_surv_02

                             Cum_risk_02i <- c(Cum_risk_02i, P.sort.unique[j]*sum(exp(pred_haz_02)%*%h_0.GK_02))

                           }
                           Cum_risk1 <- Cum_risk_01i
                           Cum_risk2 <- Cum_risk_02i

                         }
                         else{
                           Cum_risk1 <- NULL
                           Cum_risk2 <- NULL
                         }
                         list(random_effects_pred, cv_pred, Cum_risk1, Cum_risk2)

                       }
    random.effects.Predictions <- results$random.effects.Predictions
    cv.Pred <- results$cv.Pred
    Cum_risk01 <- results$Cum_risk1
    Cum_risk02 <- results$Cum_risk2

    # Fermer le cluster
    stopCluster(cl)

  }
  else{

    if('Y' %in% which){
      cv.Pred <- NULL
      for(id.boucle in 1:length(unique(data.long$id))){
      random.effects_i <- as.matrix(Objectranef[id.boucle,-1], nrow = 1)
      X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
      X_base_i <- unique(X_base_i)
      U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
      U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
      U_base_i <- unique(U_i)
      y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]
      ID.visit_i <- ID.visit[offset[id.boucle]:(offset[id.boucle+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      len_visit_i <- length(unique(ID.visit_i))



      CV <- X_base_i %*% beta + U_base_i %*% random.effects_i[1,1:(x$control$Objectlsmm$control$nb.e.a)]
      time.measures_i <- unique(time.measures[offset[id.boucle]:(offset[id.boucle+1]-1)])

      if(x$control$Objectlsmm$control$var_inter && x$control$Objectlsmm$control$var_intra){
        Varia.inter <- exp(mu.inter + random.effects_i[x$control$Objectlsmm$control$nb.e.a+1])
        Varia.intra <- exp(mu.intra + random.effects_i[x$control$Objectlsmm$control$nb.e.a+2])
      }
      else{
        if(x$control$Objectlsmm$control$var_inter){
          Varia.inter <- exp(mu.inter + random.effects_i[x$control$Objectlsmm$control$nb.e.a+1])
          Varia.intra <- sigma.epsilon.intra
        }
        else{
          if(x$control$Objectlsmm$control$var_intra){
            Varia.intra <- exp(mu.intra + random.effects_i[x$control$Objectlsmm$control$nb.e.a+1])
            Varia.inter <- sigma.epsilon.inter
          }
          else{
            Varia.intra <- sigma.epsilon.intra
            Varia.inter <- sigma.epsilon.inter
          }
        }
      }
      cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)),
                       time.measures_i, CV, rep(Varia.inter,length(time.measures_i)), rep(Varia.intra, length(time.measures_i)))
      cv.Pred <- rbind(cv.Pred,cv_pred)
    }
    }
    if('Cum' %in% which){
      for(id.boucle in 1:length(unique(data.long$id))){
        Cum_risk_02i <- c()
        Cum_risk_01i <- c()
        for(j in 1:nrow(st_calc.sort.unique)){
          pred_haz_01 <- 0
          pred_haz_02 <- 0

          if("variability inter" %in% x$control$sharedtype_02){
            pred_haz_02 <- pred_haz_02 + alpha.inter_02*exp(mu.inter+random.effects_i[3])
          }
          if("variability inter" %in% x$control$sharedtype_01){
            pred_haz_01 <- pred_haz_01 + alpha.inter_01*exp(mu.inter+random.effects_i[3])
          }


          if("variability intra" %in% x$control$sharedtype_02){
            pred_haz_02 <- pred_haz_02 + alpha.intra_02*exp(mu.intra+random.effects_i[4])

          }
          if("variability intra" %in% x$control$sharedtype_01){
            pred_haz_01 <- pred_haz_01 + alpha.intra_01*exp(mu.intra+random.effects_i[4])

          }


          if("value" %in% x$control$sharedtype_01 || "value" %in% x$control$sharedtype_02){
            list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                x$control$Objectlsmm$control$formFixed, x$control$Objectlsmm$control$formRandom,x$control$Objectlsmm$control$timeVar)
            Xs.j <- list.data.GK.current.sigma.sort.unique$Xtime
            Us.j <- list.data.GK.current.sigma.sort.unique$Utime
            current.GK <- beta%*%t(Xs.j) + random.effects_i[1:(x$control$Objectlsmm$control$nb.e.a)]%*%t(Us.j)
            if("value" %in% x$control$sharedtype_01){
              pred_haz_01 <- pred_haz_01 + alpha.current_01*current.GK
            }
            if("value" %in% x$control$sharedtype_02){
              pred_haz_02 <- pred_haz_02 + alpha.current_02*current.GK
            }
          }

          if("slope" %in% x$control$sharedtype_01 || "slope" %in% x$control$sharedtype_02){
            list.data.GK.current.sigma.sort.unique <- data.time(data.GaussKronrod.sort.unique$data.id2[(x$control$nb_pointsGK*(id.boucle-1)+1):(x$control$nb_pointsGK*id.boucle),], st_calc.sort.unique[j,],
                                                                x$control$formSlopeFixed, x$control$formSlopeRandom,x$control$Objectlsmm$control$timeVar)
            Xs.slope.j <- list.data.GK.current.sigma.sort.unique$Xtime
            Us.slope.j <- list.data.GK.current.sigma.sort.unique$Utime
            bslope <- random.effects_i[1:(x$control$Objectlsmm$control$nb.e.a)]
            bslope <- bslope[x$control$index_b_slope]
            bslope <- as.matrix(bslope, nrow = 1)
            slope.GK <- beta_slope%*%t(Xs.slope.j) + bslope%*%t(Us.slope.j)
            if("slope" %in% x$control$sharedtype_01){
              pred_haz_01 <- pred_haz_01 + alpha.slope_01*slope.GK
            }
            if("slope" %in% x$control$sharedtype_02){
              pred_haz_02 <- pred_haz_02 + alpha.slope_02*slope.GK
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

          Z_02_i <- Z_02[id.boucle,]
          if(length(Z_02_i)==0){
            pred_surv_02 <- 0
          }
          else{
            pred_surv_02 <- (alpha_02%*%Z_02_i)[1,1]
          }

          pred_haz_02 <- pred_haz_02 + pred_surv_02

          Cum_risk_02i <- c(Cum_risk_02i, P.sort.unique[j]*sum(exp(pred_haz_02)%*%h_0.GK_02))

        }
        Cum_risk01 <- rbind(Cum_risk01, Cum_risk_01i)
        Cum_risk02 <- rbind(Cum_risk02, Cum_risk_02i)







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
    colnames(cv.Pred) <- c("id", "time", "predY", "predSD_inter", "predSD_intra")
  }
  else{
    cv.Pred <- NULL
  }


  resultat <- list(predictRE = random.effects.Predictions, predictY = cv.Pred, predictCum_01 = Cum_risk01, predictCum_02 = Cum_risk02)
  class(resultat) <- c("predict.lsjm_interintraCR")
  resultat








}
