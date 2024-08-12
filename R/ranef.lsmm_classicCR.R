#A REFAIRE

#' ranef : Compute the random effects of the longitudinal submodel
#'
#' @param object A lsmm or lsjm object
#'
#' @rdname ranef
#' @export
#'

ranef.lsmm_classicCR <- function(object,...){

  x <- object
  if(!inherits(x, "lsmm_classicCR")) stop("use only \"lsmm_classicCR\" objects")
  if(x$result_step1$istop != 1|| x$result_step2$istop !=1){
    stop("The model didn't reach convergence.")
  }
  param <- x$result_step2$b
  x$control$nproc <- 1
  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(object$control$hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(object$control$hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]**2
    Gompertz.2_01 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(object$control$hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+object$control$nb.knots.splines[1]+2+1)]
    curseur <- curseur + object$control$nb.knots.splines[1]+2 + 2
  }
  ### Covariables :
  nb.alpha_01 <- object$control$nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+object$control$nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("current value" %in% object$control$sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% object$control$sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Risque 02
  if(object$control$hazard_baseline_02 == "Weibull"){
    shape_02 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(object$control$hazard_baseline_02 == "Gompertz"){
    Gompertz.1_02 <- param[curseur]**2
    Gompertz.2_02 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(object$control$hazard_baseline_02 == "Splines"){
    gamma_02 <- param[(curseur):(curseur+object$control$nb.knots.splines[2]+2+1)]
    curseur <- curseur + object$control$nb.knots.splines[2]+2 + 2
  }
  ### Covariables :
  nb.alpha_02 <- object$control$nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("current value" %in% object$control$sharedtype_02){
    alpha.current_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% object$control$sharedtype_02){
    alpha.slope_02 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+object$Objectlsmm$control$nb.beta-1)]
  curseur <- curseur+object$Objectlsmm$control$nb.beta
  sigma_epsilon <- param[curseur]
  curseur <- curseur +1
  C1 <- matrix(rep(0,(object$Objectlsmm$control$nb.e.a)**2),nrow=object$Objectlsmm$control$nb.e.a,ncol=object$Objectlsmm$control$nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
  MatCov <- C1
  MatCov <- as.matrix(MatCov)

  data.long <- object$Objectlsmm$control$data.long
  random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = object$Objectlsmm$control$nb.e.a+1+choose(n = object$Objectlsmm$control$nb.e.a, k = 2) + object$Objectlsmm$control$nb.e.a)
  binit <- matrix(0, nrow = 1, ncol = object$Objectlsmm$control$nb.e.a)

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(object$Objectlsmm$control$formGroup,object$Objectlsmm$control$formFixed, object$Objectlsmm$control$formRandom,data.long)
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  offset <- list.long$offset
  cv.Pred <- c()
  for(id.boucle in 1:length(unique(data.long$id))){
    X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
    X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
    X_base_i <- unique(X_base_i)
    U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
    U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
    U_base_i <- unique(U_i)
    y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]

    random.effects_i <- marqLevAlg(binit, fn = re_lcmm_classic, minimize = FALSE, nb.e.a = x$control$nb.e.a, Sigma.re= MatCov,beta=beta,
                                   X_base_i = X_base_i, U_base_i = U_base_i, y_i =y_i,sigma_epsilon=sigma_epsilon,
                                   nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                   file = "", blinding = FALSE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4)

    while(random.effects_i$istop !=1){
      binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
      random.effects_i <- marqLevAlg(binit, fn = re_lcmm_classic, minimize = FALSE,
                                     nb.e.a = x$control$nb.e.a, Sigma.re= MatCov,beta=beta,
                                     X_base_i = X_base_i, U_base_i = U_base_i, y_i =y_i,sigma_epsilon=sigma_epsilon,
                                     nproc = x$control$nproc, clustertype = x$control$clustertype, maxiter = x$control$maxiter, print.info = FALSE,
                                     file = "", blinding = TRUE, epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
      binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a)
    }

    #browser()
    CV <- X_base_i%*%beta + U_base_i%*%random.effects_i$b[1:(x$control$nb.e.a)]

    random.effects.Predictions[id.boucle,] <- c(data.id$id[id.boucle] ,random.effects_i$b, random.effects_i$v)
    cv.Pred <- rbind(cv.Pred, cbind(rep(data.id$id[id.boucle], length(CV)),
                                    X_base_i[,2], CV))




  }
  list(random.effects.Predictions = random.effects.Predictions, cv.Pred = cv.Pred)









}
