#' @rdname predict
#' @export
#'

predict.lsmm_classic <- function(Objectlsmm, which = "RE", Objectranef = NULL, data.long = NULL){

  if(missing(Objectlsmm)) stop("The argument Objectlsmm must be specified")
  if(!inherits((Objectlsmm),"lsmm_classic")) stop("use only \"lsmm_classic\" objects")
  if(missing(data.long)) stop("The argument data.long must be specified")
  if(!inherits((data.long),"data.frame")) stop("use only \"data.frame\" objects")
  if(missing(which)) stop("The argument which must be specified")
  if(!inherits((which),"character")) stop("The argument which must be a character object")
  #if(!is.null(Objectranef)&&!inherits((Objectranef),"ranef.lsmm_classic")) stop("The argument Objectranef must be a \"ranef.lsmm_classic\" object")

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

  #x$control$nproc <- 1
  #Manage parameter
  curseur <- 1
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+x$control$nb.beta-1)]
  curseur <- curseur+x$control$nb.beta
  sigma_epsilon <- param[curseur]
  curseur <- curseur +1
  C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
  Cholesky <- C1
  Cholesky <- as.matrix(Cholesky)

  MatCov <- Cholesky%*%t(Cholesky)
  data.long <- x$control$data.long
  random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$nb.e.a+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a)

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(x$control$formGroup,x$control$formFixed, x$control$formRandom,data.long)
  time.measures <- data.long[,x$control$timeVar]
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  offset <- list.long$offset
  cv.Pred <- c()

  if(is.null(Objectranef) || ('RE' %in% which)){
    n_cores <- min(x$control$nproc,detectCores() - 1)  # Utiliser tous les cœurs sauf 1
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    results <- foreach(id.boucle = 1:length(unique(data.long$id)),
                       .combine = function(...) {
                         lst <- list(...)
                         list(random.effects.Predictions = do.call(rbind, lapply(lst, `[[`, 1)),
                              cv.Pred = do.call(rbind, lapply(lst, `[[`, 2)))
                       },
                       .multicombine = TRUE,
                       .packages = c("mvtnorm", "marqLevAlg")) %dopar% {


                         # Extraire les données spécifiques à l'individu
                         X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle + 1] - 1),]
                         X_base_i <- matrix(X_base_i, nrow = offset[id.boucle + 1] - offset[id.boucle])
                         X_base_i <- unique(X_base_i)

                         U_i <- U_base[offset[id.boucle]:(offset[id.boucle + 1] - 1),]
                         U_i <- matrix(U_i, nrow = offset[id.boucle + 1] - offset[id.boucle])
                         U_base_i <- unique(U_i)

                         y_i <- y.new[offset[id.boucle]:(offset[id.boucle + 1] - 1)]

                         random.effects_i <- marqLevAlg(binit, fn = re_lsmm_classic, minimize = FALSE,
                                                        nb.e.a = x$control$nb.e.a, Sigma.re = MatCov, beta = beta,
                                                        X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i,
                                                        sigma_epsilon = sigma_epsilon, nproc = x$control$nproc,
                                                        clustertype = x$control$clustertype, maxiter = x$control$maxiter,
                                                        print.info = FALSE, file = "", blinding = FALSE,
                                                        epsa = 1e-4, epsb = 1e-4, epsd = 1e-4)

                         # Si l'optimisation n'a pas convergé, essayer à nouveau avec de nouveaux binit
                         while (random.effects_i$istop != 1) {
                           binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
                           random.effects_i <- marqLevAlg(binit, fn = re_lsmm_classic, minimize = FALSE,
                                                          nb.e.a = x$control$nb.e.a, Sigma.re = MatCov, beta = beta,
                                                          X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i,
                                                          sigma_epsilon = sigma_epsilon, nproc = x$control$nproc,
                                                          clustertype = x$control$clustertype, maxiter = x$control$maxiter,
                                                          print.info = FALSE, file = "", blinding = TRUE,
                                                          epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
                           binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a)
                         }
                         random_effects_pred <- c(data.id$id[id.boucle], random.effects_i$b)

                         if ('Y' %in% which) {
                           CV <- X_base_i %*% beta + U_base_i %*% random.effects_i$b[1:(x$control$nb.e.a)]
                           time.measures_i <- time.measures[offset[id.boucle]:(offset[id.boucle + 1] - 1)]
                           cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)), time.measures_i, CV, rep(sigma_epsilon,length(CV)))
                         }
                         else{
                           cv_pred <- NULL
                         }

                         # Retourner les deux résultats
                         list(random_effects_pred, cv_pred)

                       }
    # Extraire les deux tables finales
    random.effects.Predictions <- results$random.effects.Predictions
    cv.Pred <- results$cv.Pred

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
      CV <- X_base_i %*% beta + U_base_i %*% random.effects_i[1,1:(x$control$nb.e.a)]
      cv_pred <- cbind(rep(data.id$id[id.boucle], length(CV)),
                                      time.measures_i, CV, rep(sigma_epsilon, length(CV)))
      cv.Pred <- rbind(cv.Pred,cv_pred)
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



  resultat <- list(predictRE = random.effects.Predictions, predictY = cv.Pred)
  class(resultat) <- c("predict.lsmm_classic")
  resultat



}
