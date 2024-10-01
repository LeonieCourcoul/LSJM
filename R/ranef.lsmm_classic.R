#' ranef : Compute the random effects of the longitudinal submodel
#'
#' @param object A lsmm or lsjm object
#'
#' @name ranef
#' @rdname ranef
#' @export
#'

ranef.lsmm_classic <- function(object,...){

  x <- object
  if(!inherits(x, "lsmm_classic")) stop("use only \"lsmm_classic\" objects")
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
  n_cores <- min(x$control$nproc,detectCores() - 1)  # Utiliser tous les cœurs sauf 1 pour éviter de surcharger
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)



  # Parallélisation avec foreach
  random.effects.Predictions <- foreach(id.boucle = 1:length(unique(data.long$id)),
                                        .combine = 'rbind', .packages = c("mvtnorm", "marqLevAlg")) %dopar% {


                       # Extraire les données spécifiques à l'individu
                       X_base_i <- X_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                       X_base_i <- matrix(X_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                       X_base_i <- unique(X_base_i)

                       U_i <- U_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                       U_i <- matrix(U_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                       U_base_i <- unique(U_i)

                       y_i <- y.new[offset[id.boucle]:(offset[id.boucle+1]-1)]

                       # Estimation des effets aléatoires
                       random.effects_i <- marqLevAlg(binit, fn = re_lsmm_classic, minimize = FALSE,
                                                      nb.e.a = x$control$nb.e.a, Sigma.re = MatCov, beta = beta,
                                                      X_base_i = X_base_i, U_base_i = U_base_i, y_i = y_i,
                                                      sigma_epsilon = sigma_epsilon, nproc = x$control$nproc,
                                                      clustertype = x$control$clustertype, maxiter = x$control$maxiter,
                                                      print.info = FALSE, file = "", blinding = FALSE,
                                                      epsa = 1e-4, epsb = 1e-4, epsd = 1e-4)

                       # Si l'optimisation n'a pas convergé, essayer à nouveau avec de nouveaux binit
                       while(random.effects_i$istop != 1) {
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

                       # Retourner les résultats si nécessaire (par exemple random.effects_i)
                       return(c(data.id$id[id.boucle], random.effects_i$b))  # Si vous voulez collecter des résultats
                     }
  stopCluster(cl)




  random.effects.Predictions <- as.data.frame(random.effects.Predictions)
  name_b <- grep("*cov*", rownames(x$table.res), value = TRUE)
  colnames(random.effects.Predictions) <- c("id",unique(unique(gsub("\\*.*", "", gsub("__", "_", name_b)))))

  random.effects.Predictions









}
