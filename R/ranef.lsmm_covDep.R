#' @rdname ranef
#' @export
#'
ranef.lsmm_covDep <- function(object,...){

  x <- object
  if(!inherits(x, "lsmm_covDep")) stop("use only \"lsmm_covDep\" objects")
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
  omega <- param[(curseur):(curseur+x$control$nb.omega-1)]
  curseur <- curseur+x$control$nb.omega
  if(x$control$correlated_re){
    C1 <- matrix(rep(0,(x$control$nb.e.a+x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a+x$control$nb.e.a.sigma,ncol=x$control$nb.e.a+x$control$nb.e.a.sigma)
    C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
    Cholesky <- C1
    Cholesky <- as.matrix(Cholesky)
  }
  else{
    borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
    C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    borne3 <- borne1 + choose(n = x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a.sigma
    C3 <- matrix(rep(0,(x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
    C2 <- matrix(0, ncol = x$control$nb.e.a.sigma, nrow = x$control$nb.e.a)

    C4 <- matrix(0, ncol = x$control$nb.e.a, nrow = x$control$nb.e.a.sigma)
    Cholesky <- rbind(cbind(C1,C2),cbind(C4,C3))
    Cholesky <- as.matrix(Cholesky)
  }

  MatCov <- Cholesky%*%t(Cholesky)
  data.long <- x$control$data.long
  random.effects.Predictions <- matrix(NA, nrow = length(unique(data.long$id)), ncol = x$control$nb.e.a+x$control$nb.e.a.sigma+1)
  binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a+x$control$nb.e.a.sigma)

  data.id <- data.long[!duplicated(data.long$id),]
  list.long <- data.manag.long(x$control$formGroup,x$control$formFixed, x$control$formRandom,data.long)
  time.measures <- data.long[,x$control$timeVar]
  X_base <- list.long$X; U_base <- list.long$U; y.new <- list.long$y.new
  offset <- list.long$offset

  list.var <- data.manag.sigma(x$control$formGroup,x$control$formFixedVar, x$control$formRandomVar,data.long)
  O_base <- list.var$X
  O_base <- as.matrix(O_base)
  W_base <- list.var$U
  W_base <- as.matrix(W_base)
  n_cores <- min(x$control$nproc,detectCores() - 1)   # Utiliser tous les cœurs sauf 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Parallélisation de la boucle
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

                                          W_base_i <- W_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                                          W_base_i <- matrix(W_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                                          W_base_i <- unique(W_base_i)

                                          O_base_i <- O_base[offset[id.boucle]:(offset[id.boucle+1]-1),]
                                          O_base_i <- matrix(O_base_i, nrow = offset[id.boucle+1]-offset[id.boucle])
                                          O_base_i <- unique(O_base_i)

                                          # Estimation des effets aléatoires
                                          random.effects_i <- marqLevAlg(binit, fn = re_lsmm_covDep, minimize = FALSE,
                                                                         nb.e.a = x$control$nb.e.a, nb.e.a.sigma = x$control$nb.e.a.sigma,
                                                                         Sigma.re = MatCov, beta = beta, omega = omega,
                                                                         X_base_i = X_base_i, U_base_i = U_base_i, O_base_i = O_base_i,
                                                                         W_base_i = W_base_i, y_i = y_i, nproc = 1,
                                                                         clustertype = x$control$clustertype, maxiter = x$control$maxiter,
                                                                         print.info = FALSE, file = "", blinding = FALSE,
                                                                         epsa = 1e-4, epsb = 1e-4, epsd = 1e-4)

                                          # Réitérer si la convergence n'est pas atteinte
                                          while(random.effects_i$istop != 1) {
                                            binit <- mvtnorm::rmvnorm(1, mean = rep(0, ncol(MatCov)), MatCov)
                                            random.effects_i <- marqLevAlg(binit, fn = re_lsmm_covDep, minimize = FALSE,
                                                                           nb.e.a = x$control$nb.e.a, nb.e.a.sigma = x$control$nb.e.a.sigma,
                                                                           Sigma.re = MatCov, beta = beta, omega = omega,
                                                                           X_base_i = X_base_i, U_base_i = U_base_i, O_base_i = O_base_i,
                                                                           W_base_i = W_base_i, y_i = y_i, nproc = 1,
                                                                           clustertype = x$control$clustertype, maxiter = x$control$maxiter,
                                                                           print.info = FALSE, file = "", blinding = TRUE,
                                                                           epsa = 1e-4, epsb = 1e-4, epsd = 1e-4, multipleTry = 100)
                                            binit <- matrix(0, nrow = 1, ncol = x$control$nb.e.a + x$control$nb.e.a.sigma)
                                          }

                                          # Retourner les résultats pour chaque id.boucle
                                          return(c(data.id$id[id.boucle], random.effects_i$b))
                                        }

  # Fermer le cluster
  stopCluster(cl)
  random.effects.Predictions <- as.data.frame(random.effects.Predictions)
  name_b <- grep("*cov*", rownames(x$table.res), value = TRUE)
  colnames(random.effects.Predictions) <- c("id",unique(unique(gsub("\\*.*", "", gsub("__", "_", name_b)))))


  random.effects.Predictions



}
