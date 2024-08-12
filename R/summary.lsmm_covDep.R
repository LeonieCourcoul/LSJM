#' @export

summary.lsmm_covDep <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsmm_covDep")) stop("use only \"lsmm_covDep\" objects")

  if(x$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("WARNING", "\n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("We recommend to not interpret the following results and to try to reach convergence for the first step. \n")

    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("\n")
  }



  cat("Linear mixed-effect model for quantitative outcome", "\n")
  cat("with heterogenous variability and fitted by maximum likelihood method", "\n")

  #ajouter le code d'appelle à la fonction
  cat("\n")
  cat("Statistical Model:", "\n")
  cat(paste("    Number of subjects:", x$control$Ind),"\n")
  cat(paste("    Number of observations:", nrow(x$control$data.long)),"\n")

  cat("\n")
  cat("Iteration process:", "\n")

  if(x$info_conv_step2$conv==1) cat("    Convergence criteria satisfied")
  if(x$info_conv_step2$conv==2) cat("    Maximum number of iteration reached without convergence")
  if(x$info_conv_step2$conv==4) {cat("    The program stopped abnormally. No results can be displayed. \n")
  }
  else{
    cat("\n")
    cat(paste("     Number of iterations: ",x$info_conv_step2$niter), "\n")
    cat(paste("     Convergence criteria: parameters =" ,signif(x$info_conv_step2$convcrit[1],3)), "\n")
    cat(paste("                         : likelihood =" ,signif(x$info_conv_step2$convcrit[2],3)), "\n")
    cat(paste("                         : second derivatives =" ,signif(x$info_conv_step2$convcrit[3],3)), "\n")
    cat(paste("     Time of computation :" ,format(x$info_conv_step2$time)))
  }

  cat("\n")
  cat("\n")
  cat("Goodness-of-fit statistics:")
  cat("\n")
  cat(paste("    Likelihood: ", x$result_step1$fn.value),"\n")
  cat(paste("    AIC: ", 2*nrow(x$result_step1$b) - 2* x$result_step1$fn.value),"\n")

  cat("\n")
  cat("Maximum Likelihood Estimates:")
  #Manage parameters
  curseur <- 1
  param <- x$table.res$Estimation
  param.se <- x$table.res$SE
  param.names <- rownames(x$table.res)

  ## Effets fixes trend :
  beta <- param[curseur:(curseur+x$control$nb.beta-1)]
  beta.se <- param.se[curseur:(curseur+x$control$nb.beta-1)]
  beta.name <- param.names[curseur:(curseur+x$control$nb.beta-1)]
  curseur <- curseur+x$control$nb.beta
  ## Effets fixes var :
  omega <- param[curseur:(curseur+x$control$nb.omega-1)]
  omega.se <- param.se[curseur:(curseur+x$control$nb.omega-1)]
  omega.name <- param.names[curseur:(curseur+x$control$nb.omega-1)]
  curseur <- curseur + x$control$nb.omega


  ## Matrice de variance-covariance de l'ensemble des effets aléatoires :
  if(x$control$correlated_re){
    borne1 <- curseur + choose(n = x$control$nb.e.a+x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a +x$control$nb.e.a.sigma- 1
    C1 <- matrix(rep(0,(x$control$nb.e.a+x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a+x$control$nb.e.a.sigma,ncol=x$control$nb.e.a+x$control$nb.e.a.sigma)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    MatCov <- C1
    MatCov <- as.matrix(MatCov)
  }
  else{
    borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
    C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    borne3 <- borne1 + choose(n = x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a.sigma
    C3 <- matrix(rep(0,(x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
    MatCovb <- as.matrix(C1)
    MatCovSig <- as.matrix(C3)
  }



  cat("\n")
  cat("Longitudinal model:")
  cat("\n")

  cat("      Fixed effects:")

  betas_tab <- matrix(nrow = length(beta), ncol = 4)
  betas_tab[,1] <- beta
  betas_tab[,2] <- beta.se
  betas_tab[,3] <- betas_tab[,1]/betas_tab[,2]
  betas_tab[,4] <- 1 - pchisq(betas_tab[,3]**2,1)
  betas_tab <- as.data.frame(betas_tab)
  rownames(betas_tab) <- beta.name
  colnames(betas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  cat("\n")
  print(betas_tab)

  cat("\n")
  cat("     Fixed effects of the linear predictor associated with variability:")
  var_tab <- matrix(nrow = length(omega), ncol = 4)
  var_tab[,1] <- omega
  var_tab[,2] <- omega.se
  var_tab[,3] <- var_tab[,1]/var_tab[,2]
  var_tab[,4] <- 1 - pchisq(var_tab[,3]**2,1)
  var_tab <- as.data.frame(var_tab)
  rownames(var_tab) <- omega.name
  colnames(var_tab) <- c("Coeff", "SE", "Wald", "P-value")
  cat("\n")
  print(var_tab)

  cat("\n")

  if(x$control$correlated_re){
    cat("     Covariance matrix of the random effects:")
    cat("\n")
    print(MatCov%*%t(MatCov),quote=FALSE,na.print="")
    cat("\n")
  }
  else{
    cat("     Covariance matrix of the random effects of the mean:")
    cat("\n")
    print(MatCovb%*%t(MatCovb),quote=FALSE,na.print="")
    cat("\n")
    cat("     Covariance matrix of the random effects of the variance:")
    cat("\n")
    print(MatCovSig%*%t(MatCovSig),quote=FALSE,na.print="")
    cat("\n")
  }

  if(x$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("WARNING", "\n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("We recommend to not interpret the previous results and to try to reach convergence for the first step. \n")

    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("\n")
  }
}
