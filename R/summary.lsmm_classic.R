#' @importFrom stats pchisq
#' @export

summary.lsmm_classic <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsmm_classic")) stop("use only \"lsmm_classic\" objects")

  if(x$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("WARNING", "\n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("We recommend to not interpret the following results and to try to reach convergence for the first step. \n")

    cat("--------------------------------------------------------------------------------------------------------- \n")
    stop("\n")
  }



  cat("Linear mixed-effect model fitted by maximum likelihood method", "\n")

  #ajouter le code d'appelle à la fonction
  cat("\n")
  cat("Statistical Model:", "\n")
  cat(paste("    Number of subjects:", x$control$Ind),"\n")
  cat(paste("    Number of observations:", nrow(x$control$data.long)),"\n")

  cat("\n")
  cat("Iteration process:", "\n")

  if(!is.null(x$info_conv_step2)){
    if(x$info_conv_step2$conv==1) cat("    Convergence criteria satisfied")
    if(x$info_conv_step2$conv==2) cat("    Maximum number of iteration reached without convergence")
    if(x$info_conv_step2$conv==4) cat("    The program stopped abnormally. No results can be displayed. \n")
  }
  cat("\n")
  cat(paste("     Number of iterations: "), "\n")
  cat(paste("          Step 1: ",x$result_step1$ni), "\n")
  cat(paste("          Step 2: ",x$info_conv_step2$niter), "\n")
  cat(paste("     Convergence criteria (Step1): parameters =" ,signif(x$info_conv_step1$convcrit[1],3)), "\n")
  cat(paste("                                 : likelihood =" ,signif(x$info_conv_step1$convcrit[2],3)), "\n")
  cat(paste("                                 : second derivatives =" ,signif(x$info_conv_step1$convcrit[3],3)), "\n")
  cat(paste("     Time of computation :" , format(x$time.computation)),  "\n")

  cat("\n")
  cat("Goodness-of-fit statistics:")
  cat("\n")
  cat(paste("    Likelihood: ", round(x$result_step1$fn.value,3)),"\n")
  cat(paste("    AIC: ", round(2*length(x$result_step1$b) - 2* x$result_step1$fn.value,3)),"\n")

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
  sigma_epsilon <- param[curseur]
  sigma_epsilon.se <- param.se[curseur]
  sigma_epsilon.name <- param.names[curseur]
  curseur <- curseur +1

  borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
  C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
  C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
  MatCov <- C1
  MatCov <- as.matrix(MatCov)
  borne2 <- borne1 + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a
  Matcov.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))



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
  rownames(betas_tab) <- gsub("_Y", "",beta.name)
  colnames(betas_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  betas_tab <- round(betas_tab, 4)
  betas_tab$Pvalue <- ifelse(betas_tab$Pvalue < 0.001, "<0.001", round(betas_tab$Pvalue,3))
  cat("\n")
  print(betas_tab)

  cat("\n")
  cat("      Residual variability:")
  sigma_eps_tab <- matrix(nrow = 1, ncol = 4)
  sigma_eps_tab[,1] <- sigma_epsilon
  sigma_eps_tab[,2] <- sigma_epsilon.se
  sigma_eps_tab[,3] <- sigma_eps_tab[,1]/sigma_eps_tab[,2]
  sigma_eps_tab[,4] <- 1 - pchisq(sigma_eps_tab[,3]**2,1)
  sigma_eps_tab <- as.data.frame(sigma_eps_tab)
  rownames(sigma_eps_tab) <- sigma_epsilon.name
  colnames(sigma_eps_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  sigma_eps_tab <- round(sigma_eps_tab, 4)
  sigma_eps_tab$Pvalue <- ifelse(sigma_eps_tab$Pvalue < 0.001, "<0.001", round(sigma_eps_tab$Pvalue,3))
  cat("\n")
  print(sigma_eps_tab)


  cat("\n")

  cat("     Covariance matrix of the random effects:")
  cat("\n")
  Cov <- MatCov%*%t(MatCov)
  colnames(Cov) <- Matcov.name
  rownames(Cov) <- Matcov.name
  print(Cov)
  cat("\n")

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
