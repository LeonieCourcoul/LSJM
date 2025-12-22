#' @importFrom stats pchisq
#' @export

summary.lsmm_covDep <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsmm_covDep")) stop("use only \"lsmm_covDep\" objects")

  if(x$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==1) cat("    Convergence criteria satisfied")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("We recommend to not interpret the following results and to try to reach convergence for the first step. \n")

    cat("--------------------------------------------------------------------------------------------------------- \n")
    stop("\n")
  }



  cat("Location-scale linear mixed model fitted by maximum likelihood method", "\n")

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
    borne2 <- borne1 + choose(n = x$control$nb.e.a+x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a +x$control$nb.e.a.sigma
    Matcov.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))

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
    borne4 <- borne3 + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a
    borne5 <- borne4 + choose(n = x$control$nb.e.a.sigma, k = 2) + x$control$nb.e.a.sigma
    Matcovb.name.mat <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
    Matcovb.name.mat[lower.tri(Matcovb.name.mat, diag=T)] <- param.names[(borne3+1):borne4]
    Matcovb.name <- unique(gsub("_Location.*", "", unlist(regmatches(Matcovb.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Location", Matcovb.name.mat)))))
    MatcovSig.name.mat <- matrix(rep(0,(x$control$nb.e.a.sigma)**2),nrow=x$control$nb.e.a.sigma,ncol=x$control$nb.e.a.sigma)
    MatcovSig.name.mat[lower.tri(MatcovSig.name.mat, diag=T)] <- param.names[(borne4+1):borne5]
    MatcovSig.name <- unique(gsub("_Scale.*", "", unlist(regmatches(MatcovSig.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Scale", MatcovSig.name.mat)))))


    #Matcovb.name <- param.names[curseur:borne1]
    #MatcovSig.name <- param.names[(borne1+1):borne3]
  }



  cat("\n")
  cat("Longitudinal model:")
  cat("\n")

  cat("      Fixed effects of the location part:")
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
  cat("     Fixed effects of the scale part:")
  var_tab <- matrix(nrow = length(omega), ncol = 4)
  var_tab[,1] <- omega
  var_tab[,2] <- omega.se
  var_tab[,3] <- var_tab[,1]/var_tab[,2]
  var_tab[,4] <- 1 - pchisq(var_tab[,3]**2,1)
  var_tab <- as.data.frame(var_tab)
  rownames(var_tab) <- gsub("_Var", "",omega.name)
  colnames(var_tab) <- c("Coeff", "SE", "Wald", "Pvalue")
  var_tab <- round(var_tab, 4)
  var_tab$Pvalue <- ifelse(var_tab$Pvalue < 0.001, "<0.001", round(var_tab$Pvalue,3))
  cat("\n")
  print(var_tab)

  cat("\n")

  if(x$control$correlated_re){
    cat("     Covariance matrix of the random effects:")
    cat("\n")
    Cov <- MatCov%*%t(MatCov)
    colnames(Cov) <- Matcov.name
    rownames(Cov) <- Matcov.name
    print(Cov)
    #print(MatCov%*%t(MatCov),quote=FALSE,na.print="")
    cat("\n")
  }
  else{
    cat("     Covariance matrix of the location random effects:")
    cat("\n")
    Covb <- MatCovb%*%t(MatCovb)
    colnames(Covb) <- Matcovb.name
    rownames(Covb) <- Matcovb.name
    print(Covb)
    #print(MatCovb%*%t(MatCovb),quote=FALSE,na.print="")
    cat("\n")
    cat("     Covariance matrix of the scale random effects:")
    cat("\n")
    CovSig <- MatCovSig%*%t(MatCovSig)
    colnames(CovSig) <- MatcovSig.name
    rownames(CovSig) <- MatcovSig.name
    print(CovSig)
    #print(MatCovSig%*%t(MatCovSig),quote=FALSE,na.print="")
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
