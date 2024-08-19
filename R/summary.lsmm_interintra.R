#' @export

summary.lsmm_interintra <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsmm_interintra")) stop("use only \"lsmm_interintra\" objects")

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
  cat(paste("    Likelihood: ", x$result_step2$fn.value),"\n")
  cat(paste("    AIC: ", 2*length(x$result_step2$b) - 2* x$result_step2$fn.value),"\n")

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
  ## Var inter/intra
  if(x$control$var_inter){
    mu.inter <- param[curseur]
    mu.inter.se <- param.se[curseur]
    mu.inter.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.inter <- param[curseur]
    sigma.epsilon.inter.se <- param.se[curseur]
    sigma.epsilon.inter.name <- param.names[curseur]
    curseur <- curseur +1
  }
  if(x$control$var_intra){
    mu.intra <- param[curseur]
    mu.intra.se <- param.se[curseur]
    mu.intra.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.intra <- param[curseur]
    sigma.epsilon.intra.se <- param.se[curseur]
    sigma.epsilon.intra.name <- param.names[curseur]
    curseur <- curseur +1
  }


  ### Cholesky matrix for random effects
  if(x$control$var_inter && x$control$var_intra){
    if(x$control$correlated_re){

      borne1 <- curseur + choose(n = x$control$nb.e.a+2, k = 2) + x$control$nb.e.a +2- 1
      C1 <- matrix(rep(0,(x$control$nb.e.a+2)**2),nrow=x$control$nb.e.a+2,ncol=x$control$nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      MatCov <- C1
      MatCov <- as.matrix(MatCov)
    }
    else{
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      MatCovb <- as.matrix(C1)
      MatCovSig <- as.matrix(C2)
    }
  }
  else{
    if(x$control$var_inter){
      if(x$control$correlated_re){
        borne1 <- curseur + choose(n = x$control$nb.e.a+1, k = 2) + x$control$nb.e.a +1- 1
        C1 <- matrix(rep(0,(x$control$nb.e.a+1)**2),nrow=x$control$nb.e.a+1,ncol=x$control$nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        MatCov <- C1
        MatCov <- as.matrix(MatCov)
      }
      else{
        borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
        C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        MatCovb <- as.matrix(C1)
        MatCovSig <- as.matrix(C2)
      }
    }
    else{
      if(x$control$var_intra){
        if(x$control$correlated_re){
          borne1 <- curseur + choose(n = x$control$nb.e.a+1, k = 2) + x$control$nb.e.a +1- 1
          C1 <- matrix(rep(0,(x$control$nb.e.a+1)**2),nrow=x$control$nb.e.a+1,ncol=x$control$nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          MatCov <- C1
          MatCov <- as.matrix(MatCov)
        }
        else{
          borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
          C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          MatCovb <- as.matrix(C1)
          MatCovSig <- as.matrix(C2)
        }
      }
    }
    if(!x$control$var_inter && !x$control$var_intra){
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a- 1
      C1 <- matrix(rep(0,(x$control$nb.e.a+1)**2),nrow=x$control$nb.e.a+1,ncol=x$control$nb.e.a+1)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      MatCovb <- C1
      MatCovb <- as.matrix(MatCovb)
      MatCov <- C1
      MatCov <- as.matrix(MatCov)
    }

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

  if(!x$control$var_inter){
    cat("     Residual standard error for constant inter-visit variability:")
    var_inter <- matrix(nrow = length(1), ncol = 4)
    var_inter[,1] <- sigma.epsilon.inter
    var_inter[,2] <- sigma.epsilon.inter.se
    var_inter[,3] <- var_inter[,1]/var_inter[,2]
    var_inter[,4] <- 1 - pchisq(var_inter[,3]**2,1)
    var_inter <- as.data.frame(var_inter)
    rownames(var_inter) <- sigma.epsilon.inter.name
    colnames(var_inter) <- c("Coeff", "SE", "Wald", "P-value")
    cat("\n")
    print(var_inter)
    cat("\n")
  }
  if(!x$control$var_intra){
    cat("     Fixed intercept of the intra-visit variability:")
    var_intra <- matrix(nrow = length(1), ncol = 4)
    var_intra[,1] <- sigma.epsilon.intra
    var_intra[,2] <- sigma.epsilon.intra.se
    var_intra[,3] <- var_intra[,1]/var_intra[,2]
    var_intra[,4] <- 1 - pchisq(var_intra[,3]**2,1)
    var_intra <- as.data.frame(var_intra)
    rownames(var_intra) <- sigma.epsilon.intra.name
    colnames(var_intra) <- c("Coeff", "SE", "Wald", "P-value")
    cat("\n")
    print(var_intra)
    cat("\n")
  }


  if(x$control$var_inter){
    cat("     Fixed intercept of the inter-visits variability:")
    var_inter <- matrix(nrow = length(1), ncol = 4)
    var_inter[,1] <- mu.inter
    var_inter[,2] <- mu.inter.se
    var_inter[,3] <- var_inter[,1]/var_inter[,2]
    var_inter[,4] <- 1 - pchisq(var_inter[,3]**2,1)
    var_inter <- as.data.frame(var_inter)
    rownames(var_inter) <- mu.inter.name
    colnames(var_inter) <- c("Coeff", "SE", "Wald", "P-value")
    cat("\n")
    print(var_inter)
    cat("\n")
  }
  if(x$control$var_intra){
    cat("     Fixed intercept of the intra-visit variability:")
    var_intra <- matrix(nrow = length(1), ncol = 4)
    var_intra[,1] <- mu.intra
    var_intra[,2] <- mu.intra.se
    var_intra[,3] <- var_intra[,1]/var_intra[,2]
    var_intra[,4] <- 1 - pchisq(var_intra[,3]**2,1)
    var_intra <- as.data.frame(var_intra)
    rownames(var_intra) <- mu.intra.name
    colnames(var_intra) <- c("Coeff", "SE", "Wald", "P-value")
    cat("\n")
    print(var_intra)
    cat("\n")
  }


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

    if(x$control$var_inter || x$control$var_intra){
      cat("     Covariance matrix of the random effects of the variance:")
      cat("\n")
      print(MatCovSig%*%t(MatCovSig),quote=FALSE,na.print="")
      cat("\n")
    }

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
