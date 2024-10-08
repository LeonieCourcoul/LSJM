#' @export

summary.lsmm_interintra <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsmm_interintra")) stop("use only \"lsmm_interintra\" objects")

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
      borne2 <- borne1 + choose(n = x$control$nb.e.a+2, k = 2) + x$control$nb.e.a +2
      Matcov.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))
    }
    else{
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      MatCovb <- as.matrix(C1)
      MatCovSig <- as.matrix(C2)
      borne3 <- borne1 + choose(n = 2, k = 2) + 2

      borne4 <- borne3 + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a
      borne5 <- borne4 + choose(n = 2, k = 2) + 2

      Matcovb.name.mat <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      Matcovb.name.mat[lower.tri(Matcovb.name.mat, diag=T)] <- param.names[(borne3+1):borne4]
      Matcovb.name <- unique(gsub("_Location.*", "", unlist(regmatches(Matcovb.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Location", Matcovb.name.mat)))))

      MatcovSig.name.mat <- matrix(rep(0,(2)**2),nrow=2,ncol=2)
      MatcovSig.name.mat[lower.tri(MatcovSig.name.mat, diag=T)] <- param.names[(borne4+1):borne5]
      MatcovSig.name <- unique(gsub("_Scale.*", "", unlist(regmatches(MatcovSig.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Scale", MatcovSig.name.mat)))))

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
        borne2 <- borne1 + choose(n = x$control$nb.e.a+1, k = 2) + x$control$nb.e.a +1
        Matcov.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))
      }
      else{
        borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
        C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        MatCovb <- as.matrix(C1)
        MatCovSig <- as.matrix(C2)

        borne3 <- borne1 + choose(n = 1, k = 2) + 1

        borne4 <- borne3 + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a
        borne5 <- borne4 + choose(n = 1, k = 2) + 1

        Matcovb.name.mat <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
        Matcovb.name.mat[lower.tri(Matcovb.name.mat, diag=T)] <- param.names[(borne3+1):borne4]
        Matcovb.name <- unique(gsub("_Location.*", "", unlist(regmatches(Matcovb.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Location", Matcovb.name.mat)))))

        MatcovSig.name.mat <- matrix(rep(0,(1)**2),nrow=1,ncol=1)
        MatcovSig.name.mat[lower.tri(MatcovSig.name.mat, diag=T)] <- param.names[(borne4+1):borne5]
        MatcovSig.name <- unique(gsub("_Scale.*", "", unlist(regmatches(MatcovSig.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Scale", MatcovSig.name.mat)))))

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
          borne2 <- borne1 + choose(n = x$control$nb.e.a+1, k = 2) + x$control$nb.e.a +1
          Matcov.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))
        }
        else{
          borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a - 1
          C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          MatCovb <- as.matrix(C1)
          MatCovSig <- as.matrix(C2)

          borne3 <- borne1 + choose(n = 1, k = 2) + 1

          borne4 <- borne3 + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a
          borne5 <- borne4 + choose(n = 1, k = 2) + 1

          Matcovb.name.mat <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
          Matcovb.name.mat[lower.tri(Matcovb.name.mat, diag=T)] <- param.names[(borne3+1):borne4]
          Matcovb.name <- unique(gsub("_Location.*", "", unlist(regmatches(Matcovb.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Location", Matcovb.name.mat)))))

          MatcovSig.name.mat <- matrix(rep(0,(1)**2),nrow=1,ncol=1)
          MatcovSig.name.mat[lower.tri(MatcovSig.name.mat, diag=T)] <- param.names[(borne4+1):borne5]
          MatcovSig.name <- unique(gsub("_Scale.*", "", unlist(regmatches(MatcovSig.name.mat, gregexpr("\\(?[A-Za-z0-9\\.\\^]+\\)?_Scale", MatcovSig.name.mat)))))

        }
      }
    }
    if(!x$control$var_inter && !x$control$var_intra){
      borne1 <- curseur + choose(n = x$control$nb.e.a, k = 2) + x$control$nb.e.a- 1
      C1 <- matrix(rep(0,(x$control$nb.e.a)**2),nrow=x$control$nb.e.a,ncol=x$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      MatCovb <- C1
      MatCovb <- as.matrix(MatCovb)
      MatCov <- C1
      MatCov <- as.matrix(MatCov)
      borne2 <- borne1 + choose(n = x$control$nb.e.a+0, k = 2) + x$control$nb.e.a +0
      Matcovb.name <- unique( unique(gsub("\\*.*", "", gsub("__", "_", param.names[(borne1+1):borne2]))))
    }

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

  if(!x$control$var_inter){
    cat("     Residual standard error for constant inter-visit variability:")
    var_inter <- matrix(nrow = length(1), ncol = 4)
    var_inter[,1] <- sigma.epsilon.inter
    var_inter[,2] <- sigma.epsilon.inter.se
    var_inter[,3] <- var_inter[,1]/var_inter[,2]
    var_inter[,4] <- 1 - pchisq(var_inter[,3]**2,1)
    var_inter <- as.data.frame(var_inter)
    rownames(var_inter) <- sigma.epsilon.inter.name
    colnames(var_inter) <- c("Coeff", "SE", "Wald", "Pvalue")
    var_inter <- round(var_inter, 4)
    var_inter$Pvalue <- ifelse(var_inter$Pvalue < 0.001, "<0.001", round(var_inter$Pvalue,3))
    cat("\n")
    print(var_inter)
    cat("\n")
  }
  if(!x$control$var_intra){
    cat("     Residual standard error for constant intra-visit variability:")
    var_intra <- matrix(nrow = length(1), ncol = 4)
    var_intra[,1] <- sigma.epsilon.intra
    var_intra[,2] <- sigma.epsilon.intra.se
    var_intra[,3] <- var_intra[,1]/var_intra[,2]
    var_intra[,4] <- 1 - pchisq(var_intra[,3]**2,1)
    var_intra <- as.data.frame(var_intra)
    rownames(var_intra) <- sigma.epsilon.intra.name
    colnames(var_intra) <- c("Coeff", "SE", "Wald", "Pvalue")
    var_intra <- round(var_intra, 4)
    var_intra$Pvalue <- ifelse(var_intra$Pvalue < 0.001, "<0.001", round(var_intra$Pvalue,3))
    cat("\n")
    print(var_intra)
    cat("\n")
  }


  if(x$control$var_inter && x$control$var_intra){
    cat("     Fixed intercept of the scale part(inter/intra variabilities):")
    var_inter <- matrix(nrow = length(1), ncol = 4)
    var_inter[,1] <- mu.inter
    var_inter[,2] <- mu.inter.se
    var_inter[,3] <- var_inter[,1]/var_inter[,2]
    var_inter[,4] <- 1 - pchisq(var_inter[,3]**2,1)
    var_inter <- as.data.frame(var_inter)
    var_intra <- matrix(nrow = length(1), ncol = 4)
    var_intra[,1] <- mu.intra
    var_intra[,2] <- mu.intra.se
    var_intra[,3] <- var_intra[,1]/var_intra[,2]
    var_intra[,4] <- 1 - pchisq(var_intra[,3]**2,1)
    var_intra <- as.data.frame(var_intra)
    var_part <- rbind(var_inter, var_intra)
    var_part <- as.data.frame(var_part)
    rownames(var_part) <- c("inter","intra")
    colnames(var_part) <- c("Coeff", "SE", "Wald", "Pvalue")
    var_part <- round(var_part, 4)
    var_part$Pvalue <- ifelse(var_part$Pvalue < 0.001, "<0.001", round(var_part$Pvalue,3))
    cat("\n")
    print(var_part)
    cat("\n")
  }
  else{
    if(x$control$var_inter){
      cat("     Fixed intercept of the scale part (inter variability):")
      var_inter <- matrix(nrow = length(1), ncol = 4)
      var_inter[,1] <- mu.inter
      var_inter[,2] <- mu.inter.se
      var_inter[,3] <- var_inter[,1]/var_inter[,2]
      var_inter[,4] <- 1 - pchisq(var_inter[,3]**2,1)
      var_inter <- as.data.frame(var_inter)
      rownames(var_inter) <- "inter"
      colnames(var_inter) <- c("Coeff", "SE", "Wald", "Pvalue")
      var_inter <- round(var_inter, 4)
      var_inter$Pvalue <- ifelse(var_inter$Pvalue < 0.001, "<0.001", round(var_inter$Pvalue,3))
      cat("\n")
      print(var_inter)
      cat("\n")
    }
    if(x$control$var_intra){
      cat("     Fixed intercept of the scale part (inter variability):")
      var_intra <- matrix(nrow = length(1), ncol = 4)
      var_intra[,1] <- mu.intra
      var_intra[,2] <- mu.intra.se
      var_intra[,3] <- var_intra[,1]/var_intra[,2]
      var_intra[,4] <- 1 - pchisq(var_intra[,3]**2,1)
      var_intra <- as.data.frame(var_intra)
      rownames(var_intra) <- "intra"
      colnames(var_intra) <- c("Coeff", "SE", "Wald", "Pvalue")
      var_intra <- round(var_intra, 4)
      var_intra$Pvalue <- ifelse(var_intra$Pvalue < 0.001, "<0.001", round(var_intra$Pvalue,3))
      cat("\n")
      print(var_intra)
      cat("\n")
    }
  }





  cat("\n")

  if(x$control$correlated_re){
    cat("     Covariance matrix of the random effects:")
    cat("\n")
    Cov <- MatCov%*%t(MatCov)
    colnames(Cov) <- Matcov.name
    rownames(Cov) <- Matcov.name
    print(Cov)
    cat("\n")
  }
  else{
    cat("     Covariance matrix of the random effects of the mean:")
    cat("\n")
    Covb <- MatCovb%*%t(MatCovb)
    colnames(Covb) <- Matcovb.name
    rownames(Covb) <- Matcovb.name
    print(Covb)
    cat("\n")

    if(x$control$var_inter || x$control$var_intra){
      cat("     Covariance matrix of the random effects of the variance:")
      cat("\n")
      CovSig <- MatCovSig%*%t(MatCovSig)
      colnames(CovSig) <- MatcovSig.name
      rownames(CovSig) <- MatcovSig.name
      print(CovSig)
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
