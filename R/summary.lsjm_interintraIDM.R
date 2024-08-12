#' @export

summary.lsjm_interintraIDM <- function(object,...)
{
  x <- object
  if(!inherits(x, "lsjm_interintraIDM")) stop("use only \"lsjm_interintraIDM\" objects")

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

  cat("Joint location-scale model with within and between visits residual variability for quantitative outcome and an illness-death model", "\n")

  #ajouter le code d'appelle à la fonction
  cat("\n")
  cat("Statistical Model:", "\n")
  cat(paste("    Number of subjects:", x$control$Objectlsmm$control$Ind),"\n")
  cat(paste("    Number of observations:", nrow(x$control$Objectlsmm$control$data.long)),"\n")

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

  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(x$control$hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]
    shape_01.se <- param.se[curseur]
    shape_01.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]
    Gompertz.2_01 <- param[curseur+1]
    Gompertz.1_01.se <- param.se[curseur]
    Gompertz.2_01.se <- param.se[curseur+1]
    Gompertz.1_01.name <- param.names[curseur]
    Gompertz.2_01.name <- param.names[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    gamma_01.se <- param.se[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    gamma_01.name <- param.names[(curseur):(curseur+x$control$nb.knots.splines[1]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[1]+2 + 2
  }
  ### Covariables :
  nb.alpha_01 <- x$control$nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    alpha_01.se <-  param.se[(curseur):(curseur+nb.alpha_01-1)]
    alpha_01.name <-  param.names[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("current value" %in% x$control$sharedtype_01){
    alpha.current_01 <-  param[curseur]
    alpha.current_01.se <-  param.se[curseur]
    alpha.current_01.name <-  param.names[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_01){
    alpha.slope_01 <- param[curseur]
    alpha.slope_01.se <- param.se[curseur]
    alpha.slope_01.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% x$control$sharedtype_01){
    alpha.intervar_01 <- param[curseur]
    alpha.intervar_01.se <- param.se[curseur]
    alpha.intervar_01.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% x$control$sharedtype_01){
    alpha.intravar_01 <- param[curseur]
    alpha.intravar_01.se <- param.se[curseur]
    alpha.intravar_01.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  ## Risque 02
  if(x$control$hazard_baseline_02 == "Weibull"){
    shape_02 <- param[curseur]
    shape_02.se <- param.se[curseur]
    shape_02.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline_02 == "Gompertz"){
    Gompertz.1_02 <- param[curseur]
    Gompertz.2_02 <- param[curseur+1]
    Gompertz.1_02.se <- param.se[curseur]
    Gompertz.2_02.name <- param.names[curseur+1]
    Gompertz.1_02.se <- param.se[curseur]
    Gompertz.2_02.name <- param.names[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_02 == "Splines"){
    gamma_02 <- param[(curseur):(curseur+x$control$nb.knots.splines[2]+2+1)]
    gamma_02.se <- param.se[(curseur):(curseur+x$control$nb.knots.splines[2]+2+1)]
    gamma_02.name <- param.names[(curseur):(curseur+x$control$nb.knots.splines[2]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[2]+2 + 2
  }
  ### Covariables :
  nb.alpha_02 <- x$control$nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    alpha_02.se <-  param.se[(curseur):(curseur+nb.alpha_02-1)]
    alpha_02.name <-  param.names[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("current value" %in% x$control$sharedtype_02){
    alpha.current_02 <- param[curseur]
    alpha.current_02.se <- param.se[curseur]
    alpha.current_02.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_02){
    alpha.slope_02 <- param[curseur]
    alpha.slope_02.se <- param.se[curseur]
    alpha.slope_02.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% x$control$sharedtype_02){
    alpha.intervar_02 <- param[curseur]
    alpha.intervar_02.se <- param.se[curseur]
    alpha.intervar_02.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% x$control$sharedtype_02){
    alpha.intravar_02 <- param[curseur]
    alpha.intravar_02.se <- param.se[curseur]
    alpha.intravar_02.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  ## Risque 12
  if(x$control$hazard_baseline_12 == "Weibull"){
    shape_12 <- param[curseur]
    shape_12.se <- param.se[curseur]
    shape_12.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if(x$control$hazard_baseline_12 == "Gompertz"){
    Gompertz.1_12 <- param[curseur]
    Gompertz.2_12 <- param[curseur+1]
    Gompertz.1_12.se <- param.se[curseur]
    Gompertz.2_12.name <- param.names[curseur+1]
    Gompertz.1_12.se <- param.se[curseur]
    Gompertz.2_12.name <- param.names[curseur+1]
    curseur <- curseur + 2
  }
  if(x$control$hazard_baseline_12 == "Splines"){
    gamma_12 <- param[(curseur):(curseur+x$control$nb.knots.splines[3]+2+1)]
    gamma_12.se <- param.se[(curseur):(curseur+x$control$nb.knots.splines[3]+2+1)]
    gamma_12.name <- param.names[(curseur):(curseur+x$control$nb.knots.splines[3]+2+1)]
    curseur <- curseur + x$control$nb.knots.splines[3]+2 + 2
  }
  ### Covariables :
  nb.alpha_12 <- x$control$nb.alpha[3]
  if(nb.alpha_12 >=1){
    alpha_12 <-  param[(curseur):(curseur+nb.alpha_12-1)]
    alpha_12.se <-  param.se[(curseur):(curseur+nb.alpha_12-1)]
    alpha_12.name <-  param.names[(curseur):(curseur+nb.alpha_12-1)]
    curseur <- curseur+nb.alpha_12
  }
  ### Association
  if("current value" %in% x$control$sharedtype_12){
    alpha.current_12 <- param[curseur]
    alpha.current_12.se <- param.se[curseur]
    alpha.current_12.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% x$control$sharedtype_12){
    alpha.slope_12 <- param[curseur]
    alpha.slope_12.se <- param.se[curseur]
    alpha.slope_12.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% x$control$sharedtype_12){
    alpha.intervar_12 <- param[curseur]
    alpha.intervar_12.se <- param.se[curseur]
    alpha.intervar_12.name <- param.names[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% x$control$sharedtype_12){
    alpha.intravar_12 <- param[curseur]
    alpha.intravar_12.se <- param.se[curseur]
    alpha.intravar_12.name <- param.names[curseur]
    curseur <- curseur + 1
  }


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
  if(x$control$Objectlsmm$control$var_inter && x$control$Objectlsmm$control$var_intra){
    if(x$control$Objectlsmm$control$correlated_re){

      borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a+2, k = 2) + x$control$Objectlsmm$control$nb.e.a +2- 1
      C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+2)**2),nrow=x$control$Objectlsmm$control$nb.e.a+2,ncol=x$control$Objectlsmm$control$nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      MatCov <- C1
      MatCov <- as.matrix(MatCov)
    }
    else{
      borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
      C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      MatCovb <- as.matrix(C1)
      MatCovSig <- as.matrix(C2)
    }
  }
  else{
    if(x$control$Objectlsmm$control$var_inter){
      if(x$control$Objectlsmm$control$correlated_re){
        borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a+1, k = 2) + x$control$nObjectlsmm$control$b.e.a +1- 1
        C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+1)**2),nrow=x$control$Objectlsmm$control$nb.e.a+1,ncol=x$control$Objectlsmm$control$nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        MatCov <- C1
        MatCov <- as.matrix(MatCov)
      }
      else{
        borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
        C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        MatCovb <- as.matrix(C1)
        MatCovSig <- as.matrix(C2)
      }
    }
    else{
      if(x$control$Objectlsmm$control$var_intra){
        if(x$control$Objectlsmm$control$correlated_re){
          borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a+1, k = 2) + x$control$Objectlsmm$control$nb.e.a +1- 1
          C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+1)**2),nrow=x$control$Objectlsmm$control$nb.e.a+1,ncol=x$control$Objectlsmm$control$nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          MatCov <- C1
          MatCov <- as.matrix(MatCov)
        }
        else{
          borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a - 1
          C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a)**2),nrow=x$control$Objectlsmm$control$nb.e.a,ncol=x$control$Objectlsmm$control$nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          MatCovb <- as.matrix(C1)
          MatCovSig <- as.matrix(C2)
        }
      }
    }
    if(!x$control$Objectlsmm$control$var_inter && !x$control$Objectlsmm$control$var_intra){
      borne1 <- curseur + choose(n = x$control$Objectlsmm$control$nb.e.a, k = 2) + x$control$Objectlsmm$control$nb.e.a- 1
      C1 <- matrix(rep(0,(x$control$Objectlsmm$control$nb.e.a+1)**2),nrow=x$control$Objectlsmm$control$nb.e.a+1,ncol=x$control$Objectlsmm$control$nb.e.a+1)
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

  if(!x$control$Objectlsmm$control$var_inter){
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
  if(!x$control$Objectlsmm$control$var_intra){
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


  if(x$control$Objectlsmm$control$var_inter){
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
  if(x$control$Objectlsmm$control$var_intra){
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

  if(x$control$Objectlsmm$control$correlated_re){
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

    if(x$control$Objectlsmm$control$var_inter || x$control$Objectlsmm$control$var_intra){
      cat("     Covariance matrix of the random effects of the variance:")
      cat("\n")
      print(MatCovSig%*%t(MatCovSig),quote=FALSE,na.print="")
      cat("\n")
    }

  }
  cat("Survival models:")
  cat("\n")
  cat("    Transition 0-1:")
  #browser()
  e1_var_tab <- NULL
  e1_share_current_tab <- NULL
  e1_share_slope_tab <- NULL
  e1_alpha_tab <- NULL
  e1_names_tab <- c()

  if(c("current value") %in% x$control$sharedtype_01){
    e1_share_current_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_current_tab[,1] <- alpha.current_01
    e1_share_current_tab[,2] <- alpha.current_01.se
    e1_share_current_tab[,3] <- e1_share_current_tab[,1]/e1_share_current_tab[,2]
    e1_share_current_tab[,4] <- 1 - pchisq(e1_share_current_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.current_01.name)
  }
  if(c("slope") %in% x$control$sharedtype_01){
    e1_share_slope_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_slope_tab[,1] <- c(alpha.slope_01)
    e1_share_slope_tab[,2] <- c(alpha.slope_01.se)
    e1_share_slope_tab[,3] <- e1_share_slope_tab[,1]/e1_share_slope_tab[,2]
    e1_share_slope_tab[,4] <- 1 - pchisq(e1_share_slope_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.slope_01.name)
  }
  if(c("inter visit variability") %in% x$control$sharedtype_01){
    e1_share_intervar_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_intervar_tab[,1] <- c(alpha.intervar_01)
    e1_share_intervar_tab[,2] <- c(alpha.intervar_01.se)
    e1_share_intervar_tab[,3] <- e1_share_intervar_tab[,1]/e1_share_intervar_tab[,2]
    e1_share_intervar_tab[,4] <- 1 - pchisq(e1_share_intervar_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.intervar_01.name)
  }
  if(c("intra visit variability") %in% x$control$sharedtype_01){
    e1_share_intravar_tab <- matrix(nrow = 1, ncol = 4)
    e1_share_intravar_tab[,1] <- c(alpha.intravar_01)
    e1_share_intravar_tab[,2] <- c(alpha.intravar_01.se)
    e1_share_intravar_tab[,3] <- e1_share_intravar_tab[,1]/e1_share_intravar_tab[,2]
    e1_share_intravar_tab[,4] <- 1 - pchisq(e1_share_intravar_tab[,3]**2,1)
    e1_names_tab <- c(e1_names_tab, alpha.intravar_01.name)
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    if(x$control$nb.alpha[1] >=1){
      e1_alpha_tab <- matrix(nrow = length(alpha_01), ncol = 4)
      e1_alpha_tab[,1] <- alpha_01
      e1_alpha_tab[,2] <- alpha_01.se
      e1_alpha_tab[,3] <- e1_alpha_tab[,1]/e1_alpha_tab[,2]
      e1_alpha_tab[,4] <- 1 - pchisq(e1_alpha_tab[,3]**2,1)
      e1_names_tab <- c(e1_names_tab, alpha_01.name)
    }
  }
  else{
    if(x$control$nb.alpha[1] >=2){
      e1_alpha_tab <- matrix(nrow = length(alpha_01)-1, ncol = 4)
      e1_alpha_tab[,1] <- alpha_01[-1]
      e1_alpha_tab[,2] <- alpha_01.se[-1]
      e1_alpha_tab[,3] <- e1_alpha_tab[,1]/e1_alpha_tab[,2]
      e1_alpha_tab[,4] <- 1 - pchisq(e1_alpha_tab[,3]**2,1)
      e1_names_tab <- c(e1_names_tab, alpha_01.name[-1])
    }
  }




  e1_bas_tab <- NULL
  if(x$control$hazard_baseline_01 == "Exponential"){
    e1_bas_tab <- matrix(nrow = 1, ncol = 4)
    e1_bas_tab[1,1] <- alpha_01[1]
    e1_bas_tab[1,2] <- alpha_01.se[1]
    e1_bas_tab[1,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[1,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    #e1_names_tab <- c(e1_names_tab, alpha_01.name[-1])
    rownames(e1_bas_tab) <- c("intercept")
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline_01 == "Weibull"){
    e1_bas_tab <- matrix(nrow = 2, ncol = 4)
    e1_bas_tab[1,1] <- alpha_01[1]
    e1_bas_tab[1,2] <- alpha_01.se[1]
    #e1_bas_tab[1,] <- e1_alpha_tab[1,]
    #e1_alpha_tab <- e1_alpha_tab[-1,]
    #e1_names_tab <- c(e1_names_tab, alpha.name[-1])
    e1_bas_tab[2,1] <- shape_01
    e1_bas_tab[2,2] <- shape_01.se
    e1_bas_tab[,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    rownames(e1_bas_tab) <- c("intercept",shape_01.name)
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline_01 == "Splines"){
    e1_bas_tab <- matrix(nrow = length(gamma_01), ncol = 4)
    e1_bas_tab[,1] <- gamma_01
    e1_bas_tab[,2] <- gamma_01.se
    e1_bas_tab[,3] <- e1_bas_tab[,1]/e1_bas_tab[,2]
    e1_bas_tab[,4] <- 1 - pchisq(e1_bas_tab[,3]**2,1)
    rownames(e1_bas_tab) <- gamma_01.name
    colnames(e1_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  e1_surv_tab <- rbind(e1_var_tab, e1_share_current_tab, e1_share_slope_tab, e1_alpha_tab)
  rownames(e1_surv_tab) <- e1_names_tab
  colnames(e1_surv_tab) <- c("Coeff", "SE", "Wald", "P-value")

  if(nrow(e1_bas_tab)!=0){
    cat("\n")
    cat("       Regression:")
    cat("\n")
    print(e1_surv_tab)
  }
  cat("\n")
  cat(paste("     Baseline: ",x$control$hazard_baseline_01), "\n")
  cat("\n")
  print(e1_bas_tab)


  cat("\n")

  cat("    Transition 0-2:")
  e2_var_tab <- NULL
  e2_share_current_tab <- NULL
  e2_share_slope_tab <- NULL
  e2_alpha_tab <- NULL
  e2_names_tab <- c()


  if(c("current value") %in% x$control$sharedtype_02){
    e2_share_current_tab <- matrix(nrow = 1, ncol = 4)
    e2_share_current_tab[,1] <- alpha.current_02
    e2_share_current_tab[,2] <- alpha.current_02.se
    e2_share_current_tab[,3] <- e2_share_current_tab[,1]/e2_share_current_tab[,2]
    e2_share_current_tab[,4] <- 1 - pchisq(e2_share_current_tab[,3]**2,1)
    e2_names_tab <- c(e2_names_tab, alpha.current_02.name)
  }
  if(c("slope") %in% x$control$sharedtype_02){
    e2_share_slope_tab <- matrix(nrow = 1, ncol = 4)
    e2_share_slope_tab[,1] <- alpha.slope_02
    e2_share_slope_tab[,2] <- alpha.slope_02.se
    e2_share_slope_tab[,3] <- e2_share_slope_tab[,1]/e2_share_slope_tab[,2]
    e2_share_slope_tab[,4] <- 1 - pchisq(e2_share_slope_tab[,3]**2,1)
    e2_names_tab <- c(e2_names_tab, alpha.slope_02.name)
  }
  if(c("inter visit variability") %in% x$control$sharedtype_02){
    e2_share_intervar_tab <- matrix(nrow = 1, ncol = 4)
    e2_share_intervar_tab[,1] <- alpha.intervar_02
    e2_share_intervar_tab[,2] <- alpha.intervar_02.se
    e2_share_intervar_tab[,3] <- e2_share_intervar_tab[,1]/e2_share_intervar_tab[,2]
    e2_share_intervar_tab[,4] <- 1 - pchisq(e2_share_intervar_tab[,3]**2,1)
    e2_names_tab <- c(e2_names_tab, alpha.intervar_02.name)
  }
  if(c("intra visit variability") %in% x$control$sharedtype_02){
    e2_share_intravar_tab <- matrix(nrow = 1, ncol = 4)
    e2_share_intravar_tab[,1] <- alpha.intravar_02
    e2_share_intravar_tab[,2] <- alpha.intravar_02.se
    e2_share_intravar_tab[,3] <- e2_share_intravar_tab[,1]/e2_share_intravar_tab[,2]
    e2_share_intravar_tab[,4] <- 1 - pchisq(e2_share_intravar_tab[,3]**2,1)
    e2_names_tab <- c(e2_names_tab, alpha.intravar_02.name)
  }
  if(x$control$hazard_baseline_02 == "Splines"){
    if(x$control$nb.alpha[2] >=1){
      e2_alpha_tab <- matrix(nrow = length(alpha_02), ncol = 4)
      e2_alpha_tab[,1] <- alpha_02
      e2_alpha_tab[,2] <- alpha_02.se
      e2_alpha_tab[,3] <- e2_alpha_tab[,1]/e2_alpha_tab[,2]
      e2_alpha_tab[,4] <- 1 - pchisq(e2_alpha_tab[,3]**2,1)
      e2_names_tab <- c(e2_names_tab, alpha_02.name)
    }
  }
  else{
    if(x$control$nb.alpha[2] >=2){
      e2_alpha_tab <- matrix(nrow = length(alpha_02)-1, ncol = 4)
      e2_alpha_tab[,1] <- alpha_02[-1]
      e2_alpha_tab[,2] <- alpha_02.se[-1]
      e2_alpha_tab[,3] <- e2_alpha_tab[,1]/e2_alpha_tab[,2]
      e2_alpha_tab[,4] <- 1 - pchisq(e2_alpha_tab[,3]**2,1)
      e2_names_tab <- c(e2_names_tab, alpha_02.name[-1])
    }
  }




  e2_bas_tab <- NULL
  if(x$control$hazard_baseline_02 == "Exponential"){
    e2_bas_tab <- matrix(nrow = 1, ncol = 4)
    e2_bas_tab[,1] <- alpha_02[1]
    e2_bas_tab[,2] <- alpha_02.se[1]
    e2_bas_tab[,3] <- e2_bas_tab[,1]/e2_bas_tab[,2]
    e2_bas_tab[,4] <- 1 - pchisq(e2_bas_tab[,3]**2,1)
    # e2_bas_tab[1,] <- e2_alpha_tab[1,]
    # e2_alpha_tab <- e2_alpha_tab[-1,]
    # e2_names_tab <- c(e2_names_tab, alpha_02.name[-1])
    rownames(e2_bas_tab) <- c("intercept")
    colnames(e2_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline_02 == "Weibull"){
    e2_bas_tab <- matrix(nrow = 2, ncol = 4)
    #e2_bas_tab[1,] <- e2_alpha_tab[1,]
    #e2_alpha_tab <- e2_alpha_tab[-1,]
    #e2_names_tab <- c(e2_names_tab, alpha_02.name[-1])
    e2_bas_tab[,1] <- alpha_02[1]
    e2_bas_tab[,2] <- alpha_02.se[1]
    e2_bas_tab[2,1] <- shape_02
    e2_bas_tab[2,2] <- shape_02.se
    e2_bas_tab[,3] <- e2_bas_tab[,1]/e2_bas_tab[,2]
    e2_bas_tab[,4] <- 1 - pchisq(e2_bas_tab[,3]**2,1)
    rownames(e2_bas_tab) <- c("intercept",shape_02.name)
    colnames(e2_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline_02 == "Splines"){
    e2_bas_tab <- matrix(nrow = length(gamma_02), ncol = 4)
    e2_bas_tab[,1] <- gamma_02
    e2_bas_tab[,2] <- gamma_02.se
    e2_bas_tab[,3] <- e2_bas_tab[,1]/e2_bas_tab[,2]
    e2_bas_tab[,4] <- 1 - pchisq(e2_bas_tab[,3]**2,1)
    rownames(e2_bas_tab) <- gamma_02.name
    colnames(e2_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }

  e2_surv_tab <- rbind(e2_var_tab, e2_share_current_tab, e2_share_slope_tab, e2_alpha_tab)
  rownames(e2_surv_tab) <- e2_names_tab
  colnames(e2_surv_tab) <- c("Coeff", "SE", "Wald", "P-value")

  if(nrow(e2_bas_tab)!=0){
    cat("\n")
    cat("       Regression:")
    cat("\n")
    print(e2_surv_tab)
  }
  cat("\n")
  cat(paste("     Baseline: ",x$control$hazard_baseline_02), "\n")
  cat("\n")
  print(e2_bas_tab)

  cat("\n")

  cat("\n")

  cat("    Transition 1-2:")
  e12_var_tab <- NULL
  e12_share_current_tab <- NULL
  e12_share_slope_tab <- NULL
  e12_alpha_tab <- NULL
  e12_names_tab <- c()


  if(c("current value") %in% x$control$sharedtype_12){
    e12_share_current_tab <- matrix(nrow = 1, ncol = 4)
    e12_share_current_tab[,1] <- alpha.current_12
    e12_share_current_tab[,2] <- alpha.current_12.se
    e12_share_current_tab[,3] <- e12_share_current_tab[,1]/e12_share_current_tab[,2]
    e12_share_current_tab[,4] <- 1 - pchisq(e12_share_current_tab[,3]**2,1)
    e12_names_tab <- c(e12_names_tab, alpha.current_12.name)
  }
  if(c("slope") %in% x$control$sharedtype_12){
    e12_share_slope_tab <- matrix(nrow = 1, ncol = 4)
    e12_share_slope_tab[,1] <- alpha.slope_12
    e12_share_slope_tab[,2] <- alpha.slope_12.se
    e12_share_slope_tab[,3] <- e12_share_slope_tab[,1]/e12_share_slope_tab[,2]
    e12_share_slope_tab[,4] <- 1 - pchisq(e12_share_slope_tab[,3]**2,1)
    e12_names_tab <- c(e12_names_tab, alpha.slope_12.name)
  }
  if(c("inter visit variability") %in% x$control$sharedtype_12){
    e12_share_intervar_tab <- matrix(nrow = 1, ncol = 4)
    e12_share_intervar_tab[,1] <- alpha.intervar_12
    e12_share_intervar_tab[,2] <- alpha.intervar_12.se
    e12_share_intervar_tab[,3] <- e12_share_intervar_tab[,1]/e12_share_intervar_tab[,2]
    e12_share_intervar_tab[,4] <- 1 - pchisq(e12_share_intervar_tab[,3]**2,1)
    e12_names_tab <- c(e12_names_tab, alpha.intervar_12.name)
  }
  if(c("intra visit variability") %in% x$control$sharedtype_12){
    e12_share_intravar_tab <- matrix(nrow = 1, ncol = 4)
    e12_share_intravar_tab[,1] <- alpha.intravar_12
    e12_share_intravar_tab[,2] <- alpha.intravar_12.se
    e12_share_intravar_tab[,3] <- e12_share_intravar_tab[,1]/e12_share_intravar_tab[,2]
    e12_share_intravar_tab[,4] <- 1 - pchisq(e12_share_intravar_tab[,3]**2,1)
    e12_names_tab <- c(e12_names_tab, alpha.intravar_12.name)
  }
  if(x$control$hazard_baseline_12 == "Splines"){
    if(x$control$nb.alpha[3] >=1){
      e12_alpha_tab <- matrix(nrow = length(alpha_12), ncol = 4)
      e12_alpha_tab[,1] <- alpha_12
      e12_alpha_tab[,2] <- alpha_12.se
      e12_alpha_tab[,3] <- e12_alpha_tab[,1]/e12_alpha_tab[,2]
      e12_alpha_tab[,4] <- 1 - pchisq(e12_alpha_tab[,3]**2,1)
      e12_names_tab <- c(e12_names_tab, alpha_12.name)
    }
  }
  else{
    if(x$control$nb.alpha[3] >=2){
      e12_alpha_tab <- matrix(nrow = length(alpha_12)-1, ncol = 4)
      e12_alpha_tab[,1] <- alpha_12[-1]
      e12_alpha_tab[,2] <- alpha_12.se[-1]
      e12_alpha_tab[,3] <- e12_alpha_tab[,1]/e12_alpha_tab[,2]
      e12_alpha_tab[,4] <- 1 - pchisq(e12_alpha_tab[,3]**2,1)
      e12_names_tab <- c(e12_names_tab, alpha_12.name[-1])
    }
  }




  e12_bas_tab <- NULL
  if(x$control$hazard_baseline_12 == "Exponential"){
    e12_bas_tab <- matrix(nrow = 1, ncol = 4)
    e12_bas_tab[,1] <- alpha_12[1]
    e12_bas_tab[,2] <- alpha_12.se[1]
    e12_bas_tab[,3] <- e12_bas_tab[,1]/e12_bas_tab[,2]
    e12_bas_tab[,4] <- 1 - pchisq(e12_bas_tab[,3]**2,1)
    # e2_bas_tab[1,] <- e2_alpha_tab[1,]
    # e2_alpha_tab <- e2_alpha_tab[-1,]
    # e2_names_tab <- c(e2_names_tab, alpha_02.name[-1])
    rownames(e12_bas_tab) <- c("intercept")
    colnames(e12_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline_12 == "Weibull"){
    e12_bas_tab <- matrix(nrow = 2, ncol = 4)
    #e2_bas_tab[1,] <- e2_alpha_tab[1,]
    #e2_alpha_tab <- e2_alpha_tab[-1,]
    #e2_names_tab <- c(e2_names_tab, alpha_02.name[-1])
    e12_bas_tab[,1] <- alpha_12[1]
    e12_bas_tab[,2] <- alpha_12.se[1]
    e12_bas_tab[2,1] <- shape_12
    e12_bas_tab[2,2] <- shape_12.se
    e12_bas_tab[,3] <- e12_bas_tab[,1]/e12_bas_tab[,2]
    e12_bas_tab[,4] <- 1 - pchisq(e12_bas_tab[,3]**2,1)
    rownames(e12_bas_tab) <- c("intercept",shape_12.name)
    colnames(e12_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }
  if(x$control$hazard_baseline_12 == "Splines"){
    e12_bas_tab <- matrix(nrow = length(gamma_12), ncol = 4)
    e12_bas_tab[,1] <- gamma_12
    e12_bas_tab[,2] <- gamma_12.se
    e12_bas_tab[,3] <- e12_bas_tab[,1]/e12_bas_tab[,2]
    e12_bas_tab[,4] <- 1 - pchisq(e12_bas_tab[,3]**2,1)
    rownames(e12_bas_tab) <- gamma_12.name
    colnames(e12_bas_tab) <- c("Coeff", "SE", "Wald", "P-value")
  }

  e12_surv_tab <- rbind(e12_var_tab, e12_share_current_tab, e12_share_slope_tab, e12_alpha_tab)
  rownames(e12_surv_tab) <- e12_names_tab
  colnames(e12_surv_tab) <- c("Coeff", "SE", "Wald", "P-value")

  if(nrow(e12_bas_tab)!=0){
    cat("\n")
    cat("       Regression:")
    cat("\n")
    print(e12_surv_tab)
  }
  cat("\n")
  cat(paste("     Baseline: ",x$control$hazard_baseline_12), "\n")
  cat("\n")
  print(e12_bas_tab)

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
