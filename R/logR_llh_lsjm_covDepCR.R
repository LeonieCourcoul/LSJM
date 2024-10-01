logR_llh_lsjm_covDepCR <- function(param,hazard_baseline_01, sharedtype_01,
                                   hazard_baseline_02, sharedtype_02,
                                   ord.splines, nb.beta, Zq, nb_pointsGK,
                                   nb.e.a, S, wk, rep_wk, sk_GK, nb.alpha,
                                   Matrices,
                                   left_trunc,
                                   knots.hazard_baseline.splines_01,
                                   knots.hazard_baseline.splines_02,
                                   index_beta_slope, index_b_slope, Ind,
                                   nb.omega , nb.e.a.sigma,
                                   correlated_re){


  # Initialiser certains paramètres
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0; Gompertz.1_02 <- 0; Gompertz.2_02 <- 0;
  shape_01 <- 0; shape_02 <- 0;
  alpha.current_01 <- 0; alpha.current_02 <- 0;
  alpha.slope_01 <- 0; alpha.slope_02 <- 0;
  alpha.var_01 <- 0; alpha.var_02 <- 0;
  alpha_01 <- c(0); alpha_02 <- c(0);
  gamma_01 <- c(0,3); gamma_02 <- c(0);
  beta_slope <- c(0); b_y_slope <- as.matrix(1);
  Xslope_T_i <- c(0); Uslope_T_i <- c(0); Xslope_GK_T_i <- as.matrix(1); Uslope_GK_T_i <- as.matrix(1);
  X_GK_T0_i<- as.matrix(1); U_GK_T0_i<- as.matrix(1); Xslope_GK_T0_i<- as.matrix(1); Uslope_GK_T0_i<- as.matrix(1);
  st_T0_i <- c(0); B_T_i_01 <- c(0); B_T_i_02 <- c(0);
  Bs_T_i_01 <- as.matrix(1);Bs_T_i_02 <- as.matrix(1);
  Bs_T0_i_01 <- as.matrix(1); Bs_T0_i_02 <- as.matrix(1); Time_T0_i <- 0;
  st_T_i <- c(0); st_T0_i <- c(0);
  O_T <- as.matrix(0); W_T <- as.matrix(0);
  O_GK_T <- as.matrix(0); W_GK_T <- as.matrix(0); alpha_b_01 <- c(0); alpha_b_02 <- c(0)
  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]**2
    Gompertz.2_01 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+ord.splines[1]+1)]
    curseur <- curseur + ord.splines[1] + 2
  }
  ### Covariables :
  nb.alpha_01 <- nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association

  if("random effects" %in% sharedtype_01){
    alpha_b_01 <- param[curseur:(curseur+nb.e.a-1)]
    curseur <- curseur + nb.e.a
  }
  if("value" %in% sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% sharedtype_01){
    alpha.var_01 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Risque 02
  if(hazard_baseline_02 == "Weibull"){
    shape_02 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline_02 == "Gompertz"){
    Gompertz.1_02 <- param[curseur]**2
    Gompertz.2_02 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_02 == "Splines"){
    gamma_02 <- param[(curseur):(curseur+ord.splines[2]+1)]
    curseur <- curseur + ord.splines[2] + 2
  }
  ### Covariables :
  nb.alpha_02 <- nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("random effects" %in% sharedtype_02){
    alpha_b_02 <- param[curseur:(curseur+nb.e.a-1)]
    curseur <- curseur + nb.e.a
  }
  if("value" %in% sharedtype_02){
    alpha.current_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_02){
    alpha.slope_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability" %in% sharedtype_02){
    alpha.var_02 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  curseur <- curseur+nb.beta
  omega <- param[(curseur):(curseur+nb.omega-1)]
  curseur <- curseur+nb.omega
  if(correlated_re){
    C1 <- matrix(rep(0,(nb.e.a+nb.e.a.sigma)**2),nrow=nb.e.a+nb.e.a.sigma,ncol=nb.e.a+nb.e.a.sigma)
    C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
    MatCov <- C1
    MatCov <- as.matrix(MatCov)
    random.effects <- Zq%*%t(MatCov)
    b_al <- random.effects[,1:nb.e.a]
    b_al <- matrix(b_al, ncol = nb.e.a)
    b_om <- random.effects[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]
    b_om <- matrix(b_om, ncol = nb.e.a.sigma)
  }
  else{
    borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
    C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
    C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
    borne3 <- borne1 + choose(n = nb.e.a.sigma, k = 2) + nb.e.a.sigma
    C3 <- matrix(rep(0,(nb.e.a.sigma)**2),nrow=nb.e.a.sigma,ncol=nb.e.a.sigma)
    C3[lower.tri(C3, diag=T)] <- param[(borne1+1):borne3]
    MatCovb <- as.matrix(C1)
    MatCovSig <- as.matrix(C3)
    b_al <- Zq[,1:nb.e.a]%*%t(MatCovb)
    b_al <- matrix(b_al, ncol = nb.e.a)
    b_om <- Zq[,(nb.e.a+1):(nb.e.a+nb.e.a.sigma)]%*%t(MatCovSig)
    b_om <- matrix(b_om, ncol = nb.e.a.sigma)
  }
  if("slope" %in% sharedtype_01 || "slope" %in% sharedtype_02 ){
    b_y_slope <- as.matrix(b_al[,index_b_slope], ncol = length(index_b_slope))
    beta_slope <- beta[index_beta_slope]
  }
  ll_glob <- rep(NA, Ind)

  # Creations entrees rcpp
  sharedtype <- c("value" %in% sharedtype_01, "slope" %in% sharedtype_01, "variability" %in% sharedtype_01,
                  "value" %in% sharedtype_02, "slope" %in% sharedtype_02, "variability" %in% sharedtype_02,
                  "random effects" %in% sharedtype_01, "random effects" %in% sharedtype_02
  )
  HB <- list(hazard_baseline_01, hazard_baseline_02)
  Weibull <- c(shape_01, shape_02)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01, Gompertz.1_02, Gompertz.2_02)
  alpha_y_slope <- c(alpha.current_01,alpha.current_02, alpha.slope_01,alpha.slope_02)
  alpha_var <- c(alpha.var_01, alpha.var_02)
  alpha_z <- list(alpha_01, alpha_02)
  nb_points_integral <- c(S, nb_pointsGK)
  gamma_z0 <- list(gamma_01, gamma_02)

  delta1 <- Matrices[["delta1"]];delta2 <- Matrices[["delta2"]];  Z_01 <- Matrices[["Z_01"]]; Z_02 <- Matrices[["Z_02"]]
  Time_T <- Matrices[["Time_T"]];
  st_T <- Matrices[["st_T"]];
  X_GK_T <- Matrices[["X_GK_T"]];U_GK_T <- Matrices[["U_GK_T"]]
  Xslope_GK_T <- Matrices[["Xslope_GK_T"]];Uslope_GK_T <- Matrices[["Uslope_GK_T"]]
  X_T <- Matrices[["X_T"]]; U_T <- Matrices[["U_T"]]; Xslope_T <- Matrices[["Xslope_T"]]; Uslope_T <- Matrices[["Uslope_T"]]
  X_base <- Matrices[["X_base"]]; U_base <- Matrices[["U_base"]]; y.new <- Matrices[["y.new"]]
  ID.visit <- Matrices[["ID.visit"]]; offset <- Matrices[["offset"]];
  B_T_01 <- Matrices[["B_T_01"]]; Bs_T_01 <- Matrices[["Bs_T_01"]];
  B_T_02 <- Matrices[["B_T_02"]]; Bs_T_02 <- Matrices[["Bs_T_02"]];

  O_base <- Matrices[["O_base"]];  W_base <- Matrices[["W_base"]];
  O_T <- Matrices[["O_T"]];  W_T <- Matrices[["W_T"]]; O_GK_T <- Matrices[["O_GK_T"]];  W_GK_T <- Matrices[["W_GK_T"]];


    Time_T0 <- Matrices[["Time_T0"]]; st_T0 <- Matrices[["st_T0"]]; Bs_T0_01 <- Matrices[["Bs_T0_01"]]; Bs_T0_02 <- Matrices[["Bs_T0_02"]]
    X_GK_T0 <- Matrices[["X_GK_T0"]];U_GK_T0 <- Matrices[["U_GK_T0"]]
    Xslope_GK_T0 <- Matrices[["Xslope_GK_T0"]];Uslope_GK_T0 <- Matrices[["Uslope_GK_T0"]]
    O_GK_T0 <- Matrices[["O_GK_T0"]];  W_GK_T0 <- Matrices[["W_GK_T0"]];

  ll_glob <- log_llh_lsjm_covDepCR(sharedtype,  HB,  Gompertz,  Weibull,
                                    nb_points_integral,
                                    alpha_y_slope, alpha_var, t(alpha_b_01), t(alpha_b_02), alpha_z,  gamma_z0,  beta,  beta_slope, omega,
                                    b_al,  b_y_slope, b_om, wk,
                                    delta1, delta2,  Z_01,  Z_02,  X_T,  U_T,
                                    Xslope_T,  Uslope_T,  X_GK_T,  U_GK_T,  Xslope_GK_T,
                                    Uslope_GK_T,
                                    X_GK_T0,  U_GK_T0,  Xslope_GK_T0,  Uslope_GK_T0,
                                    Time_T,  Time_T0, st_T,  st_T0,
                                    Bs_T0_01,  Bs_T0_02, left_trunc,
                                    X_base,  U_base,   y.new,  Ind,
                                    offset, B_T_01 , Bs_T_01 ,
                                    B_T_02 , Bs_T_02,

                                    O_base, W_base, O_T, W_T, O_GK_T, W_GK_T,
                                    O_GK_T0, W_GK_T0


  )


  ll_glob2 <- sum(ll_glob)
  if(is.na(ll_glob2) || ll_glob2>0){
    #print(param)
    #print(ll_glob2)
    ll_glob2 <- -1E09
  }
  ll_glob2

}
