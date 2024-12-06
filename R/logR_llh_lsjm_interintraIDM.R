logR_llh_lsjm_interintraIDM <- function(param,hazard_baseline_01, sharedtype_01,
                                       hazard_baseline_02, sharedtype_02,
                                       hazard_baseline_12, sharedtype_12,
                                       ord.splines, nb.beta, Zq, nb_pointsGK,
                                       nb.e.a, S, wk, rep_wk, sk_GK, nb.alpha,
                                       variability_inter_visit, variability_intra_visit,
                                       correlated_re, Case1, Case1bis, Case2, Case3,
                                       nbCase1, nbCase1bis, nbCase2, nbCase3, left_trunc,
                                       knots.hazard_baseline.splines_01,
                                       knots.hazard_baseline.splines_02,
                                       knots.hazard_baseline.splines_12,
                                       index_beta_slope, index_b_slope,
                                       control){


  # Initialiser certains paramètres
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0; Gompertz.1_02 <- 0; Gompertz.2_02 <- 0; Gompertz.1_12 <- 0; Gompertz.2_12 <- 0
  shape_01 <- 0; shape_02 <- 0; shape_12 <- 0
  alpha.inter_01 <- 0; alpha.inter_02 <- 0; alpha.inter_12 <- 0
  alpha.intra_01 <- 0; alpha.intra_02 <- 0; alpha.intra_12 <- 0
  alpha.current_01 <- 0; alpha.current_02 <- 0; alpha.current_12 <- 0
  alpha.slope_01 <- 0; alpha.slope_02 <- 0; alpha.slope_12 <- 0
  alpha_01 <- c(0); alpha_02 <- c(0); alpha_12 <- c(0);
  gamma_01 <- c(0,3); gamma_02 <- c(0); gamma_12 <- c(0);
  beta_slope <- c(0); b_y_slope <- as.matrix(1); sigma_inter <- c(0); sigma_intra <- c(0);
  Xslope_T_i <- c(0); Uslope_T_i <- c(0); Xslope_GK_T_i <- as.matrix(1); Uslope_GK_T_i <- as.matrix(1);
  Xslope_L_i <- c(0); Uslope_L_i <- c(0); Xslope_GK_L_i <- as.matrix(1); Uslope_GK_L_i <- as.matrix(1);
  Xslope_GK_L_R_i <- as.matrix(1); Uslope_GK_L_R_i <- as.matrix(1);
  Xslope_GK_L_T_i <- as.matrix(1); Uslope_GK_L_T_i <- as.matrix(1);
  Xslope_GK_0_LR_i <- as.matrix(1); Uslope_GK_0_LR_i <- as.matrix(1); X_GK_T0_i<- as.matrix(1); U_GK_T0_i<- as.matrix(1); Xslope_GK_T0_i<- as.matrix(1); Uslope_GK_T0_i<- as.matrix(1);
  Xslope_GK_0_LT_i <- as.matrix(1); Uslope_GK_0_LT_i <- as.matrix(1);
  st_T0_i <- c(0); B_T_i_01 <- c(0); B_T_i_02 <- c(0); B_T_i_12 <- c(0); B_L_i_02 <- c(0); B_L_i_12 <- c(0);
  Bs_T_i_01 <- as.matrix(1);Bs_T_i_02 <- as.matrix(1); Bs_T_i_12 <- as.matrix(1); Bs_L_i_01 <- as.matrix(1);Bs_L_i_02 <- as.matrix(1); Bs_L_i_12 <- as.matrix(1);
  Bs_0_LR_i_01 <- as.matrix(1); Bs_0_LR_i_02 <- as.matrix(1); Bs_0_LR_i_12 <- as.matrix(1);
  Bs_0_LT_i_01 <- as.matrix(1); Bs_0_LT_i_02 <- as.matrix(1); Bs_0_LT_i_12 <- as.matrix(1);
  Bs_L_R_i_01 <- as.matrix(1); Bs_L_R_i_02 <- as.matrix(1); Bs_L_R_i_12<- as.matrix(1);
  Bs_L_T_i_01 <- as.matrix(1); Bs_L_T_i_02 <- as.matrix(1); Bs_L_T_i_12<- as.matrix(1);
  Bs_T0_i_01 <- as.matrix(1); Bs_T0_i_02 <- as.matrix(1); Bs_T0_i_12 <- as.matrix(1); Time_T0_i <- 0;
  st_T_i <- c(0); st_L_i <- c(0); st_0_LR_i <- as.matrix(1); st_L_R_i <- c(0); st_T0_i <- c(0); st_0_LT_i <- as.matrix(1); st_L_T_i <- c(0);
  B_L_i_01 <- c(0);  alpha_b_01 <- c(0); alpha_b_02 <- c(0); alpha_b_12 <- c(0)
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
  if("variability inter" %in% sharedtype_01){
    alpha.inter_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability intra" %in% sharedtype_01){
    alpha.intra_01 <- param[curseur]
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
  if("variability inter" %in% sharedtype_02){
    alpha.inter_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability intra" %in% sharedtype_02){
    alpha.intra_02 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Risque 12
  if(hazard_baseline_12 == "Weibull"){
    shape_12 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline_12 == "Gompertz"){
    Gompertz.1_12 <- param[curseur]**2
    Gompertz.2_12 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_12 == "Splines"){
    gamma_12 <- param[(curseur):(curseur+ord.splines[3]+1)]
    curseur <- curseur + ord.splines[3] + 2
  }
  ### Covariables :
  nb.alpha_12 <- nb.alpha[3]
  if(nb.alpha_12 >=1){
    alpha_12 <- param[(curseur):(curseur+nb.alpha_12-1)]
    curseur <- curseur+nb.alpha_12
  }
  ### Association
  if("random effects" %in% sharedtype_12){
    alpha_b_12 <- param[curseur:(curseur+nb.e.a-1)]
    curseur <- curseur + nb.e.a
  }
  if("value" %in% sharedtype_12){
    alpha.current_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_12){
    alpha.slope_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability inter" %in% sharedtype_12){
    alpha.inter_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("variability intra" %in% sharedtype_12){
    alpha.intra_12 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Marker
  ### Fixed effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  if( "slope" %in% sharedtype_01 || "slope" %in% sharedtype_02 || "slope" %in% sharedtype_12){
    beta_slope <- beta[index_beta_slope]
  }
  curseur <- curseur+nb.beta
  if(variability_inter_visit){
    mu.inter <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.inter <- param[curseur]
    curseur <- curseur +1
  }
  if(variability_intra_visit){
    mu.intra <- param[curseur]
    curseur <- curseur + 1
  }
  else{
    sigma.epsilon.intra <- param[curseur]
    curseur <- curseur +1
  }
  ### Cholesky matrix for random effects
  if(variability_inter_visit && variability_intra_visit){
    if(correlated_re){

      C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }
    else{
      borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      C3 <- matrix(rep(0,2*nb.e.a), ncol = nb.e.a)
      C4 <- matrix(rep(0,2*nb.e.a), nrow = nb.e.a)
      Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
      Cholesky <- as.matrix(Cholesky)
    }
  }
  else{
    if(variability_inter_visit){
      if(correlated_re){
        C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
        Cholesky <- C1
      }
      else{
        borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
        C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        C3 <- matrix(rep(0,nb.e.a), ncol = 1)
        C4 <- matrix(rep(0,nb.e.a), nrow = 1)
        Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
        Cholesky <- as.matrix(Cholesky)
      }
    }
    else{
      if(variability_intra_visit){
        if(correlated_re){
          C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
          C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
          Cholesky <- C1
        }
        else{
          borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
          C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
          C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
          C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
          C3 <- matrix(rep(0,nb.e.a), ncol = 1)
          C4 <- matrix(rep(0,nb.e.a), nrow = 1)
          Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
          Cholesky <- as.matrix(Cholesky)
        }
      }
    }
    if(!variability_inter_visit && !variability_intra_visit){
      C1 <- matrix(rep(0,(length(param)-curseur)**2),nrow=length(param)-curseur,ncol=length(param)-curseur)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }

  }
  # Manage random effects
  random.effects <- Zq%*%t(Cholesky)
  b_y <- random.effects[,1:nb.e.a]
  b_y <- matrix(b_y, ncol = nb.e.a)
  if("slope" %in% sharedtype_01 || "slope" %in% sharedtype_02 || "slope" %in% sharedtype_12 ){
    b_y_slope <- as.matrix(b_y[,index_b_slope], ncol = length(index_b_slope))
    beta_slope <- beta[index_beta_slope]
  }
  if(variability_inter_visit){
    b_inter <- random.effects[,nb.e.a+1]
    sigma_inter <- exp(mu.inter + b_inter)
    var.inter <- sigma_inter**2
  }
  else{
    sigma_inter <- rep(sigma.epsilon.inter,S)
    var.inter <- rep(sigma.epsilon.inter**2,S)
  }
  if(variability_intra_visit){
    b_intra <- random.effects[,ncol(random.effects)]
    sigma_intra <- exp(mu.intra + b_intra)
    var.intra <- sigma_intra**2
  }
  else{
    sigma_intra <- rep(sigma.epsilon.intra,S)
    var.intra <- rep(sigma.epsilon.intra**2,S)
  }
  sigma_inter_intra <- list(sigma_inter, sigma_intra, var.inter+var.intra, var.inter, var.intra, var.intra*(2*var.inter+var.intra))
  ll_glob <- rep(NA, nbCase1 + nbCase1bis + nbCase2 + nbCase3)

  # Creations entrees rcpp
  sharedtype <- c("value" %in% sharedtype_01, "slope" %in% sharedtype_01, "variability inter" %in% sharedtype_01, "variability intra" %in% sharedtype_01,
                  "value" %in% sharedtype_02, "slope" %in% sharedtype_02, "variability inter" %in% sharedtype_02, "variability intra" %in% sharedtype_02,
                  "value" %in% sharedtype_12, "slope" %in% sharedtype_12, "variability inter" %in% sharedtype_12, "variability intra" %in% sharedtype_12,
                  "random effects" %in% sharedtype_01, "random effects" %in% sharedtype_02, "random effects" %in% sharedtype_12)
  HB <- list(hazard_baseline_01, hazard_baseline_02, hazard_baseline_12)
  Weibull <- c(shape_01, shape_02, shape_12)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01, Gompertz.1_02, Gompertz.2_02, Gompertz.1_12, Gompertz.2_12)
  alpha_inter_intra <- c(alpha.inter_01,alpha.inter_02,alpha.inter_12, alpha.intra_01,alpha.intra_02,alpha.intra_12)
  alpha_y_slope <- c(alpha.current_01,alpha.current_02,alpha.current_12, alpha.slope_01,alpha.slope_02,alpha.slope_12)
  alpha_z <- list(alpha_01, alpha_02, alpha_12)
  nb_points_integral <- c(S, nb_pointsGK)
  gamma_z0 <- list(gamma_01, gamma_02, gamma_12)


  ll_glob <- rep(NA, nbCase1 + nbCase1bis + nbCase2 + nbCase3)
  alpha_b <- list(t(alpha_b_01), t(alpha_b_02), t(alpha_b_12))


  if(nbCase1 != 0){
    delta2 <- Case1[["delta2"]]; Z_12 <- Case1[["Z_12"]]; Time_T <- Case1[["Time_T"]]; st_T <- Case1[["st_T"]]
    X_GK_T <- Case1[["X_GK_T"]]; U_GK_T <- Case1[["U_GK_T"]]
    Xslope_GK_T <- Case1[["Xslope_GK_T"]];Uslope_GK_T <- Case1[["Uslope_GK_T"]]
    X_T <- Case1[["X_T"]]; U_T <- Case1[["U_T"]]; Xslope_T <- Case1[["Xslope_T"]]; Uslope_T <- Case1[["Uslope_T"]]
    X_base <- Case1[["X_base"]]; U_base <- Case1[["U_base"]]; y.new <- Case1[["y.new"]];
    ID.visit <- Case1[["ID.visit"]]; offset <- Case1[["offset"]];
    B_T_01 <- Case1[["B_T_01"]]; Bs_T_01 <- Case1[["Bs_T_01"]]
    B_T_02 <- Case1[["B_T_02"]]; Bs_T_02 <- Case1[["Bs_T_02"]]
    B_T_12 <- Case1[["B_T_12"]]; Bs_T_12 <- Case1[["Bs_T_12"]]
    if(left_trunc){
      Time_T0 <- Case1[["Time_T0"]]; st_T0 <- Case1[["st_T0"]]; Bs_T0_01 <- Case1[["Bs_T0_01"]]
      Bs_T0_02 <- Case1[["Bs_T0_02"]];
      X_GK_T0 <- Case1[["X_GK_T0"]];U_GK_T0 <- Case1[["U_GK_T0"]]
      Xslope_GK_T0 <- Case1[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case1[["Uslope_GK_T0"]]
    }
    ### Intégrale L_i to R_i :
    Z_01 <- Case1[["Z_01"]]; Z_02 <- Case1[["Z_02"]];
    X_GK_L_R <- Case1[["X_GK_L_R"]]; U_GK_L_R <- Case1[["U_GK_L_R"]]
    Xslope_GK_L_R <- Case1[["Xslope_GK_L_R"]];Uslope_GK_L_R <- Case1[["Uslope_GK_L_R"]]
    Bs_L_R_01 <- Case1[["Bs_L_R_01"]]; Bs_L_R_02 <- Case1[["Bs_L_R_02"]]; Bs_L_R_12 <- Case1[["Bs_L_R_12"]];
    Time_L_R <- Case1[["Time_L_R"]]; st_L_R <- Case1[["st_L_R"]] #Time_L_R = R-L # st_L_R = (xk+1)/2*(R-L)+L
    ### Integrale de 0 à u (u dans L_i-R_i)
    st_0_LR <- Case1[["st_0_LR"]]; X_0_LR <- Case1[["X_0_LR"]]; U_0_LR <- Case1[["U_0_LR"]];
    Xslope_0_LR <- Case1[["Xslope_0_LR"]]; Uslope_0_LR <- Case1[["Uslope_0_LR"]]; Time_L <- Case1[["Time_L"]];
    Bs_0_LR_01 <- Case1[["Bs_0_LR_01"]]; Bs_0_LR_02 <- Case1[["Bs_0_LR_02"]]; Bs_0_LR_12 <- Case1[["Bs_0_LR_12"]];
    offset_ID <- Case1[["offset_ID"]]; len_visit <-  Case1[["len_visit"]]; offset_position <- Case1[["offset_position"]]
    ck = list(sk_GK = sk_GK, Time_L = Time_L, wk = wk, rep_wk = rep_wk)

    ll_glob[1:nbCase1] <- log_llh_lsjm_interintraIDM_C1(sharedtype, HB,  Gompertz,  Weibull,
                                                nb_points_integral,  alpha_inter_intra,
                                                alpha_y_slope, alpha_b, alpha_z,  gamma_z0,  beta, beta_slope,
                                                b_y,  b_y_slope,  sigma_inter_intra,
                                                delta2,  Z_01,  Z_02,  Z_12, X_T,  U_T,
                                                Xslope_T,  Uslope_T,  X_GK_T,  U_GK_T,  Xslope_GK_T,
                                                Uslope_GK_T,  X_GK_L_R,  U_GK_L_R,  Xslope_GK_L_R,  Uslope_GK_L_R,
                                                X_0_LR,  U_0_LR,  Xslope_0_LR,  Uslope_0_LR,
                                                X_GK_T0,  U_GK_T0,  Xslope_GK_T0,  Uslope_GK_T0,
                                                Time_T,  Time_L_R,  Time_T0, st_T,  st_0_LR,  st_L_R,  st_T0,
                                                B_T_12,
                                                Bs_T_12,
                                                Bs_0_LR_01,  Bs_0_LR_02,  Bs_0_LR_12,
                                                Bs_L_R_01,
                                                Bs_T0_01,  Bs_T0_02,   left_trunc,
                                                len_visit,  X_base,  U_base,   y.new,  offset_ID,  nbCase1,
                                                offset,  offset_position,  ck
    )


  }
  if(nbCase1bis != 0){
    delta2 <- Case1bis[["delta2"]]; Z_12 <- Case1bis[["Z_12"]]; Z_01 <- Case1bis[["Z_01"]]; Z_02 <- Case1bis[["Z_02"]]
    Time_T <- Case1bis[["Time_T"]];Time_L <- Case1bis[["Time_L"]]
    st_T <- Case1bis[["st_T"]];st_L <- Case1bis[["st_L"]]
    X_GK_L <- Case1bis[["X_GK_L"]];U_GK_L <- Case1bis[["U_GK_L"]]
    Xslope_GK_L <- Case1bis[["Xslope_GK_L"]];Uslope_GK_L <- Case1bis[["Uslope_GK_L"]]
    X_GK_T <- Case1bis[["X_GK_T"]];U_GK_T <- Case1bis[["U_GK_T"]]
    Xslope_GK_T <- Case1bis[["Xslope_GK_T"]];Uslope_GK_T <- Case1bis[["Uslope_GK_T"]]
    X_T <- Case1bis[["X_T"]]; U_T <- Case1bis[["U_T"]]; Xslope_T <- Case1bis[["Xslope_T"]]; Uslope_T <- Case1bis[["Uslope_T"]]
    X_L <- Case1bis[["X_L"]]; U_L <- Case1bis[["U_L"]]; Xslope_L <- Case1bis[["Xslope_L"]]; Uslope_L <- Case1bis[["Uslope_L"]]
    X_base <- Case1bis[["X_base"]]; U_base <- Case1bis[["U_base"]]; y.new <- Case1bis[["y.new"]]
    ID.visit <- Case1bis[["ID.visit"]]; offset <- Case1bis[["offset"]];
    B_T_01 <- Case1bis[["B_T_01"]]; Bs_T_01 <- Case1bis[["Bs_T_01"]]; B_L_01 <- Case1bis[["B_L_01"]]; Bs_L_01 <- Case1bis[["Bs_L_01"]]
    B_T_02 <- Case1bis[["B_T_02"]]; Bs_T_02 <- Case1bis[["Bs_T_02"]]; B_L_02 <- Case1bis[["B_L_02"]]; Bs_L_02 <- Case1bis[["Bs_L_02"]]
    B_T_12 <- Case1bis[["B_T_12"]]; Bs_T_12 <- Case1bis[["Bs_T_12"]]; B_L_12 <- Case1bis[["B_L_12"]]; Bs_L_12 <- Case1bis[["Bs_L_12"]]
    if(left_trunc){
      Time_T0 <- Case1bis[["Time_T0"]]; st_T0 <- Case1bis[["st_T0"]]; Bs_T0_01 <- Case1bis[["Bs_T0_01"]]; Bs_T0_02 <- Case1bis[["Bs_T0_02"]]
      X_GK_T0 <- Case1bis[["X_GK_T0"]];U_GK_T0 <- Case1bis[["U_GK_T0"]]
      Xslope_GK_T0 <- Case1bis[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case1bis[["Uslope_GK_T0"]]
    }

    offset_ID <- Case1bis[["offset_ID"]]; len_visit <-  Case1bis[["len_visit"]]; offset_position <- Case1bis[["offset_position"]]

    ll_glob[(nbCase1+1):(nbCase1 + nbCase1bis)] <- log_llh_lsjm_interintraIDM_C1bis( sharedtype,  HB,  Gompertz,  Weibull,
                                                                             nb_points_integral,  alpha_inter_intra,
                                                                             alpha_y_slope, t(alpha_b_01), t(alpha_b_02), t(alpha_b_12),  alpha_z,  gamma_z0,  beta,  beta_slope,
                                                                             b_y,  b_y_slope,  wk,  sigma_inter_intra,
                                                                             delta2,  Z_01,  Z_02,  Z_12,  X_T,  U_T,
                                                                             Xslope_T,  Uslope_T,  X_GK_T,  U_GK_T,  Xslope_GK_T,
                                                                             Uslope_GK_T,  X_L,  U_L,
                                                                             Xslope_L,  Uslope_L,  X_GK_L,  U_GK_L,  Xslope_GK_L,
                                                                             Uslope_GK_L,
                                                                             X_GK_T0,  U_GK_T0,  Xslope_GK_T0,  Uslope_GK_T0,
                                                                             Time_T,  Time_L,  Time_T0, st_T,  st_L,  st_T0,
                                                                             B_T_12,
                                                                             Bs_T_12,
                                                                             B_L_01,
                                                                             Bs_L_01,  Bs_L_02,  Bs_L_12,
                                                                             Bs_T0_01,  Bs_T0_02, left_trunc,
                                                                             len_visit,  X_base,  U_base,   y.new,  offset_ID,  nbCase1bis,
                                                                             offset,  offset_position
    )

  }

  if(nbCase2 != 0){
    delta2 <- Case2[["delta2"]]; Z_01 <- Case2[["Z_01"]]; Z_02 <- Case2[["Z_02"]]
    Time_T <- Case2[["Time_T"]]; st_T <- Case2[["st_T"]];
    X_GK_T <- Case2[["X_GK_T"]];U_GK_T <- Case2[["U_GK_T"]]
    Xslope_GK_T <- Case2[["Xslope_GK_T"]];Uslope_GK_T <- Case2[["Uslope_GK_T"]]
    X_T <- Case2[["X_T"]]; U_T <- Case2[["U_T"]]; Xslope_T <- Case2[["Xslope_T"]]; Uslope_T <- Case2[["Uslope_T"]]
    X_base <- Case2[["X_base"]]; U_base <- Case2[["U_base"]]; y.new <- Case2[["y.new"]]
    ID.visit <- Case2[["ID.visit"]]; offset <- Case2[["offset"]];
    B_T_01 <- Case2[["B_T_01"]]; Bs_T_01 <- Case2[["Bs_T_01"]]
    B_T_02 <- Case2[["B_T_02"]]; Bs_T_02 <- Case2[["Bs_T_02"]]
    if(left_trunc){
      Time_T0 <- Case2[["Time_T0"]]; st_T0 <- Case2[["st_T0"]]; Bs_T0_01 <- Case2[["Bs_T0_01"]]; Bs_T0_02 <- Case2[["Bs_T0_02"]]
      X_GK_T0 <- Case2[["X_GK_T0"]];U_GK_T0 <- Case2[["U_GK_T0"]]
      Xslope_GK_T0 <- Case2[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case2[["Uslope_GK_T0"]]
    }
    offset_ID <- Case2[["offset_ID"]]; len_visit <-  Case2[["len_visit"]]; offset_position <- Case2[["offset_position"]]

    ll_glob[(nbCase1 + nbCase1bis + 1):(nbCase1 + nbCase1bis + nbCase2)] <- log_llh_lsjm_interintraIDM_C2( sharedtype,  HB,  Gompertz,  Weibull,
                                                                                                   nb_points_integral,  alpha_inter_intra,
                                                                                                   alpha_y_slope, t(alpha_b_01), t(alpha_b_02),  alpha_z,  gamma_z0,  beta,  beta_slope,
                                                                                                   b_y,  b_y_slope,  wk,  sigma_inter_intra,
                                                                                                   delta2,  Z_01,  Z_02,  X_T,  U_T,
                                                                                                   Xslope_T,  Uslope_T,  X_GK_T,  U_GK_T,  Xslope_GK_T,
                                                                                                   Uslope_GK_T,
                                                                                                   X_GK_T0,  U_GK_T0,  Xslope_GK_T0,  Uslope_GK_T0,
                                                                                                   Time_T,   Time_T0, st_T,   st_T0,
                                                                                                   B_T_02,
                                                                                                   Bs_T_01,  Bs_T_02,
                                                                                                   Bs_T0_01,  Bs_T0_02,  left_trunc,
                                                                                                   len_visit,  X_base,  U_base,   y.new,  offset_ID,  nbCase2,
                                                                                                   offset,  offset_position
    )
  }


  if(nbCase3 != 0){
    delta2 <- Case3[["delta2"]]; Z_01 <- Case3[["Z_01"]]; Z_02 <- Case3[["Z_02"]]; Z_12 <- Case3[["Z_12"]]
    Time_T <- Case3[["Time_T"]]; st_T <- Case3[["st_T"]]; Time_L <- Case3[["Time_L"]]
    X_GK_T <- Case3[["X_GK_T"]];U_GK_T <- Case3[["U_GK_T"]]
    Xslope_GK_T <- Case3[["Xslope_GK_T"]];Uslope_GK_T <- Case3[["Uslope_GK_T"]]
    X_T <- Case3[["X_T"]]; U_T <- Case3[["U_T"]]; Xslope_T <- Case3[["Xslope_T"]]; Uslope_T <- Case3[["Uslope_T"]]
    X_base <- Case3[["X_base"]]; U_base <- Case3[["U_base"]]; y.new <- Case3[["y.new"]]
    ID.visit <- Case3[["ID.visit"]]; offset <- Case3[["offset"]];
    B_T_01 <- Case3[["B_T_01"]]; Bs_T_01 <- Case3[["Bs_T_01"]]
    B_T_02 <- Case3[["B_T_02"]]; Bs_T_02 <- Case3[["Bs_T_02"]]
    B_T_12 <- Case3[["B_T_12"]]; Bs_T_12 <- Case3[["Bs_T_12"]]
    if(left_trunc){
      Time_T0 <- Case3[["Time_T0"]]; st_T0 <- Case3[["st_T0"]]; Bs_T0_01 <- Case3[["Bs_T0_01"]]; Bs_T0_02 <- Case3[["Bs_T0_02"]]
      X_GK_T0 <- Case3[["X_GK_T0"]];U_GK_T0 <- Case3[["U_GK_T0"]]
      Xslope_GK_T0 <- Case3[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case3[["Uslope_GK_T0"]]
    }
    ### Partie dépendante de l'intégrale :
    X_GK_L_T <- Case3[["X_GK_L_T"]]; U_GK_L_T <- Case3[["U_GK_L_T"]]
    Xslope_GK_L_T <- Case3[["Xslope_GK_L_T"]];Uslope_GK_L_T<- Case3[["Uslope_GK_L_T"]]
    Bs_L_T_01 <- Case3[["Bs_L_T_01"]]; Time_L_T <- Case3[["Time_L_T"]]; st_L_T <- Case3[["st_L_T"]] #Time_L_R = R-L # st_L_R = (xk+1)/2*(R-L)+L
    ### Integrale de 0 à u (u dans L_i-T_i)
    st_0_LT <- Case3[["st_0_LT"]]; X_0_LT <- Case3[["X_0_LT"]]; U_0_LT <- Case3[["U_0_LT"]];
    Xslope_0_LT <- Case3[["Xslope_0_LT"]]; Uslope_0_LT <- Case3[["Uslope_0_LT"]];
    Bs_0_LT_01 <- Case3[["Bs_0_LT_01"]]; Bs_0_LT_02 <- Case3[["Bs_0_LT_02"]]; Bs_0_LT_12 <- Case3[["Bs_0_LT_12"]];

    offset_ID <- Case3[["offset_ID"]]; len_visit <-  Case3[["len_visit"]]; offset_position <- Case3[["offset_position"]]
    ck = list(sk_GK = sk_GK, Time_L = Time_L, nbCase3 = nbCase3, wk = wk, rep_wk = rep_wk)
    list_Times = list(Time_T, Time_L_T, Time_T0)
    ll_glob[(nbCase1 + nbCase1bis + nbCase2 + 1):(nbCase1 + nbCase1bis + nbCase2 + nbCase3)] <- log_llh_lsjm_interintraIDM_C3( sharedtype,  HB,  Gompertz,  Weibull,
                                                                                                                       nb_points_integral,  alpha_inter_intra,
                                                                                                                       alpha_y_slope,  alpha_b, alpha_z,  gamma_z0,  beta,  beta_slope,
                                                                                                                       b_y,  b_y_slope,  sigma_inter_intra,
                                                                                                                       delta2,  Z_01,  Z_02,  Z_12,  X_T,  U_T,
                                                                                                                       Xslope_T,  Uslope_T,  X_GK_T,  U_GK_T,  Xslope_GK_T,
                                                                                                                       Uslope_GK_T,  X_GK_L_T,  U_GK_L_T,  Xslope_GK_L_T,  Uslope_GK_L_T,
                                                                                                                       X_0_LT,  U_0_LT,  Xslope_0_LT,  Uslope_0_LT,
                                                                                                                       X_GK_T0,  U_GK_T0,  Xslope_GK_T0,  Uslope_GK_T0,
                                                                                                                       list_Times,  st_T,  st_0_LT,  st_L_T,  st_T0,
                                                                                                                       B_T_02,  B_T_12,
                                                                                                                       Bs_T_01,  Bs_T_02,  Bs_T_12,
                                                                                                                       Bs_0_LT_01,  Bs_0_LT_02,  Bs_0_LT_12,
                                                                                                                       Bs_L_T_01,
                                                                                                                       Bs_T0_01,  Bs_T0_02,  left_trunc,
                                                                                                                       len_visit,  X_base,  U_base,   y.new,  offset_ID,
                                                                                                                       offset,  offset_position,  ck
    )

  }

  ll_glob2 <- sum(ll_glob)
  if(is.na(ll_glob2) || ll_glob2>0){
    #print(param)
    #print(ll_glob2)
    ll_glob2 <- -1E09
  }
  ll_glob2

}
