#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <stdlib.h> /* srand, rand */
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]

arma::vec log_llh_lsjm_covDepIDM_C3(arma::vec sharedtype, List HB, arma::vec W_G,
                                     arma::vec alpha_y_slope_var, List alpha_z, List gamma_B, List fixed_par,
                                     arma::mat b_y, arma::mat b_y_slope, arma::mat b_om,
                                     arma::mat Z_01, arma::mat Z_02, arma::mat Z_12, arma::mat X_T, arma::mat U_T,
                                     arma::mat Xslope_T, arma::mat Uslope_T, arma::mat X_GK_T, arma::mat U_GK_T, arma::mat Xslope_GK_T,
                                     arma::mat Uslope_GK_T, arma::mat X_GK_L_T, arma::mat U_GK_L_T, arma::mat Xslope_GK_L_T, arma::mat Uslope_GK_L_T,
                                     arma::mat X_GK_0_LT, arma::mat U_GK_0_LT, arma::mat Xslope_GK_0_LT, arma::mat Uslope_GK_0_LT,
                                     arma::mat X_GK_T0, arma::mat U_GK_T0, arma::mat Xslope_GK_T0, arma::mat Uslope_GK_T0,
                                     List list_Times, arma::mat st_T, arma::mat st_0_LT, arma::mat st_L_T, arma::mat st_T0,
                                     arma::mat B_T_02, arma::mat B_T_12,
                                     arma::mat Bs_T_01, arma::mat Bs_T_02, arma::mat Bs_T_12,
                                     arma::mat Bs_0_LT_01, arma::mat Bs_0_LT_02, arma::mat Bs_0_LT_12,
                                     arma::mat Bs_L_T_01,
                                     arma::mat Bs_T0_01, arma::mat Bs_T0_02,
                                     arma::mat X_base, arma::mat U_base,  List longi, List ck,
                                     arma::mat O_T,  arma::mat W_T, arma::mat O_GK_T, arma::mat W_GK_T,
                                     arma::mat O_GK_L_T, arma::mat W_GK_L_T,
                                     arma::mat O_GK_0_LT, arma::mat W_GK_0_LT, arma::mat O_GK_T0, arma::mat W_GK_T0,
                                     arma::mat O_base, arma::mat W_base
){

  //Rcout << "The value of v : \n" << 1 << "\n";
  // parameters

  bool dep_cv_01 = sharedtype[0];
  bool dep_slope_01 = sharedtype[1];
  bool dep_var_01 = sharedtype[2];
  bool dep_cv_02 = sharedtype[3];
  bool dep_slope_02 = sharedtype[4];
  bool dep_var_02 = sharedtype[5];
  bool dep_cv_12 = sharedtype[6];
  bool dep_slope_12 = sharedtype[7];
  bool dep_var_12 = sharedtype[8];
  const std::string& hazard_baseline_01 = HB[0];
  const std::string& hazard_baseline_02 = HB[1];
  const std::string& hazard_baseline_12 = HB[2];
  bool left_trunc = HB[3];
  double shape_01 = W_G[0];
  double shape_02 = W_G[1];
  double shape_12 = W_G[2];
  double Gompertz_1_01 = W_G[3];
  double Gompertz_2_01 = W_G[4];
  double Gompertz_1_02 = W_G[5];
  double Gompertz_2_02 = W_G[6];
  double Gompertz_1_12 = W_G[7];
  double Gompertz_2_12 = W_G[8];


  double alpha_y_01 = alpha_y_slope_var[0];
  double alpha_y_02 = alpha_y_slope_var[1];
  double alpha_y_12 = alpha_y_slope_var[2];
  double alpha_slope_01 = alpha_y_slope_var[3];
  double alpha_slope_02 = alpha_y_slope_var[4];
  double alpha_slope_12 = alpha_y_slope_var[5];
  double alpha_var_01 = alpha_y_slope_var[6];
  double alpha_var_02 = alpha_y_slope_var[7];
  double alpha_var_12 = alpha_y_slope_var[8];
  arma::vec alpha_z_01 = alpha_z[0];
  arma::vec alpha_z_02 = alpha_z[1];
  arma::vec alpha_z_12 = alpha_z[2];
  //Rcout << "The value of v : \n" << 3 << "\n";
  arma::vec gamma_01 = gamma_B[0];
  arma::vec gamma_02 = gamma_B[1];
  arma::vec gamma_12 = gamma_B[2];

  arma::vec beta = fixed_par[0];
  arma::vec beta_slope = fixed_par[1];
  arma::vec omega = fixed_par[2];

  arma::vec y_new = longi[0];
  arma::vec offset = longi[1];
  int S = longi[2];
  int nb_pointsGK = longi[3];
  int nbCase3 = longi[4];

  arma::vec Time_T = list_Times[0];
  arma::vec Time_L = list_Times[1];
  arma::vec Time_L_T = list_Times[2];
  arma::vec Time_T0 = list_Times[3];
  arma::vec delta2 = list_Times[4];

  arma::vec sk_GK = ck[0];
  arma::vec wk = ck[1];
  arma::vec rep_wk = ck[2];

  arma::vec ll_glob(nbCase3,fill::ones);

  // Survival part
  ///// h
  arma::vec h_02_T_i_com(S,fill::ones);
  arma::vec etaBaseline_01_T_i_com(S,fill::zeros);
  arma::vec etaBaseline_02_T_i_com(S,fill::zeros);
  arma::vec etaBaseline_02_T0_i_com(S,fill::zeros);
  arma::vec etaBaseline_01_T0_i_com(S,fill::zeros);

  arma::vec h_12_T_i_com(S,fill::ones);
  arma::vec etaBaseline_12_T_i_com(S,fill::zeros);
  arma::vec etaBaseline_01_L_T_i_com(S,fill::zeros);
  arma::vec etaBaseline_01_0_LT_i_com(S,fill::zeros);
  arma::vec etaBaseline_02_0_LT_i_com(S,fill::zeros);
  arma::vec etaBaseline_12_0_LT_i_com(S,fill::zeros);



  for(int i_provCase3 =0; i_provCase3 < nbCase3; ++i_provCase3){

    arma::vec h_02_T_i=h_02_T_i_com;
    arma::vec etaBaseline_01_T_i=etaBaseline_01_T_i_com;
    arma::mat survLong_01_T_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_02_T_i=etaBaseline_02_T_i_com;
    arma::mat survLong_02_T_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_02_T0_i=etaBaseline_02_T0_i_com;
    arma::mat survLong_02_T0_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_01_T0_i=etaBaseline_01_T0_i_com;
    arma::mat survLong_01_T0_i(S,nb_pointsGK,fill::zeros);

    arma::vec h_12_T_i=h_12_T_i_com;
    arma::vec etaBaseline_12_T_i=etaBaseline_12_T_i_com;
    arma::mat survLong_12_T_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_01_L_T_i=etaBaseline_01_L_T_i_com;
    arma::mat survLong_01_L_T_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_01_0_LT_i=etaBaseline_01_0_LT_i_com;
    arma::mat survLong_01_0_LT_i(S,nb_pointsGK*nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_02_0_LT_i=etaBaseline_02_0_LT_i_com;
    arma::mat survLong_02_0_LT_i(S,nb_pointsGK*nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_12_0_LT_i=etaBaseline_12_0_LT_i_com;
    arma::mat survLong_12_0_LT_i(S,nb_pointsGK*nb_pointsGK,fill::zeros);

    arma::mat CV_T;
    arma::mat current_GK_T;
    arma::mat slope_T;
    arma::mat slope_GK_T;
    arma::mat current_GK_T0;
    arma::mat slope_GK_T0;

    arma::mat current_GK_L_T;
    arma::mat current_GK_0_LT;
    arma::mat slope_GK_L_T;
    arma::mat slope_GK_0_LT;

    arma::mat sigma_T;
    arma::mat sigma_GK_T;
    arma::mat sigma_GK_T0;
    arma::mat sigma_GK_L_T;
    arma::mat sigma_GK_0_LT;

    if(dep_cv_01 || dep_cv_02 || dep_cv_12){
      arma::rowvec X_T_i = X_T.row(i_provCase3);
      arma::vec U_T_i = U_T.col(i_provCase3);
      CV_T = arma::dot(beta, X_T_i) + b_y*U_T_i;
      arma::mat X_GK_T_i = X_GK_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat U_GK_T_i = U_GK_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      current_GK_T = arma::repmat(beta.t()*X_GK_T_i.t(),S,1)+b_y*U_GK_T_i.t();
      arma::mat X_GK_L_T_i = X_GK_L_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat U_GK_L_T_i = U_GK_L_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      current_GK_L_T = arma::repmat(beta.t()*X_GK_L_T_i.t(),S,1)+b_y*U_GK_L_T_i.t();
      arma::mat X_GK_0_LT_i = X_GK_0_LT.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      arma::mat U_GK_0_LT_i = U_GK_0_LT.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      current_GK_0_LT = arma::repmat(beta.t()*X_GK_0_LT_i.t(),S,1)+b_y*U_GK_0_LT_i.t();
      if(left_trunc){
        arma::mat X_GK_T0_i = X_GK_T0.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        arma::mat U_GK_T0_i = U_GK_T0.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        current_GK_T0 = arma::repmat(beta.t()*X_GK_T0_i.t(),S,1)+b_y*U_GK_T0_i.t();
      }
      if(dep_cv_01){
        survLong_01_T_i = survLong_01_T_i + alpha_y_01*current_GK_T;
        survLong_01_L_T_i = survLong_01_L_T_i + alpha_y_01*current_GK_L_T;
        survLong_01_0_LT_i = survLong_01_0_LT_i + alpha_y_01*current_GK_0_LT;
        if(left_trunc){
          survLong_01_T0_i = survLong_01_T0_i + alpha_y_01*current_GK_T0;
        }
      }
      if(dep_cv_02){
        h_02_T_i = h_02_T_i%exp(alpha_y_02*CV_T);
        survLong_02_T_i = survLong_02_T_i + alpha_y_02*current_GK_T;
        survLong_02_0_LT_i = survLong_02_0_LT_i + alpha_y_02*current_GK_0_LT;
        if(left_trunc){
          survLong_02_T0_i = survLong_02_T0_i + alpha_y_02*current_GK_T0;
        }
      }
      if(dep_cv_12){
        survLong_12_0_LT_i = survLong_12_0_LT_i + alpha_y_12*current_GK_0_LT;
        survLong_12_T_i = survLong_12_T_i + alpha_y_12*current_GK_T;
        h_12_T_i = h_12_T_i%exp(alpha_y_12*CV_T);
      }
    }

    if(dep_slope_01 || dep_slope_02 || dep_slope_12){
      arma::rowvec Xslope_T_i = Xslope_T.row(i_provCase3);
      arma::vec Uslope_T_i = Uslope_T.col(i_provCase3);
      slope_T = arma::dot(beta_slope, Xslope_T_i)+b_y_slope*Uslope_T_i;
      arma::mat Xslope_GK_T_i = Xslope_GK_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Uslope_GK_T_i = Uslope_GK_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      slope_GK_T = arma::repmat(beta_slope.t()*Xslope_GK_T_i.t(),S,1) + b_y_slope*Uslope_GK_T_i.t();
      arma::mat Xslope_GK_L_T_i = Xslope_GK_L_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Uslope_GK_L_T_i = Uslope_GK_L_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      slope_GK_L_T = arma::repmat(beta_slope.t()*Xslope_GK_L_T_i.t(),S,1) + b_y_slope*Uslope_GK_L_T_i.t();
      arma::mat Xslope_GK_0_LT_i = Xslope_GK_0_LT.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Uslope_GK_0_LT_i = Uslope_GK_0_LT.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      slope_GK_0_LT = arma::repmat(beta_slope.t()*Xslope_GK_0_LT_i.t(),S,1) + b_y_slope*Uslope_GK_0_LT_i.t();
      if(left_trunc){
        arma::mat Xslope_GK_T0_i = Xslope_GK_T0.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        arma::mat Uslope_GK_T0_i = Uslope_GK_T0.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        slope_GK_T0 = arma::repmat(beta_slope.t()*Xslope_GK_T0_i.t(),S,1)+b_y_slope*Uslope_GK_T0_i.t();
      }
      if(dep_slope_01){
        survLong_01_T_i = survLong_01_T_i + alpha_slope_01*slope_GK_T;
        survLong_01_L_T_i = survLong_01_L_T_i + alpha_slope_01*slope_GK_L_T;
        survLong_01_0_LT_i = survLong_01_0_LT_i + alpha_slope_01*slope_GK_0_LT;
        if(left_trunc){
          survLong_01_T0_i = survLong_01_T0_i + alpha_slope_01*slope_GK_T0;
        }
      }
      if(dep_slope_02){
        h_02_T_i = h_02_T_i%exp(alpha_slope_02*slope_T);
        survLong_02_T_i = survLong_02_T_i + alpha_slope_02*slope_GK_T;
        survLong_02_0_LT_i = survLong_02_0_LT_i + alpha_slope_02*slope_GK_0_LT;
        if(left_trunc){
          survLong_02_T0_i = survLong_02_T0_i + alpha_slope_02*slope_GK_T0;
        }
      }
      if(dep_slope_12){
        survLong_12_0_LT_i = survLong_12_0_LT_i + alpha_slope_12*slope_GK_0_LT;
        survLong_12_T_i = survLong_12_T_i + alpha_slope_12*slope_GK_T;
        h_12_T_i = h_12_T_i%exp(alpha_slope_12*slope_T);
      }
    }

    if(dep_var_01 || dep_var_02 || dep_var_12){
      arma::rowvec O_T_i = O_T.row(i_provCase3);
      arma::vec W_T_i = W_T.col(i_provCase3);
      sigma_T = exp(arma::dot(omega, O_T_i) + b_om*W_T_i);
      arma::mat O_GK_T_i = O_GK_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat W_GK_T_i = W_GK_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      sigma_GK_T = exp(arma::repmat(omega.t()*O_GK_T_i.t(),S,1)+b_om*W_GK_T_i.t());
      arma::mat O_GK_L_T_i = O_GK_L_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat W_GK_L_T_i = W_GK_L_T.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      sigma_GK_L_T = exp(arma::repmat(omega.t()*O_GK_L_T_i.t(),S,1)+b_om*W_GK_L_T_i.t());
      arma::mat O_GK_0_LT_i = O_GK_0_LT.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      arma::mat W_GK_0_LT_i = W_GK_0_LT.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      sigma_GK_0_LT = exp(arma::repmat(omega.t()*O_GK_0_LT_i.t(),S,1)+b_om*W_GK_0_LT_i.t());
      if(left_trunc){
        arma::mat O_GK_T0_i = O_GK_T0.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        arma::mat W_GK_T0_i = W_GK_T0.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        sigma_GK_T0 = exp(arma::repmat(omega.t()*O_GK_T0_i.t(),S,1)+b_om*W_GK_T0_i.t());
      }
      if(dep_var_01){
        survLong_01_T_i = survLong_01_T_i + alpha_var_01*sigma_GK_T;
        survLong_01_L_T_i = survLong_01_L_T_i + alpha_var_01*sigma_GK_L_T;
        survLong_01_0_LT_i = survLong_01_0_LT_i + alpha_var_01*sigma_GK_0_LT;
        if(left_trunc){
          survLong_01_T0_i = survLong_01_T0_i + alpha_var_01*sigma_GK_T0;
        }
      }
      if(dep_var_02){
        h_02_T_i = h_02_T_i%exp(alpha_var_02*sigma_T);
        survLong_02_T_i = survLong_02_T_i + alpha_var_02*sigma_GK_T;
        survLong_02_0_LT_i = survLong_02_0_LT_i + alpha_var_02*sigma_GK_0_LT;
        if(left_trunc){
          survLong_02_T0_i = survLong_02_T0_i + alpha_var_02*sigma_GK_T0;
        }
      }
      if(dep_var_12){
        survLong_12_0_LT_i = survLong_12_0_LT_i + alpha_var_12*sigma_GK_0_LT;
        survLong_12_T_i = survLong_12_T_i + alpha_var_12*sigma_GK_T;
        h_12_T_i = h_12_T_i%exp(alpha_var_12*sigma_T);
      }
    }

    ///// h0
    double Time_T_i = Time_T(i_provCase3);
    double Time_L_T_i = Time_L_T(i_provCase3);
    ///////// 1-2
    double h_0_12_T_i;
    arma::vec h_0_GK_12_T_i;
    arma::mat h_0_GK_12_0_LT_i;

    if(hazard_baseline_12 == "Exponential"){
      h_0_12_T_i = 1;
      h_0_GK_12_T_i = wk;
      h_0_GK_12_0_LT_i = rep_wk.t();
    }
    if(hazard_baseline_12 == "Weibull"){
      arma::vec st_T_i = st_T.col(i_provCase3);
      arma::mat st_0_LT_i = st_0_LT.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      h_0_12_T_i = shape_12*(pow(Time_T_i,(shape_12-1)));
      h_0_GK_12_T_i = shape_12*(pow(st_T_i,shape_12-1))%wk;
      h_0_GK_12_0_LT_i = shape_12*(pow(st_0_LT_i,shape_12-1));
      arma::rowvec  h_0_GK_12_0_LT_i2 = vectorise(h_0_GK_12_0_LT_i,1);
      h_0_GK_12_0_LT_i = h_0_GK_12_0_LT_i2%rep_wk.t();
    }
    if(hazard_baseline_12 == "Gompertz"){
      arma::vec st_T_i = st_T.col(i_provCase3);
      arma::mat st_0_LT_i = st_0_LT.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      h_0_12_T_i = Gompertz_1_12*exp(Gompertz_2_12*Time_T_i);
      h_0_GK_12_T_i = Gompertz_1_12*exp(Gompertz_2_12*st_T_i)%wk;
      h_0_GK_12_0_LT_i = Gompertz_1_12*exp(Gompertz_2_12*st_0_LT_i);
      arma::rowvec  h_0_GK_12_0_LT_i2 = vectorise(h_0_GK_12_0_LT_i,1);
      h_0_GK_12_0_LT_i = h_0_GK_12_0_LT_i2%rep_wk.t();
    }
    if(hazard_baseline_12 == "Splines"){
      arma::vec B_T_i_12 = B_T_12.col(i_provCase3);
      arma::mat Bs_T_i_12 = Bs_T_12.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Bs_0_LT_i_12 = Bs_0_LT_12.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      h_0_12_T_i = exp(arma::dot(gamma_12,B_T_i_12));
      h_0_GK_12_T_i = wk%exp(Bs_T_i_12*gamma_12);
      h_0_GK_12_0_LT_i = (exp(Bs_0_LT_i_12*gamma_12)%rep_wk).t();
    }
    double predsurv_12;
    arma::rowvec Z_12_i = Z_12.row(i_provCase3);
    if(Z_12_i.is_empty()){
      predsurv_12 = 0;
    }
    else{
      predsurv_12 = arma::dot(alpha_z_12, Z_12_i);
    }

    h_12_T_i = h_0_12_T_i*exp(predsurv_12)*h_12_T_i;

    etaBaseline_12_T_i = etaBaseline_12_T_i + predsurv_12;
    survLong_12_T_i = exp(survLong_12_T_i)*h_0_GK_12_T_i;
    arma::vec A_12_T_i;
    A_12_T_i = (exp(etaBaseline_12_T_i)%survLong_12_T_i*(Time_T_i/2)); //% : multiplication element par element

    etaBaseline_12_0_LT_i = exp(etaBaseline_12_0_LT_i + predsurv_12);
    survLong_12_0_LT_i = exp(survLong_12_0_LT_i);
    survLong_12_0_LT_i = survLong_12_0_LT_i%arma::repelem(h_0_GK_12_0_LT_i,S,1);
    arma::mat A_12_0_LT_i;
    A_12_0_LT_i = arma::repelem(etaBaseline_12_0_LT_i,1,nb_pointsGK*nb_pointsGK)%survLong_12_0_LT_i;

    ///////// 0-1
    arma::vec h_0_GK_01_L_T_i;
    arma::mat h_0_GK_01_0_LT_i;
    arma::vec h_0_GK_01_T0_i;
    arma::vec h_0_GK_01_T_i;
    if(hazard_baseline_01 == "Exponential"){
      h_0_GK_01_L_T_i = wk;
      h_0_GK_01_T_i = wk;
      h_0_GK_01_0_LT_i = rep_wk.t();
      if(left_trunc){
        h_0_GK_01_T0_i = wk;
      }
    }
    if(hazard_baseline_01 == "Weibull"){
      arma::vec st_T_i = st_T.col(i_provCase3);
      arma::vec st_L_T_i = st_L_T.col(i_provCase3);
      arma::mat st_0_LT_i = st_0_LT.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      h_0_GK_01_T_i = shape_01*(pow(st_T_i,shape_01-1))%wk;
      h_0_GK_01_L_T_i = shape_01*(pow(st_L_T_i,shape_01-1))%wk;
      h_0_GK_01_0_LT_i = shape_01*(pow(st_0_LT_i,shape_01-1));
      arma::rowvec  h_0_GK_01_0_LT_i2 = vectorise(h_0_GK_01_0_LT_i,1);
      h_0_GK_01_0_LT_i = h_0_GK_01_0_LT_i2%rep_wk.t();
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase3);
        h_0_GK_01_T0_i = shape_01*(pow(st_T0_i,shape_01-1))%wk;
      }
    }
    if(hazard_baseline_01 == "Gompertz"){
      arma::vec st_T_i = st_T.col(i_provCase3);
      arma::vec st_L_T_i = st_L_T.col(i_provCase3);
      arma::mat st_0_LT_i = st_0_LT.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      h_0_GK_01_T_i = Gompertz_1_01*exp(Gompertz_2_01*st_T_i)%wk;
      h_0_GK_01_L_T_i = Gompertz_1_01*exp(Gompertz_2_01*st_L_T_i)%wk;

      h_0_GK_01_0_LT_i = Gompertz_1_01*exp(Gompertz_2_01*st_0_LT_i);
      arma::rowvec  h_0_GK_01_0_LT_i2 = vectorise(h_0_GK_01_0_LT_i,1);
      h_0_GK_01_0_LT_i = h_0_GK_01_0_LT_i2%rep_wk.t();
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase3);
        h_0_GK_01_T0_i = Gompertz_1_01*exp(st_T0_i*Gompertz_2_01)%wk;
      }
    }
    if(hazard_baseline_01 == "Splines"){
      arma::mat Bs_T_i_01 = Bs_T_01.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Bs_L_T_i_01 = Bs_L_T_01.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Bs_0_LT_i_01 = Bs_0_LT_01.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      h_0_GK_01_T_i = wk%exp(Bs_T_i_01*gamma_01);
      h_0_GK_01_L_T_i = wk%exp(Bs_L_T_i_01*gamma_01);
      h_0_GK_01_0_LT_i = (exp(Bs_0_LT_i_01*gamma_01)%rep_wk).t();
      if(left_trunc){
        arma::mat Bs_T0_i_01 = Bs_T0_01.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        h_0_GK_01_T0_i = wk%exp(Bs_T0_i_01*gamma_01);
      }
    }
    double predsurv_01;
    arma::rowvec Z_01_i = Z_01.row(i_provCase3);
    if(Z_01_i.is_empty()){
      predsurv_01 = 0;
    }
    else{
      predsurv_01 = arma::dot(alpha_z_01, Z_01_i);
    }
    etaBaseline_01_L_T_i = exp(etaBaseline_01_L_T_i + predsurv_01);
    survLong_01_L_T_i = exp(survLong_01_L_T_i)%arma::repelem(h_0_GK_01_L_T_i.t(),S,1);
    arma::mat A_01_L_T_i;
    A_01_L_T_i = arma::repelem(etaBaseline_01_L_T_i,1,nb_pointsGK)%survLong_01_L_T_i*(Time_L_T_i/2);
    etaBaseline_01_T_i = exp(etaBaseline_01_T_i + predsurv_01);
    survLong_01_T_i = exp(survLong_01_T_i)*h_0_GK_01_T_i;
    arma::mat A_01_T_i;
    A_01_T_i = etaBaseline_01_T_i%survLong_01_T_i*(Time_T_i/2);

    etaBaseline_01_0_LT_i = exp(etaBaseline_01_0_LT_i + predsurv_01);
    survLong_01_0_LT_i = exp(survLong_01_0_LT_i);
    survLong_01_0_LT_i = survLong_01_0_LT_i%arma::repelem(h_0_GK_01_0_LT_i,S,1);
    arma::mat A_01_0_LT_i;
    A_01_0_LT_i = arma::repelem(etaBaseline_01_0_LT_i,1,nb_pointsGK*nb_pointsGK)%survLong_01_0_LT_i;

    arma::vec A_01_T0_i;
    if(left_trunc){
      double Time_T0_i = Time_T0(i_provCase3);
      etaBaseline_01_T0_i = etaBaseline_01_T0_i + predsurv_01;
      survLong_01_T0_i = exp(survLong_01_T0_i)*h_0_GK_01_T0_i;
      A_01_T0_i = (exp(etaBaseline_01_T0_i)%survLong_01_T0_i*(Time_T0_i/2));
    }

    ///////// 0-2
    arma::mat h_0_GK_02_0_LT_i;
    arma::vec h_0_GK_02_T0_i;
    arma::vec h_0_GK_02_T_i;
    double h_0_02_T_i;
    if(hazard_baseline_02 == "Exponential"){
      h_0_02_T_i = 1;
      h_0_GK_02_T_i = wk;
      h_0_GK_02_0_LT_i = rep_wk.t();
      if(left_trunc){
        h_0_GK_02_T0_i = wk;
      }
    }
    if(hazard_baseline_02 == "Weibull"){
      arma::vec st_T_i = st_T.col(i_provCase3);
      arma::mat st_0_LT_i = st_0_LT.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      h_0_GK_02_0_LT_i = shape_02*(pow(st_0_LT_i,shape_02-1));
      arma::rowvec  h_0_GK_02_0_LT_i2 = vectorise(h_0_GK_02_0_LT_i,1);
      h_0_GK_02_0_LT_i = h_0_GK_02_0_LT_i2%rep_wk.t();
      h_0_02_T_i = shape_02*(pow(Time_T_i,(shape_02-1)));
      h_0_GK_02_T_i = shape_02*(pow(st_T_i,shape_02-1))%wk;

      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase3);
        h_0_GK_02_T0_i = shape_02*(pow(st_T0_i,shape_02-1))%wk;
      }
    }
    if(hazard_baseline_02 == "Gompertz"){
      arma::vec st_T_i = st_T.col(i_provCase3);
      arma::mat st_0_LT_i = st_0_LT.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      h_0_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*Time_T_i);
      h_0_GK_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*st_T_i)%wk;
      h_0_GK_02_0_LT_i = Gompertz_1_02*exp(Gompertz_2_02*st_0_LT_i);
      arma::rowvec  h_0_GK_02_0_LT_i2 = vectorise(h_0_GK_02_0_LT_i,1);
      h_0_GK_02_0_LT_i = h_0_GK_02_0_LT_i2%rep_wk.t();
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase3);
        h_0_GK_02_T0_i = Gompertz_1_02*exp(st_T0_i*Gompertz_2_02)%wk;
      }
    }
    if(hazard_baseline_02 == "Splines"){
      arma::mat Bs_T_i_02= Bs_T_02.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
      arma::mat Bs_0_LT_i_02 = Bs_0_LT_02.rows((nb_pointsGK*nb_pointsGK*i_provCase3),(nb_pointsGK*nb_pointsGK*(i_provCase3+1)-1));
      arma::vec B_T_i_02 = B_T_02.col(i_provCase3);
      h_0_02_T_i = exp(arma::dot(gamma_02,B_T_i_02));
      h_0_GK_02_T_i = wk%exp(Bs_T_i_02*gamma_02);
      h_0_GK_02_0_LT_i = (exp(Bs_0_LT_i_02*gamma_02)%rep_wk).t();
      if(left_trunc){
        arma::mat Bs_T0_i_02 = Bs_T0_02.rows((nb_pointsGK*i_provCase3),(nb_pointsGK*(i_provCase3+1)-1));
        h_0_GK_02_T0_i = wk%exp(Bs_T0_i_02*gamma_02);
      }
    }
    double predsurv_02;
    arma::rowvec Z_02_i = Z_02.row(i_provCase3);
    if(Z_02_i.is_empty()){
      predsurv_02 = 0;
    }
    else{
      predsurv_02 = arma::dot(alpha_z_02, Z_02_i);
    }


    etaBaseline_02_T_i = exp(etaBaseline_02_T_i + predsurv_02);
    survLong_02_T_i = exp(survLong_02_T_i)*h_0_GK_02_T_i;
    arma::mat A_02_T_i;
    A_02_T_i = etaBaseline_02_T_i%survLong_02_T_i*(Time_T_i/2);

    etaBaseline_02_0_LT_i = exp(etaBaseline_02_0_LT_i + predsurv_02);
    survLong_02_0_LT_i = exp(survLong_02_0_LT_i);
    survLong_02_0_LT_i = survLong_02_0_LT_i%arma::repelem(h_0_GK_02_0_LT_i,S,1);
    arma::mat A_02_0_LT_i;
    A_02_0_LT_i = arma::repelem(etaBaseline_02_0_LT_i,1,nb_pointsGK*nb_pointsGK)%survLong_02_0_LT_i;

    h_02_T_i = h_0_02_T_i*exp(predsurv_02)*h_02_T_i;

    arma::vec A_02_T0_i;
    if(left_trunc){
      double Time_T0_i = Time_T0(i_provCase3);
      etaBaseline_02_T0_i = etaBaseline_02_T0_i + predsurv_02;
      survLong_02_T0_i = exp(survLong_02_T0_i)*h_0_GK_02_T0_i;
      A_02_T0_i = (exp(etaBaseline_02_T0_i)%survLong_02_T0_i*(Time_T0_i/2));
    }

    arma::mat A_0_LT = -A_01_0_LT_i - A_02_0_LT_i + A_12_0_LT_i;
    arma::mat A_0_LT_red(S,nb_pointsGK,fill::zeros);
    for(int column = 0; column < nb_pointsGK; ++column){
      A_0_LT_red.col(column) = sum(A_0_LT.cols(column*nb_pointsGK, (column+1)*nb_pointsGK-1),  1);
    }

    arma::vec ck = ((sk_GK+1)/4)*Time_L_T_i + Time_L(i_provCase3)/2;
    A_0_LT_red = A_0_LT_red%arma::repelem(ck.t(),S,1);
    arma::mat A_12_T_i_rep = arma::repelem(A_12_T_i,1,nb_pointsGK);
    int delta2_i = delta2(i_provCase3);

    arma::vec SurvTotCase3 = log(exp(-A_01_T_i - A_02_T_i)%(pow(h_02_T_i,delta2_i))+(pow(h_12_T_i,delta2_i))%sum(A_01_L_T_i%exp(A_0_LT_red-A_12_T_i_rep),1));
    /// Longitudinal part
    /// Longitudinal part

    arma::mat X_base_i = X_base.rows((offset(i_provCase3)-1),(offset(i_provCase3+1)-2));
    arma::mat U_base_i = U_base.rows((offset(i_provCase3)-1),(offset(i_provCase3+1)-2));
    arma::vec y_i = y_new.subvec(offset(i_provCase3)-1,offset(i_provCase3+1)-2);

    arma::mat O_base_i = O_base.rows((offset(i_provCase3)-1),(offset(i_provCase3+1)-2));
    arma::mat W_base_i = W_base.rows((offset(i_provCase3)-1),(offset(i_provCase3+1)-2));

    arma::vec f_Y_b_sigma(S,fill::zeros);
    arma::vec sigma_long;
    arma::vec CV;
    int n_rows_X = X_base_i.n_rows;
    for(int k=0; k<n_rows_X; k++){
      sigma_long = exp(dot(omega,O_base_i.row(k)) + b_om*W_base_i.row(k).t());
      CV = dot(beta,X_base_i.row(k)) + b_y*U_base_i.row(k).t();
      f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i(k)-CV)/sigma_long,2);
    }




    arma::vec log_dens_int;
    double Clogexp;
    double log_dens;
    log_dens_int = f_Y_b_sigma + SurvTotCase3;
    Clogexp = max(log_dens_int) - 500;
    log_dens_int = log_dens_int - Clogexp;
    log_dens = Clogexp + log(sum(exp(log_dens_int))) - log(S);
    //Rcout << "The value of log_dens : \n" << log_dens << "\n";
   // double den = 0;
   // if(left_trunc){
   //   //Rcout << "The value of A_01_T0_i : \n" << A_01_T0_i(0) << "\n";
   //   //Rcout << "The value of A_02_T0_i : \n" << A_02_T0_i(0) << "\n";
   //   den = log(sum(exp(-A_01_T0_i - A_02_T0_i)))-log(S);
   //   //Rcout << "The value of den : \n" << den << "\n";
//
   //   log_dens = log_dens - den;
   //   //Rcout << "The value of 3 : \n" << log_dens << "\n";
//
   // }

    ll_glob(i_provCase3) = log_dens;

  }

  // Rcout << "The value of v : \n" << A_01_L_i << "\n";
  return ll_glob;
}
