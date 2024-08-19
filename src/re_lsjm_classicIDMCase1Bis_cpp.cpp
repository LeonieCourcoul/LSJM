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

double re_lsjm_classicIDMCase1Bis_cpp(arma::vec sharedtype, List HB, arma::vec Gompertz, arma::vec Weibull,
                                         double nb_pointsGK,
                                         arma::vec alpha_y_slope, List alpha_z, List gamma, arma::vec beta, arma::vec beta_slope,
                                         arma::mat b_y, arma::mat b_y_slope, arma::vec wk, double sigma_epsilon,
                                         int delta2_i, arma::rowvec Z_01_i, arma::rowvec Z_02_i, arma::rowvec Z_12_i, arma::rowvec X_T_i, arma::rowvec U_T_i,
                                         arma::rowvec Xslope_T_i, arma::rowvec Uslope_T_i, arma::mat X_GK_T_i, arma::mat U_GK_T_i, arma::mat Xslope_GK_T_i,
                                         arma::mat Uslope_GK_T_i, arma::rowvec X_L_i, arma::rowvec U_L_i,
                                         arma::rowvec Xslope_L_i, arma::rowvec Uslope_L_i, arma::mat X_GK_L_i, arma::mat U_GK_L_i, arma::mat Xslope_GK_L_i,
                                         arma::mat Uslope_GK_L_i,
                                         arma::mat X_GK_T0_i, arma::mat U_GK_T0_i, arma::mat Xslope_GK_T0_i, arma::mat Uslope_GK_T0_i,
                                         double Time_T_i, double Time_L_i, double Time_T0_i,arma::vec st_T_i, arma::vec st_L_i, arma::vec st_T0_i,
                                         arma::vec B_T_i_12,
                                         arma::mat Bs_T_i_12,
                                         arma::vec B_L_i_01,
                                         arma::mat Bs_L_i_01, arma::mat Bs_L_i_02, arma::mat Bs_L_i_12,
                                         arma::mat Bs_T0_i_01, arma::mat Bs_T0_i_02, bool left_trunc,
                                         arma::mat X_base_i, arma::mat U_base_i,  arma::vec y_i
){


  // parameters
  bool dep_cv_01 = sharedtype[0];
  bool dep_slope_01 = sharedtype[1];

  bool dep_cv_02 = sharedtype[2];
  bool dep_slope_02 = sharedtype[3];

  bool dep_cv_12 = sharedtype[4];
  bool dep_slope_12 = sharedtype[5];

  const std::string& hazard_baseline_01 = HB[0];
  const std::string& hazard_baseline_02 = HB[1];
  const std::string& hazard_baseline_12 = HB[2];
  double Gompertz_1_01 = Gompertz[0];
  double Gompertz_2_01 = Gompertz[1];
  double Gompertz_1_02 = Gompertz[2];
  double Gompertz_2_02 = Gompertz[3];
  double Gompertz_1_12 = Gompertz[4];
  double Gompertz_2_12 = Gompertz[5];
  double shape_01 = Weibull[0];
  double shape_02 = Weibull[1];
  double shape_12 = Weibull[2];
  double alpha_y_01 = alpha_y_slope[0];
  double alpha_y_02 = alpha_y_slope[1];
  double alpha_y_12 = alpha_y_slope[2];
  double alpha_slope_01 = alpha_y_slope[3];
  double alpha_slope_02 = alpha_y_slope[4];
  double alpha_slope_12 = alpha_y_slope[5];
  arma::vec alpha_z_01 = alpha_z[0];
  arma::vec alpha_z_02 = alpha_z[1];
  arma::vec alpha_z_12 = alpha_z[2];

  arma::vec gamma_01 = gamma[0];
  arma::vec gamma_02 = gamma[1];
  arma::vec gamma_12 = gamma[2];


  // Survival part
  ///// h
  int S = 1;
  arma::vec h_12_T_i(S,fill::ones);
  arma::vec h_01_L_i(S,fill::ones);
  arma::vec etaBaseline_12_T_i(S,fill::zeros);
  arma::mat survLong_12_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_01_L_i(S,fill::zeros);
  arma::mat survLong_01_L_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_L_i(S,fill::zeros);
  arma::mat survLong_02_L_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_12_L_i(S,fill::zeros);
  arma::mat survLong_12_L_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_T0_i(S,fill::zeros);
  arma::mat survLong_02_T0_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_01_T0_i(S,fill::zeros);
  arma::mat survLong_01_T0_i(S,nb_pointsGK,fill::zeros);
  arma::mat CV_T;
  arma::mat current_GK_T;
  arma::mat CV_L;
  arma::mat current_GK_L;
  arma::mat slope_T;
  arma::mat slope_GK_T;
  arma::mat slope_L;
  arma::mat slope_GK_L;
  arma::mat current_GK_T0;
  arma::mat slope_GK_T0;

  if(dep_cv_01 || dep_cv_02 || dep_cv_12){
    CV_T = arma::dot(beta, X_T_i) + U_T_i*b_y;
    current_GK_T = X_GK_T_i*beta+U_GK_T_i*b_y;
    CV_L = arma::dot(beta, X_L_i) + U_L_i*b_y;
    current_GK_L = X_GK_L_i*beta+U_GK_L_i*b_y;;
    if(left_trunc){
      current_GK_T0 = X_GK_T0_i*beta+U_GK_T0_i*b_y;
    }
    if(dep_cv_01){
      h_01_L_i = h_01_L_i%exp(alpha_y_01*CV_L);
      survLong_01_L_i = survLong_01_L_i + alpha_y_01*current_GK_L.t();
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + alpha_y_01*current_GK_T0.t();
      }
    }
    if(dep_cv_02){
      survLong_02_L_i = survLong_02_L_i + alpha_y_02*current_GK_L.t();
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + alpha_y_02*current_GK_T0.t();
      }
    }
    if(dep_cv_12){
      h_12_T_i = h_12_T_i%exp(alpha_y_12*CV_T);
      survLong_12_L_i = survLong_12_L_i + alpha_y_12*current_GK_L.t();
      survLong_12_T_i = survLong_12_T_i + alpha_y_12*current_GK_T.t();
    }
  }
  if(dep_slope_01 || dep_slope_02 || dep_slope_12){
    slope_T = arma::dot(beta_slope, Xslope_T_i)+Uslope_T_i*b_y_slope;
    slope_GK_T = Xslope_GK_T_i*beta_slope+Uslope_GK_T_i*b_y_slope;
    slope_L = arma::dot(beta_slope, Xslope_L_i)+Uslope_L_i*b_y_slope;
    slope_GK_L = Xslope_GK_L_i*beta_slope+Uslope_GK_L_i*b_y_slope;
    if(left_trunc){
      slope_GK_T0 = Xslope_GK_T0_i*beta_slope+Uslope_GK_T0_i*b_y_slope;
    }
    if(dep_slope_01){
      h_01_L_i = h_01_L_i%exp(alpha_slope_01*slope_L);
      survLong_01_L_i = survLong_01_L_i + alpha_slope_01*slope_GK_L.t();
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + alpha_slope_01*slope_GK_T0.t();
      }
    }
    if(dep_slope_02){
      survLong_02_L_i = survLong_02_L_i + alpha_slope_02*slope_GK_L.t();
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + alpha_slope_02*slope_GK_T0.t();
      }
    }
    if(dep_slope_12){
      h_12_T_i = h_12_T_i%exp(alpha_slope_12*slope_T);
      survLong_12_L_i = survLong_12_L_i + alpha_slope_12*slope_GK_L.t();
      survLong_12_T_i = survLong_12_T_i + alpha_slope_12*slope_GK_T.t();
    }
  }

  ///// h0
  ///////// 1-2
  double h_0_12_T_i;
  arma::vec h_0_GK_12_T_i;
  arma::vec h_0_GK_12_L_i;
  if(hazard_baseline_12 == "Exponential"){
    h_0_12_T_i = 1;
    h_0_GK_12_T_i = wk;
    h_0_GK_12_L_i = wk;
  }
  if(hazard_baseline_12 == "Weibull"){
    h_0_12_T_i = shape_12*(pow(Time_T_i,(shape_12-1)));
    h_0_GK_12_T_i = shape_12*(pow(st_T_i,shape_12-1))%wk;
    h_0_GK_12_L_i = shape_12*(pow(st_L_i,shape_12-1))%wk;
  }
  if(hazard_baseline_12 == "Gompertz"){
    h_0_12_T_i = Gompertz_1_12*exp(Gompertz_2_12*Time_T_i);
    h_0_GK_12_T_i = Gompertz_1_12*exp(Gompertz_2_12*st_T_i)%wk;
    h_0_GK_12_L_i = Gompertz_1_12*exp(Gompertz_2_12*st_L_i)%wk;
  }
  if(hazard_baseline_12 == "Splines"){
    h_0_12_T_i = exp(arma::dot(gamma_12,B_T_i_12));
    h_0_GK_12_T_i = wk%exp(Bs_T_i_12*gamma_12);
    h_0_GK_12_L_i = wk%exp(Bs_L_i_12*gamma_12);
  }
  double predsurv_12;
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

  etaBaseline_12_L_i = etaBaseline_12_L_i + predsurv_12;
  survLong_12_L_i = exp(survLong_12_L_i)*h_0_GK_12_L_i;
  arma::vec A_12_L_i;
  A_12_L_i = (exp(etaBaseline_12_L_i)%survLong_12_L_i*(Time_L_i/2));

  ///////// 0-1
  double h_0_01_L_i;
  arma::vec h_0_GK_01_L_i;
  arma::vec h_0_GK_01_T0_i;
  if(hazard_baseline_01 == "Exponential"){
    h_0_01_L_i = 1;
    h_0_GK_01_L_i = wk;
    if(left_trunc){
      h_0_GK_01_T0_i = wk;
    }
  }
  if(hazard_baseline_01 == "Weibull"){
    h_0_01_L_i = shape_01*(pow(Time_L_i,(shape_01-1)));
    h_0_GK_01_L_i = shape_01*(pow(st_L_i,shape_01-1))%wk;
    if(left_trunc){
      h_0_GK_01_T0_i = shape_01*(pow(st_T0_i,shape_01-1))%wk;
    }
  }
  if(hazard_baseline_01 == "Gompertz"){
    h_0_01_L_i = Gompertz_1_01*exp(Gompertz_2_01*Time_L_i);
    h_0_GK_01_L_i = Gompertz_1_01*exp(Gompertz_2_01*st_L_i)%wk;
    if(left_trunc){
      h_0_GK_01_T0_i = Gompertz_1_01*exp(st_T0_i*Gompertz_2_01)%wk;
    }
  }
  if(hazard_baseline_01 == "Splines"){
    h_0_01_L_i = exp(arma::dot(gamma_01,B_L_i_01));
    h_0_GK_01_L_i = wk%exp(Bs_L_i_01*gamma_01);
    if(left_trunc){
      h_0_GK_01_T0_i = wk%exp(Bs_T0_i_01*gamma_01);
    }
  }
  double predsurv_01;
  if(Z_01_i.is_empty()){
    predsurv_01 = 0;
  }
  else{
    predsurv_01 = arma::dot(alpha_z_01, Z_01_i);
  }

  h_01_L_i = h_0_01_L_i*exp(predsurv_01)*h_01_L_i;

  etaBaseline_01_L_i = etaBaseline_01_L_i + predsurv_01;
  survLong_01_L_i = exp(survLong_01_L_i)*h_0_GK_01_L_i;
  arma::vec A_01_L_i;
  A_01_L_i = (exp(etaBaseline_01_L_i)%survLong_01_L_i*(Time_L_i/2));

  arma::vec A_01_T0_i;
  if(left_trunc){
    etaBaseline_01_T0_i = etaBaseline_01_T0_i + predsurv_01;
    survLong_01_T0_i = exp(survLong_01_T0_i)*h_0_GK_01_T0_i;
    A_01_T0_i = (exp(etaBaseline_01_T0_i)%survLong_01_T0_i*(Time_T0_i/2));
  }
  ///////// 0-2
  arma::vec h_0_GK_02_L_i;
  arma::vec h_0_GK_02_T0_i;
  if(hazard_baseline_02 == "Exponential"){
    h_0_GK_02_L_i = wk;
    if(left_trunc){
      h_0_GK_02_T0_i = wk;
    }
  }
  if(hazard_baseline_02 == "Weibull"){
    h_0_GK_02_L_i = shape_02*(pow(st_L_i,shape_02-1))%wk;
    if(left_trunc){
      h_0_GK_02_T0_i = shape_02*(pow(st_T0_i,shape_02-1))%wk;
    }
  }
  if(hazard_baseline_02 == "Gompertz"){
    h_0_GK_02_L_i = Gompertz_1_02*exp(Gompertz_2_02*st_L_i)%wk;
    if(left_trunc){
      h_0_GK_02_T0_i = Gompertz_1_02*exp(st_T0_i*Gompertz_2_02)%wk;
    }
  }
  if(hazard_baseline_02 == "Splines"){
    h_0_GK_02_L_i = wk%exp(Bs_L_i_02*gamma_02);
    if(left_trunc){
      h_0_GK_02_T0_i = wk%exp(Bs_T0_i_02*gamma_02);
    }
  }
  double predsurv_02;
  if(Z_02_i.is_empty()){
    predsurv_02 = 0;
  }
  else{
    predsurv_02 = arma::dot(alpha_z_02, Z_02_i);
  }

  etaBaseline_02_L_i = etaBaseline_02_L_i + predsurv_02;
  survLong_02_L_i = exp(survLong_02_L_i)*h_0_GK_02_L_i;
  arma::vec A_02_L_i;
  A_02_L_i = (exp(etaBaseline_02_L_i)%survLong_02_L_i*(Time_L_i/2));
  arma::vec A_02_T0_i;
  if(left_trunc){
    etaBaseline_02_T0_i = etaBaseline_02_T0_i + predsurv_02;
    survLong_02_T0_i = exp(survLong_02_T0_i)*h_0_GK_02_T0_i;
    A_02_T0_i = (exp(etaBaseline_02_T0_i)%survLong_02_T0_i*(Time_T0_i/2));
  }

  arma::vec SurvTotCase1bis =  -A_01_L_i - A_02_L_i + log(h_01_L_i) - A_12_T_i + A_12_L_i + log(pow(h_12_T_i,delta2_i));

  // Rcout << "The value of v : \n" << 8 << "\n";
  arma::vec f_Y_b_sigma(1,fill::zeros);
  arma::vec CV;
  int n_rows_X = X_base_i.n_rows;
  double sigma_long;
  sigma_long = sigma_epsilon;
  for(int k=0; k<n_rows_X; k++){
    CV = dot(beta,X_base_i.row(k)) + b_y*U_base_i.row(k).t();
    f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i(k)-CV)/sigma_long,2);
  }


  // Rcout << "The value of v : \n" << 9 << "\n";

  arma::vec log_dens_int;
  double Clogexp;
  double log_dens;
  log_dens_int = f_Y_b_sigma + SurvTotCase1bis;
  Clogexp = max(log_dens_int) - 500;
  log_dens_int = log_dens_int - Clogexp;
  log_dens = Clogexp + log(sum(exp(log_dens_int)));
  double den = 0;
  if(left_trunc){
    den = log(sum(exp(-A_01_T0_i - A_02_T0_i)));
    log_dens = log_dens - den;
  }

  return log_dens;
}



