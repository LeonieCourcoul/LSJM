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

arma::vec log_llh_lsjm_covDepIDM_C2(arma::vec sharedtype, List HB, arma::vec W_G,
                                     arma::vec nb_points_integral,
                                     arma::vec alpha_y_slope_var, List alpha_b, List alpha_z, List gamma, List fixed_par,
                                     arma::mat b_y, arma::mat b_y_slope, arma::mat b_om, arma::vec wk,
                                     arma::mat Z_01, arma::mat Z_02, arma::mat X_T, arma::mat U_T,
                                     arma::mat Xslope_T, arma::mat Uslope_T, arma::mat X_GK_T, arma::mat U_GK_T, arma::mat Xslope_GK_T,
                                     arma::mat Uslope_GK_T,
                                     arma::mat X_GK_T0, arma::mat U_GK_T0, arma::mat Xslope_GK_T0, arma::mat Uslope_GK_T0,
                                     List list_Times,arma::mat st_T,  arma::mat st_T0,
                                     arma::mat B_T_02,
                                     arma::mat Bs_T_01, arma::mat Bs_T_02,
                                     arma::mat Bs_T0_01, arma::mat Bs_T0_02,
                                     arma::mat X_base, arma::mat U_base,  List longi, arma::mat O_T,
                                     arma::mat W_T, arma::mat O_GK_T, arma::mat W_GK_T, arma::mat O_GK_T0, arma::mat W_GK_T0, arma::mat O_base, arma::mat W_base
){
  // parameters
  bool dep_cv_01 = sharedtype[0];
  bool dep_slope_01 = sharedtype[1];
  bool dep_var_01 = sharedtype[2];
  bool dep_cv_02 = sharedtype[3];
  bool dep_slope_02 = sharedtype[4];
  bool dep_var_02 = sharedtype[5];
  bool dep_re_01 = sharedtype[9];
  bool dep_re_02 = sharedtype[10];

  const std::string& hazard_baseline_01 = HB[0];
  const std::string& hazard_baseline_02 = HB[1];

  bool left_trunc = HB[3];
  double shape_01 = W_G[0];
  double shape_02 = W_G[1];

  double Gompertz_1_01 = W_G[3];
  double Gompertz_2_01 = W_G[4];
  double Gompertz_1_02 = W_G[5];
  double Gompertz_2_02 = W_G[6];

  int S = nb_points_integral[0];
  int nb_pointsGK = nb_points_integral[1];
  int nbCase2 = nb_points_integral[2];

  double alpha_y_01 = alpha_y_slope_var[0];
  double alpha_y_02 = alpha_y_slope_var[1];

  double alpha_slope_01 = alpha_y_slope_var[3];
  double alpha_slope_02 = alpha_y_slope_var[4];

  double alpha_var_01 = alpha_y_slope_var[6];
  double alpha_var_02 = alpha_y_slope_var[7];

  arma::vec alpha_z_01 = alpha_z[0];
  arma::vec alpha_z_02 = alpha_z[1];

  arma::vec gamma_01 = gamma[0];
  arma::vec gamma_02 = gamma[1];

  arma::vec alpha_b_01 = alpha_b[0];
  arma::vec alpha_b_02 = alpha_b[1];

  arma::vec beta = fixed_par[0];
  arma::vec beta_slope = fixed_par[1];
  arma::vec omega = fixed_par[2];

  arma::vec y_new = longi[0];
  arma::vec offset = longi[1];

  arma::vec Time_T = list_Times[0];
  arma::vec Time_T0 = list_Times[1];
  arma::vec delta2 = list_Times[2];

  arma::vec ll_glob(nbCase2,fill::ones);

  // Survival part
  ///// h
  arma::vec h_02_T_i_com(S,fill::ones);
  arma::vec etaBaseline_01_T_i_com(S,fill::zeros);
  arma::vec etaBaseline_02_T_i_com(S,fill::zeros);
  arma::vec etaBaseline_02_T0_i_com(S,fill::zeros);
  arma::vec etaBaseline_01_T0_i_com(S,fill::zeros);


  for(int i_provCase2 =0; i_provCase2 < nbCase2; ++i_provCase2){

    arma::vec h_02_T_i=h_02_T_i_com;
    arma::vec etaBaseline_01_T_i=etaBaseline_01_T_i_com;
    arma::mat survLong_01_T_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_02_T_i=etaBaseline_02_T_i_com;
    arma::mat survLong_02_T_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_02_T0_i=etaBaseline_02_T0_i_com;
    arma::mat survLong_02_T0_i(S,nb_pointsGK,fill::zeros);
    arma::vec etaBaseline_01_T0_i=etaBaseline_01_T0_i_com;
    arma::mat survLong_01_T0_i(S,nb_pointsGK,fill::zeros);
    arma::mat CV_T;
    arma::mat current_GK_T;
    arma::mat slope_T;
    arma::mat slope_GK_T;
    arma::mat current_GK_T0;
    arma::mat slope_GK_T0;
    arma::mat sigma_T;
    arma::mat sigma_GK_T;
    arma::mat sigma_GK_T0;


    if(dep_re_01){
      survLong_01_T_i = survLong_01_T_i + arma::repmat(b_y*alpha_b_01,1,nb_pointsGK);
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + arma::repmat(b_y*alpha_b_01,1,nb_pointsGK);
      }
    }

    if(dep_re_02){
      survLong_02_T_i = survLong_02_T_i + arma::repmat(b_y*alpha_b_02,1,nb_pointsGK);
      h_02_T_i = h_02_T_i%exp(b_y*alpha_b_02);
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + arma::repmat(b_y*alpha_b_02,1,nb_pointsGK);
      }
    }


    if(dep_cv_01 || dep_cv_02){
      arma::rowvec X_T_i = X_T.row(i_provCase2);
      arma::vec U_T_i = U_T.col(i_provCase2);
      CV_T = arma::dot(beta, X_T_i) + b_y*U_T_i;
      arma::mat X_GK_T_i = X_GK_T.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      arma::mat U_GK_T_i = U_GK_T.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      current_GK_T = arma::repmat(beta.t()*X_GK_T_i.t(),S,1)+b_y*U_GK_T_i.t();
      if(left_trunc){
        arma::mat X_GK_T0_i = X_GK_T0.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        arma::mat U_GK_T0_i = U_GK_T0.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        current_GK_T0 = arma::repmat(beta.t()*X_GK_T0_i.t(),S,1)+b_y*U_GK_T0_i.t();
      }
      if(dep_cv_01){
        survLong_01_T_i = survLong_01_T_i + alpha_y_01*current_GK_T;
        if(left_trunc){
          survLong_01_T0_i = survLong_01_T0_i + alpha_y_01*current_GK_T0;
        }
      }
      if(dep_cv_02){
        h_02_T_i = h_02_T_i%exp(alpha_y_02*CV_T);
        survLong_02_T_i = survLong_02_T_i + alpha_y_02*current_GK_T;
        if(left_trunc){
          survLong_02_T0_i = survLong_02_T0_i + alpha_y_02*current_GK_T0;
        }
      }
    }

    if(dep_slope_01 || dep_slope_02){
      arma::rowvec Xslope_T_i = Xslope_T.row(i_provCase2);
      arma::vec Uslope_T_i = Uslope_T.col(i_provCase2);
      slope_T = arma::dot(beta_slope, Xslope_T_i)+b_y_slope*Uslope_T_i;
      arma::mat Xslope_GK_T_i = Xslope_GK_T.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      arma::mat Uslope_GK_T_i = Uslope_GK_T.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      slope_GK_T = arma::repmat(beta_slope.t()*Xslope_GK_T_i.t(),S,1) + b_y_slope*Uslope_GK_T_i.t();
      if(left_trunc){
        arma::mat Xslope_GK_T0_i = Xslope_GK_T0.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        arma::mat Uslope_GK_T0_i = Uslope_GK_T0.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        slope_GK_T0 = arma::repmat(beta_slope.t()*Xslope_GK_T0_i.t(),S,1)+b_y_slope*Uslope_GK_T0_i.t();
      }
      if(dep_slope_01){
        survLong_01_T_i = survLong_01_T_i + alpha_slope_01*slope_GK_T;
        if(left_trunc){
          survLong_01_T0_i = survLong_01_T0_i + alpha_slope_01*slope_GK_T0;
        }
      }
      if(dep_slope_02){
        h_02_T_i = h_02_T_i%exp(alpha_slope_02*slope_T);
        survLong_02_T_i = survLong_02_T_i + alpha_slope_02*slope_GK_T;
        if(left_trunc){
          survLong_02_T0_i = survLong_02_T0_i + alpha_slope_02*slope_GK_T0;
        }
      }
    }

    if(dep_var_01 || dep_var_02){
      arma::rowvec O_T_i = O_T.row(i_provCase2);
      arma::vec W_T_i = W_T.col(i_provCase2);
      sigma_T = exp(arma::dot(omega, O_T_i) + b_om*W_T_i);
      arma::mat O_GK_T_i = O_GK_T.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      arma::mat W_GK_T_i = W_GK_T.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      sigma_GK_T = exp(arma::repmat(omega.t()*O_GK_T_i.t(),S,1)+b_om*W_GK_T_i.t());
      if(left_trunc){
        arma::mat O_GK_T0_i = O_GK_T0.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        arma::mat W_GK_T0_i = W_GK_T0.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        sigma_GK_T0 = exp(arma::repmat(omega.t()*O_GK_T0_i.t(),S,1)+b_om*W_GK_T0_i.t());
      }
      if(dep_var_01){
        survLong_01_T_i = survLong_01_T_i + alpha_var_01*sigma_GK_T;
        if(left_trunc){
          survLong_01_T0_i = survLong_01_T0_i + alpha_var_01*sigma_GK_T0;
        }
      }
      if(dep_var_02){
        h_02_T_i = h_02_T_i%exp(alpha_var_02*sigma_T);
        survLong_02_T_i = survLong_02_T_i + alpha_var_02*sigma_GK_T;
        if(left_trunc){
          survLong_02_T0_i = survLong_02_T0_i + alpha_var_02*sigma_GK_T0;
        }
      }
    }

    ///// h0
    ///////// 0-1
    arma::vec h_0_GK_01_T_i;
    arma::vec h_0_GK_01_T0_i;
    if(hazard_baseline_01 == "Exponential"){
      h_0_GK_01_T_i = wk;
      if(left_trunc){
        h_0_GK_01_T0_i = wk;
      }
    }
    if(hazard_baseline_01 == "Weibull"){
      arma::vec st_T_i = st_T.col(i_provCase2);
      h_0_GK_01_T_i = shape_01*(pow(st_T_i,shape_01-1))%wk;
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase2);
        h_0_GK_01_T0_i = shape_01*(pow(st_T0_i,shape_01-1))%wk;
      }
    }
    if(hazard_baseline_01 == "Gompertz"){
      arma::vec st_T_i = st_T.col(i_provCase2);
      h_0_GK_01_T_i = Gompertz_1_01*exp(Gompertz_2_01*st_T_i)%wk;
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase2);
        h_0_GK_01_T0_i = Gompertz_1_01*exp(st_T0_i*Gompertz_2_01)%wk;
      }
    }
    if(hazard_baseline_01 == "Splines"){
      arma::mat Bs_T_i_01 = Bs_T_01.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      h_0_GK_01_T_i = wk%exp(Bs_T_i_01*gamma_01);
      if(left_trunc){
        arma::mat Bs_T0_i_01 = Bs_T0_01.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        h_0_GK_01_T0_i = wk%exp(Bs_T0_i_01*gamma_01);
      }
    }
    double predsurv_01;
    arma::rowvec Z_01_i = Z_01.row(i_provCase2);
    if(Z_01_i.is_empty()){
      predsurv_01 = 0;
    }
    else{
      predsurv_01 = arma::dot(alpha_z_01, Z_01_i);
    }

    etaBaseline_01_T_i = etaBaseline_01_T_i + predsurv_01;
    survLong_01_T_i = exp(survLong_01_T_i)*h_0_GK_01_T_i;
    arma::vec A_01_T_i;
    double Time_T_i = Time_T(i_provCase2);
    A_01_T_i = (exp(etaBaseline_01_T_i)%survLong_01_T_i*(Time_T_i/2));

    arma::vec A_01_T0_i;
    if(left_trunc){
      double Time_T0_i = Time_T0(i_provCase2);
      etaBaseline_01_T0_i = etaBaseline_01_T0_i + predsurv_01;
      survLong_01_T0_i = exp(survLong_01_T0_i)*h_0_GK_01_T0_i;
      A_01_T0_i = (exp(etaBaseline_01_T0_i)%survLong_01_T0_i*(Time_T0_i/2));
    }

    ///////// 0-2
    double h_0_02_T_i;
    arma::vec h_0_GK_02_T_i;
    arma::vec h_0_GK_02_T0_i;
    if(hazard_baseline_02 == "Exponential"){
      h_0_02_T_i = 1;
      h_0_GK_02_T_i = wk;
      if(left_trunc){
        h_0_GK_02_T0_i = wk;
      }
    }
    if(hazard_baseline_02 == "Weibull"){
      arma::vec st_T_i = st_T.col(i_provCase2);
      h_0_02_T_i = shape_02*(pow(Time_T_i,(shape_02-1)));
      h_0_GK_02_T_i = shape_02*(pow(st_T_i,shape_02-1))%wk;
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase2);
        h_0_GK_02_T0_i = shape_02*(pow(st_T0_i,shape_02-1))%wk;
      }
    }
    if(hazard_baseline_02 == "Gompertz"){
      arma::vec st_T_i = st_T.col(i_provCase2);
      h_0_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*Time_T_i);
      h_0_GK_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*st_T_i)%wk;
      if(left_trunc){
        arma::vec st_T0_i = st_T0.col(i_provCase2);
        h_0_GK_02_T0_i = Gompertz_1_02*exp(st_T0_i*Gompertz_2_02)%wk;
      }
    }
    if(hazard_baseline_02 == "Splines"){
      arma::vec B_T_i_02 = B_T_02.col(i_provCase2);
      arma::mat Bs_T_i_02 = Bs_T_02.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
      h_0_02_T_i = exp(arma::dot(gamma_02,B_T_i_02));
      h_0_GK_02_T_i = wk%exp(Bs_T_i_02*gamma_02);
      if(left_trunc){
        arma::mat Bs_T0_i_02 = Bs_T0_02.rows((nb_pointsGK*i_provCase2),(nb_pointsGK*(i_provCase2+1)-1));
        h_0_GK_02_T0_i = wk%exp(Bs_T0_i_02*gamma_02);
      }
    }
    double predsurv_02;
    arma::rowvec Z_02_i = Z_02.row(i_provCase2);
    if(Z_02_i.is_empty()){
      predsurv_02 = 0;
    }
    else{
      predsurv_02 = arma::dot(alpha_z_02, Z_02_i);
    }

    h_02_T_i = h_0_02_T_i*exp(predsurv_02)*h_02_T_i;

    etaBaseline_02_T_i = etaBaseline_02_T_i + predsurv_02;
    survLong_02_T_i = exp(survLong_02_T_i)*h_0_GK_02_T_i;
    arma::vec A_02_T_i;
    A_02_T_i = (exp(etaBaseline_02_T_i)%survLong_02_T_i*(Time_T_i/2));

    arma::vec A_02_T0_i;
    if(left_trunc){
      double Time_T0_i = Time_T0(i_provCase2);
      etaBaseline_02_T0_i = etaBaseline_02_T0_i + predsurv_02;
      survLong_02_T0_i = exp(survLong_02_T0_i)*h_0_GK_02_T0_i;
      A_02_T0_i = (exp(etaBaseline_02_T0_i)%survLong_02_T0_i*(Time_T0_i/2));
    }

    int delta2_i = delta2(i_provCase2);
    arma::vec SurvTotCase2 =  -A_01_T_i - A_02_T_i + log(pow(h_02_T_i,delta2_i));

    /// Longitudinal part

    arma::mat X_base_i = X_base.rows((offset(i_provCase2)-1),(offset(i_provCase2+1)-2));
    arma::mat U_base_i = U_base.rows((offset(i_provCase2)-1),(offset(i_provCase2+1)-2));
    arma::vec y_i = y_new.subvec(offset(i_provCase2)-1,offset(i_provCase2+1)-2);

    arma::mat O_base_i = O_base.rows((offset(i_provCase2)-1),(offset(i_provCase2+1)-2));
    arma::mat W_base_i = W_base.rows((offset(i_provCase2)-1),(offset(i_provCase2+1)-2));

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
    log_dens_int = f_Y_b_sigma + SurvTotCase2;
    Clogexp = max(log_dens_int) - 500;
    log_dens_int = log_dens_int - Clogexp;
    log_dens = Clogexp + log(sum(exp(log_dens_int))) - log(S);
 //   double den = 0;
   // if(left_trunc){
   //   den = log(sum(exp(-A_01_T0_i - A_02_T0_i)))-log(S);
   //   log_dens = log_dens - den;
   // }

    ll_glob(i_provCase2) = log_dens;

  }

  // Rcout << "The value of v : \n" << A_01_L_i << "\n";
  return ll_glob;
}

