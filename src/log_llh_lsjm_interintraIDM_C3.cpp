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

arma::vec log_llh_lsjm_interintraIDM_C3(arma::vec sharedtype, List HB, arma::vec Gompertz, arma::vec Weibull,
                                arma::vec nb_points_integral, arma::vec alpha_inter_intra,
                                arma::vec alpha_y_slope, List alpha_z, List gamma_B, arma::vec beta, arma::vec beta_slope,
                                arma::mat b_y, arma::mat b_y_slope, arma::vec wk, arma::vec rep_wk,  List sigma_inter_intra,
                                arma::vec delta2, arma::mat Z_01, arma::mat Z_02, arma::mat Z_12, arma::mat X_T, arma::mat U_T,
                                arma::mat Xslope_T, arma::mat Uslope_T, arma::mat X_GK_T, arma::mat U_GK_T, arma::mat Xslope_GK_T,
                                arma::mat Uslope_GK_T, arma::mat X_GK_L_T, arma::mat U_GK_L_T, arma::mat Xslope_GK_L_T, arma::mat Uslope_GK_L_T,
                                arma::mat X_GK_0_LT, arma::mat U_GK_0_LT, arma::mat Xslope_GK_0_LT, arma::mat Uslope_GK_0_LT,
                                arma::mat X_GK_T0, arma::mat U_GK_T0, arma::mat Xslope_GK_T0, arma::mat Uslope_GK_T0,
                                List list_Times, arma::mat st_T, arma::mat st_0_LT, arma::mat st_L_T, arma::mat st_T0,
                                arma::mat B_T_02, arma::mat B_T_12,
                                arma::mat Bs_T_01, arma::mat Bs_T_02, arma::mat Bs_T_12,
                                arma::mat Bs_0_LT_01, arma::mat Bs_0_LT_02, arma::mat Bs_0_LT_12,
                                arma::mat Bs_L_T_01,
                                arma::mat Bs_T0_01, arma::mat Bs_T0_02, bool left_trunc,
                                arma::vec len_visit, arma::mat X_base, arma::mat U_base,  arma::vec y_new, arma::vec offset_ID,
                                arma::vec offset, arma::vec offset_position, List ck
){

  //Rcout << "The value of v : \n" << 1 << "\n";
  // parameters
  bool dep_cv_01 = sharedtype[0];
  bool dep_slope_01 = sharedtype[1];
  bool dep_var_inter_01 = sharedtype[2];
  bool dep_var_intra_01 = sharedtype[3];
  bool dep_cv_02 = sharedtype[4];
  bool dep_slope_02 = sharedtype[5];
  bool dep_var_inter_02 = sharedtype[6];
  bool dep_var_intra_02 = sharedtype[7];
  bool dep_cv_12 = sharedtype[8];
  bool dep_slope_12 = sharedtype[9];
  bool dep_var_inter_12= sharedtype[10];
  bool dep_var_intra_12 = sharedtype[11];
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
  int S = nb_points_integral[0];
  int nb_pointsGK = nb_points_integral[1];
  double alpha_inter_01 = alpha_inter_intra[0];
  double alpha_inter_02 = alpha_inter_intra[1];
  double alpha_inter_12 = alpha_inter_intra[2];
  double alpha_intra_01 = alpha_inter_intra[3];
  double alpha_intra_02 = alpha_inter_intra[4];
  double alpha_intra_12 = alpha_inter_intra[5];
  double alpha_y_01 = alpha_y_slope[0];
  double alpha_y_02 = alpha_y_slope[1];
  double alpha_y_12 = alpha_y_slope[2];
  double alpha_slope_01 = alpha_y_slope[3];
  double alpha_slope_02 = alpha_y_slope[4];
  double alpha_slope_12 = alpha_y_slope[5];
  arma::vec alpha_z_01 = alpha_z[0];
  arma::vec alpha_z_02 = alpha_z[1];
  arma::vec alpha_z_12 = alpha_z[2];
  //Rcout << "The value of v : \n" << 3 << "\n";
  arma::vec gamma_01 = gamma_B[0];
  arma::vec gamma_02 = gamma_B[1];
  arma::vec gamma_12 = gamma_B[2];
  arma::vec sigma_inter = sigma_inter_intra[0];
  arma::vec sigma_intra = sigma_inter_intra[1];
  arma::vec sigma_long = sigma_inter_intra[2];
  arma::vec var_inter = sigma_inter_intra[3];
  arma::vec var_intra = sigma_inter_intra[4];
  arma::vec corr_intra_inter = sigma_inter_intra[5];
  arma::vec Time_T = list_Times[0];
  arma::vec Time_L_T = list_Times[1];
  arma::vec Time_T0 = list_Times[2];
  arma::vec sk_GK = ck[0];
  arma::vec Time_L = ck[1];
  int nbCase3 = ck[2];
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


  if(dep_var_inter_12){
    h_12_T_i_com = h_12_T_i_com%exp(alpha_inter_12*sigma_inter);
    etaBaseline_12_T_i_com = etaBaseline_12_T_i_com + alpha_inter_12*sigma_inter;
    etaBaseline_12_0_LT_i_com = etaBaseline_12_0_LT_i_com + alpha_inter_12*sigma_inter;
  }
  if(dep_var_inter_02){
    h_02_T_i_com = h_02_T_i_com%exp(alpha_inter_02*sigma_inter);
    etaBaseline_02_T_i_com = etaBaseline_02_T_i_com + alpha_inter_02*sigma_inter;
    etaBaseline_02_0_LT_i_com = etaBaseline_02_0_LT_i_com + alpha_inter_02*sigma_inter;
    if(left_trunc){
      etaBaseline_02_T0_i_com = etaBaseline_02_T0_i_com + alpha_inter_02*sigma_inter;
    }
  }
  if(dep_var_inter_01){
    etaBaseline_01_T_i_com = etaBaseline_01_T_i_com + alpha_inter_01*sigma_inter;
    etaBaseline_01_L_T_i_com = etaBaseline_01_L_T_i_com + alpha_inter_01*sigma_inter;
    etaBaseline_01_0_LT_i_com = etaBaseline_01_0_LT_i_com + alpha_inter_01*sigma_inter;
    if(left_trunc){
      etaBaseline_01_T0_i_com = etaBaseline_01_T0_i_com + alpha_inter_01*sigma_inter;
    }
  }
  if(dep_var_intra_12){
    h_12_T_i_com = h_12_T_i_com%exp(alpha_intra_12*sigma_intra);
    etaBaseline_12_T_i_com = etaBaseline_12_T_i_com + alpha_intra_12*sigma_intra;
    etaBaseline_12_0_LT_i_com = etaBaseline_12_0_LT_i_com + alpha_intra_12*sigma_intra;
  }
  if(dep_var_intra_02){
    h_02_T_i_com = h_02_T_i_com%exp(alpha_intra_02*sigma_intra);
    etaBaseline_02_T_i_com = etaBaseline_02_T_i_com + alpha_intra_02*sigma_intra;
    etaBaseline_02_0_LT_i_com = etaBaseline_02_0_LT_i_com + alpha_intra_02*sigma_intra;
    if(left_trunc){
      etaBaseline_02_T0_i_com = etaBaseline_02_T0_i_com + alpha_intra_02*sigma_intra;
    }
  }
  if(dep_var_intra_01){
    etaBaseline_01_T_i_com = etaBaseline_01_T_i_com + alpha_intra_01*sigma_intra;
    etaBaseline_01_L_T_i_com = etaBaseline_01_L_T_i_com + alpha_intra_01*sigma_intra;
    etaBaseline_01_0_LT_i_com = etaBaseline_01_0_LT_i_com + alpha_intra_01*sigma_intra;
    if(left_trunc){
      etaBaseline_01_T0_i_com = etaBaseline_01_T0_i_com + alpha_intra_01*sigma_intra;
    }
  }


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

    arma::mat X_base_i = X_base.rows((offset(i_provCase3)-1),(offset(i_provCase3+1)-2));
    arma::mat U_base_i = U_base.rows((offset(i_provCase3)-1),(offset(i_provCase3+1)-2));
    arma::vec y_i = y_new.subvec(offset(i_provCase3)-1,offset(i_provCase3+1)-2);
    double len_visit_i = len_visit(i_provCase3+1);
    arma::vec offset_ID_i = offset_ID.subvec(offset_position(i_provCase3)-1, offset_position(i_provCase3+1)-2);

    arma::vec f_Y_b_sigma(S,fill::zeros);
    arma::mat X_base_i_id_visit;
    arma::mat U_base_i_id_visit;
    arma::vec y_i_id_visit;
    arma::mat CV_long;
    for(int idvisit = 0; idvisit < len_visit_i; ++idvisit ){
      X_base_i_id_visit = X_base_i.row(offset_ID_i(idvisit)-1);
      U_base_i_id_visit = U_base_i.row(offset_ID_i(idvisit)-1);
      y_i_id_visit = y_i.subvec(offset_ID_i(idvisit)-1,offset_ID_i(idvisit+1)-2);
      int n_ij = y_i_id_visit.n_elem;
      CV_long = dot(beta,X_base_i_id_visit) + b_y*U_base_i_id_visit.t() ;
      if( n_ij == 1){
        f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i_id_visit(0)-CV_long)/sigma_long,2);
      }
      else{
        if(n_ij == 2){
          f_Y_b_sigma = f_Y_b_sigma + log(1/((pow(2*M_PI,n_ij/2))*sqrt(corr_intra_inter))) -
            (1/(2*corr_intra_inter))%((pow((y_i_id_visit(0)-CV_long),2)%sigma_long)-2*(var_inter%(y_i_id_visit(0)-CV_long)%(y_i_id_visit(1)-CV_long)) + (pow(y_i_id_visit(1)-CV_long,2)%(sigma_long)));
        }
        else{
          if(n_ij == 3){
            f_Y_b_sigma = f_Y_b_sigma + log(1/((pow(2*M_PI,n_ij/2.0))*var_intra%sqrt((3*var_inter+var_intra)))) -
              (1/(2*pow(var_intra,2)%(var_intra+3*var_inter)))%((corr_intra_inter%(pow(y_i_id_visit(0)-CV_long,2) + pow(y_i_id_visit(1)-CV_long,2) + pow(y_i_id_visit(2)-CV_long,2)))-
              2*var_inter%var_intra%((y_i_id_visit(0)-CV_long)%(y_i_id_visit(1)-CV_long) + (y_i_id_visit(0)-CV_long)%(y_i_id_visit(2)-CV_long) + (y_i_id_visit(1)-CV_long)%(y_i_id_visit(2)-CV_long)));
          }
          else{
            arma::vec somme1(S,fill::zeros);
            arma::vec somme2(S,fill::zeros);
            for(int k_somme = 0; k_somme < n_ij; ++k_somme){
              somme1 = somme1 + pow(y_i_id_visit(k_somme)-CV_long,2);
              if(k_somme != n_ij){
                for(int l_somme = k_somme+1; l_somme < n_ij; ++l_somme){
                  somme2 = somme2 + (y_i_id_visit(k_somme)-CV_long)%(y_i_id_visit(l_somme)-CV_long);
                }
              }
            }
            f_Y_b_sigma = f_Y_b_sigma + log(1/((pow(2*M_PI,n_ij/2.0))*sqrt(pow(var_intra,(n_ij-1))%(var_intra+n_ij*var_inter)))) -
              (1/(2*pow(var_intra,(n_ij-1))%(var_intra+n_ij*var_inter)))%(pow(var_intra,(n_ij-2))%(var_intra+(n_ij-1)*var_inter)%somme1 -
              2*var_inter%pow(var_intra,(n_ij-2))%somme2);
          }
        }
      }

    }





    arma::vec log_dens_int;
    double Clogexp;
    double log_dens;
    log_dens_int = f_Y_b_sigma + SurvTotCase3;
    Clogexp = max(log_dens_int) - 500;
    log_dens_int = log_dens_int - Clogexp;
    log_dens = Clogexp + log(sum(exp(log_dens_int))) - log(S);
    //Rcout << "The value of log_dens : \n" << log_dens << "\n";
//    double den = 0;
  //  if(left_trunc){
  //    //Rcout << "The value of A_01_T0_i : \n" << A_01_T0_i(0) << "\n";
  //    //Rcout << "The value of A_02_T0_i : \n" << A_02_T0_i(0) << "\n";
  //    den = log(sum(exp(-A_01_T0_i - A_02_T0_i)))-log(S);
  //    //Rcout << "The value of den : \n" << den << "\n";
//
  //    log_dens = log_dens - den;
  //    //Rcout << "The value of 3 : \n" << log_dens << "\n";
//
  //  }

    ll_glob(i_provCase3) = log_dens;

  }

  // Rcout << "The value of v : \n" << A_01_L_i << "\n";
  return ll_glob;
}
