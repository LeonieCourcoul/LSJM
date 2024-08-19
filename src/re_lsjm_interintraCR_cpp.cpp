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

double re_lsjm_interintraCR_cpp(arma::vec sharedtype, List HB, arma::vec Gompertz, arma::vec Weibull,
                             double nb_pointsGK ,
                             arma::vec alpha_y_slope, arma::vec alpha_inter_intra, List alpha_z, List gamma, arma::vec beta, arma::vec beta_slope,
                             arma::mat b_y, arma::mat b_y_slope, arma::vec wk, List sigma_inter_intra,
                             int delta1_i, int delta2_i, arma::rowvec Z_01_i, arma::rowvec Z_02_i, arma::rowvec X_T_i, arma::rowvec U_T_i,
                             arma::rowvec Xslope_T_i, arma::rowvec Uslope_T_i, arma::mat X_GK_T_i, arma::mat U_GK_T_i, arma::mat Xslope_GK_T_i,
                             arma::mat Uslope_GK_T_i,
                             arma::mat X_GK_T0_i, arma::mat U_GK_T0_i, arma::mat Xslope_GK_T0_i, arma::mat Uslope_GK_T0_i,
                             double Time_T_i,  double Time_T0_i,arma::vec st_T_i,  arma::vec st_T0_i,
                             arma::vec B_T_i_01, arma::vec B_T_i_02,
                             arma::mat Bs_T_i_01, arma::mat Bs_T_i_02,
                             arma::mat Bs_T0_i_01, arma::mat Bs_T0_i_02, bool left_trunc,
                             int len_visit_i, arma::mat X_base_i, arma::mat U_base_i,  arma::vec y_i, arma::vec offset_ID_i
){
  // parameters
  bool dep_cv_01 = sharedtype[0];
  bool dep_slope_01 = sharedtype[1];
  bool dep_var_inter_01 = sharedtype[2];
  bool dep_var_intra_01 = sharedtype[3];

  bool dep_cv_02 = sharedtype[4];
  bool dep_slope_02 = sharedtype[5];
  bool dep_var_inter_02 = sharedtype[6];
  bool dep_var_intra_02 = sharedtype[7];

  const std::string& hazard_baseline_01 = HB[0];
  const std::string& hazard_baseline_02 = HB[1];
  double Gompertz_1_01 = Gompertz[0];
  double Gompertz_2_01 = Gompertz[1];
  double Gompertz_1_02 = Gompertz[2];
  double Gompertz_2_02 = Gompertz[3];
  double shape_01 = Weibull[0];
  double shape_02 = Weibull[1];

  double alpha_y_01 = alpha_y_slope[0];
  double alpha_y_02 = alpha_y_slope[1];
  double alpha_slope_01 = alpha_y_slope[2];
  double alpha_slope_02 = alpha_y_slope[3];
  double alpha_inter_01 = alpha_inter_intra[0];
  double alpha_inter_02 = alpha_inter_intra[1];
  double alpha_intra_01 = alpha_inter_intra[2];
  double alpha_intra_02 = alpha_inter_intra[3];
  arma::vec alpha_z_01 = alpha_z[0];
  arma::vec alpha_z_02 = alpha_z[1];
  arma::vec gamma_01 = gamma[0];
  arma::vec gamma_02 = gamma[1];

  arma::vec sigma_inter = sigma_inter_intra[0];
  arma::vec sigma_intra = sigma_inter_intra[1];
  arma::vec sigma_long = sigma_inter_intra[2];
  arma::vec var_inter = sigma_inter_intra[3];
  arma::vec var_intra = sigma_inter_intra[4];
  arma::vec corr_intra_inter = sigma_inter_intra[5];



  // Survival part
  ///// h
  int S = 1;
  arma::vec h_02_T_i(S,fill::ones);
  arma::vec h_01_T_i(S,fill::ones);
  arma::vec etaBaseline_01_T_i(S,fill::zeros);
  arma::mat survLong_01_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_T_i(S,fill::zeros);
  arma::mat survLong_02_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_T0_i(S,fill::zeros);
  arma::mat survLong_02_T0_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_01_T0_i(S,fill::zeros);
  arma::mat survLong_01_T0_i(S,nb_pointsGK,fill::zeros);
  arma::mat CV_T;
  arma::mat current_GK_T;
  arma::mat slope_T;
  arma::mat slope_GK_T;
  arma::mat current_GK_T0;
  arma::mat slope_GK_T0;


  if(dep_var_inter_01){
    h_01_T_i = h_01_T_i%exp(alpha_inter_01*sigma_inter);
    etaBaseline_01_T_i = etaBaseline_01_T_i + alpha_inter_01*sigma_inter;
    if(left_trunc){
      etaBaseline_01_T0_i = etaBaseline_01_T0_i + alpha_inter_01*sigma_inter;
    }
  }
  if(dep_var_inter_02){
    h_02_T_i = h_02_T_i%exp(alpha_inter_02*sigma_inter);
    etaBaseline_02_T_i = etaBaseline_02_T_i + alpha_inter_02*sigma_inter;
    if(left_trunc){
      etaBaseline_02_T0_i = etaBaseline_02_T0_i + alpha_inter_02*sigma_inter;
    }
  }

  if(dep_var_intra_01){
    h_01_T_i = h_01_T_i%exp(alpha_intra_01*sigma_intra);
    etaBaseline_01_T_i = etaBaseline_01_T_i + alpha_intra_01*sigma_intra;
    if(left_trunc){
      etaBaseline_01_T0_i = etaBaseline_01_T0_i + alpha_intra_01*sigma_intra;
    }
  }
  if(dep_var_intra_02){
    h_02_T_i = h_02_T_i%exp(alpha_intra_02*sigma_intra);
    etaBaseline_02_T_i = etaBaseline_02_T_i + alpha_intra_02*sigma_intra;
    if(left_trunc){
      etaBaseline_02_T0_i = etaBaseline_02_T0_i + alpha_intra_02*sigma_intra;
    }
  }

  if(dep_cv_01 || dep_cv_02){
    CV_T = arma::dot(beta, X_T_i) + U_T_i*b_y;
    current_GK_T = X_GK_T_i*beta+U_GK_T_i*b_y;
    if(left_trunc){
      current_GK_T0 = X_GK_T0_i*beta+U_GK_T0_i*b_y;
    }
    if(dep_cv_01){
      h_01_T_i = h_01_T_i%exp(alpha_y_01*CV_T);
      survLong_01_T_i = survLong_01_T_i + alpha_y_01*current_GK_T.t();
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + alpha_y_01*current_GK_T0.t();
      }
    }
    if(dep_cv_02){
      h_02_T_i = h_02_T_i%exp(alpha_y_02*CV_T);
      survLong_02_T_i = survLong_02_T_i + alpha_y_02*current_GK_T.t();
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + alpha_y_02*current_GK_T0.t();
      }
    }
  }
  if(dep_slope_01 || dep_slope_02){
    slope_T = arma::dot(beta_slope, Xslope_T_i)+Uslope_T_i*b_y_slope;
    slope_GK_T = Xslope_GK_T_i*beta_slope+Uslope_GK_T_i*b_y_slope;
    if(left_trunc){
      slope_GK_T0 = Xslope_GK_T0_i*beta_slope+Uslope_GK_T0_i*b_y_slope;
    }
    if(dep_slope_01){
      h_01_T_i = h_01_T_i%exp(alpha_slope_01*slope_T);
      survLong_01_T_i = survLong_01_T_i + alpha_slope_01*slope_GK_T.t();
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + alpha_slope_01*slope_GK_T0.t();
      }
    }
    if(dep_slope_02){
      h_02_T_i = h_02_T_i%exp(alpha_slope_02*slope_T);
      survLong_02_T_i = survLong_02_T_i + alpha_slope_02*slope_GK_T.t();
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + alpha_slope_02*slope_GK_T0.t();
      }
    }
  }

  ///// h0
  ///////// 0-1
  double h_0_01_T_i;
  arma::vec h_0_GK_01_T_i;
  arma::vec h_0_GK_01_T0_i;
  if(hazard_baseline_01 == "Exponential"){
    h_0_01_T_i = 1;
    h_0_GK_01_T_i = wk;
    if(left_trunc){
      h_0_GK_01_T0_i = wk;
    }
  }
  if(hazard_baseline_01 == "Weibull"){

    h_0_01_T_i = shape_01*(pow(Time_T_i,(shape_01-1)));
    h_0_GK_01_T_i = shape_01*(pow(st_T_i,shape_01-1))%wk;
    if(left_trunc){
      h_0_GK_01_T0_i = shape_01*(pow(st_T0_i,shape_01-1))%wk;
    }
  }
  if(hazard_baseline_01 == "Gompertz"){
    h_0_01_T_i = Gompertz_1_01*exp(Gompertz_2_01*Time_T_i);
    h_0_GK_01_T_i = Gompertz_1_01*exp(Gompertz_2_01*st_T_i)%wk;
    if(left_trunc){
      h_0_GK_01_T0_i = Gompertz_1_01*exp(st_T0_i*Gompertz_2_01)%wk;
    }
  }
  if(hazard_baseline_01 == "Splines"){
    h_0_01_T_i = exp(arma::dot(gamma_01,B_T_i_01));
    h_0_GK_01_T_i = wk%exp(Bs_T_i_01*gamma_01);
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

  h_01_T_i = h_0_01_T_i*exp(predsurv_01)*h_01_T_i;

  etaBaseline_01_T_i = etaBaseline_01_T_i + predsurv_01;
  survLong_01_T_i = exp(survLong_01_T_i)*h_0_GK_01_T_i;
  arma::vec A_01_T_i;
  A_01_T_i = (exp(etaBaseline_01_T_i)%survLong_01_T_i*(Time_T_i/2));

  arma::vec A_01_T0_i;
  if(left_trunc){
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
    h_0_02_T_i = shape_02*(pow(Time_T_i,(shape_02-1)));
    h_0_GK_02_T_i = shape_02*(pow(st_T_i,shape_02-1))%wk;
    if(left_trunc){
      h_0_GK_02_T0_i = shape_02*(pow(st_T0_i,shape_02-1))%wk;
    }
  }
  if(hazard_baseline_02 == "Gompertz"){
    h_0_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*Time_T_i);
    h_0_GK_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*st_T_i)%wk;
    if(left_trunc){
      h_0_GK_02_T0_i = Gompertz_1_02*exp(st_T0_i*Gompertz_2_02)%wk;
    }
  }
  if(hazard_baseline_02 == "Splines"){
    h_0_02_T_i = exp(arma::dot(gamma_02,B_T_i_02));
    h_0_GK_02_T_i = wk%exp(Bs_T_i_02*gamma_02);
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

  h_02_T_i = h_0_02_T_i*exp(predsurv_02)*h_02_T_i;

  etaBaseline_02_T_i = etaBaseline_02_T_i + predsurv_02;
  survLong_02_T_i = exp(survLong_02_T_i)*h_0_GK_02_T_i;
  arma::vec A_02_T_i;
  A_02_T_i = (exp(etaBaseline_02_T_i)%survLong_02_T_i*(Time_T_i/2));

  arma::vec A_02_T0_i;
  if(left_trunc){
    etaBaseline_02_T0_i = etaBaseline_02_T0_i + predsurv_02;
    survLong_02_T0_i = exp(survLong_02_T0_i)*h_0_GK_02_T0_i;
    A_02_T0_i = (exp(etaBaseline_02_T0_i)%survLong_02_T0_i*(Time_T0_i/2));
  }


  arma::vec SurvTotCase2 =  -A_01_T_i - A_02_T_i + log(pow(h_02_T_i,delta2_i))+ log(pow(h_01_T_i,delta1_i));

  arma::vec f_Y_b_sigma(S,fill::zeros);
  arma::mat X_base_i_id_visit;
  arma::mat U_base_i_id_visit;
  arma::vec y_i_id_visit;
  arma::mat CV_long;
  for(int idvisit = 0; idvisit < len_visit_i; ++idvisit ){
    X_base_i_id_visit = X_base_i.row(idvisit);
    U_base_i_id_visit = U_base_i.row(idvisit);
    y_i_id_visit = y_i.subvec(offset_ID_i(idvisit)-1,offset_ID_i(idvisit+1)-2);
    int n_ij = y_i_id_visit.n_elem;
    CV_long = dot(beta,X_base_i_id_visit) + U_base_i_id_visit*b_y ;
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



  // Rcout << "The value of v : \n" << 9 << "\n";

  arma::vec log_dens_int;
  double Clogexp;
  double log_dens;
  log_dens_int = f_Y_b_sigma + SurvTotCase2;
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
