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

double re_lcmm_interintra_cpp(arma::vec beta,
                           arma::mat b_y, List sigma_inter_intra,
                           int len_visit_i, arma::mat X_base_i, arma::mat U_base_i,  arma::vec y_i, arma::vec offset_ID_i
){
  // parameters
  arma::vec sigma_inter = sigma_inter_intra[0];
  arma::vec sigma_intra = sigma_inter_intra[1];
  arma::vec sigma_long = sigma_inter_intra[2];
  arma::vec var_inter = sigma_inter_intra[3];
  arma::vec var_intra = sigma_inter_intra[4];
  arma::vec corr_intra_inter = sigma_inter_intra[5];

  arma::vec f_Y_b_sigma(1,fill::zeros);
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
          arma::vec somme1(1,fill::zeros);
          arma::vec somme2(1,fill::zeros);
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
  log_dens_int = f_Y_b_sigma;
  Clogexp = max(log_dens_int) - 500;
  log_dens_int = log_dens_int - Clogexp;
  log_dens = Clogexp + log(sum(exp(log_dens_int)));

  return log_dens;
}
