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

arma::vec log_llh_lsmm_classic_cpp(int S,  arma::vec beta, arma::mat b_al,
                                   arma::mat X_base,  arma::mat U_base, arma::vec y_new, int Ind,
                                   arma::vec offset, double sigma_epsilon){



  arma::vec ll_glob(Ind,fill::ones);

  for(int i_index =0; i_index < Ind; ++i_index){

    arma::mat X_base_i = X_base.rows((offset(i_index)-1),(offset(i_index+1)-2));
    arma::mat U_base_i = U_base.rows((offset(i_index)-1),(offset(i_index+1)-2));
    arma::vec y_i = y_new.subvec(offset(i_index)-1,offset(i_index+1)-2);

    arma::vec f_Y_b_sigma(S,fill::zeros);
    arma::vec CV;
    int n_rows_X = X_base_i.n_rows;
    double sigma_long;
    sigma_long = sigma_epsilon;
    for(int k=0; k<n_rows_X; k++){
      CV = dot(beta,X_base_i.row(k)) + b_al*U_base_i.row(k).t();
      f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i(k)-CV)/sigma_long,2);
    }
    //  if(n_row_X == 0){
    //    CV  = dot(beta,X_base_i) + b_al*U_i ;
    //    f_Y_b_sigma = log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i-CV)/sigma_long,2);
    //  }
    //  else{
    //    for(int k=0; k<n_row_X; k++){
    //      CV = dot(beta,X_base_i.row(k)) + b_al*U_i.row(k).t();
    //      f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i(k)-CV)/sigma_long,2);
    //    }
    arma::vec log_dens_int;
    double Clogexp;
    double log_dens;
    log_dens_int = f_Y_b_sigma;
    Clogexp = max(log_dens_int) - 500;
    log_dens_int = log_dens_int - Clogexp;
    log_dens = Clogexp + log(sum(exp(log_dens_int))) - log(S);

    ll_glob(i_index) = log_dens;
  }
  return ll_glob;
}
