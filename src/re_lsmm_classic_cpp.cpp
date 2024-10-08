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

double re_lsmm_classic_cpp(arma::vec beta, arma::mat b_y,
                          arma::mat X_base_i, arma::mat U_base_i,   arma::vec y_i, double sigma_epsilon
){
  arma::vec f_Y_b_sigma(1,fill::zeros);
  arma::vec CV;
  int n_rows_X = X_base_i.n_rows;
  double sigma_long;
  sigma_long = sigma_epsilon;
  for(int k=0; k<n_rows_X; k++){
    CV = dot(beta,X_base_i.row(k)) + b_y*U_base_i.row(k).t();
    f_Y_b_sigma = f_Y_b_sigma + log(1.0 / (sqrt(2.0*M_PI)*sigma_long)) - 0.5*pow((y_i(k)-CV)/sigma_long,2);
  }

  arma::vec log_dens_int;
  double Clogexp;
  double log_dens;
  log_dens_int = f_Y_b_sigma;
  Clogexp = max(log_dens_int) - 500;
  log_dens_int = log_dens_int - Clogexp;
  log_dens = Clogexp + log(sum(exp(log_dens_int)));

  return log_dens;
}
