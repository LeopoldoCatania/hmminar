#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double mapC(double x , double L,double U) {
  double map =  L + ( U - L ) / (1.0 + exp( - x ));
  return map;

}

//[[Rcpp::export]]
double UnMapC(double x , double L,double U) {
  double unmap =  log((x-L)/(U-x));;
  return unmap;

}

//[[Rcpp::export]]
arma::vec mapC_v(arma::vec x , double L,double U) {

  int K=x.size();
  arma::vec map(K);

  for(int i = 0;i<K;i++){
    map(i) = mapC(x(i) , L, U);
  }
  return map;

}

//[[Rcpp::export]]
arma::vec UnmapC_v(arma::vec x , double L,double U) {

  int K=x.size();
  arma::vec unmap(K);

  for(int i = 0;i<K;i++){
    unmap(i) = UnMapC(x(i) , L, U);
  }
  return unmap;

}

//[[Rcpp::export]]
arma::vec SimplexMapping(arma::vec vOmega_tilde, int M){

  int j;

  arma::vec vOmega(M);

  double Scale = 1.0;

  for (j = 0; j < M - 1; j++) {

    vOmega(j) = mapC(vOmega_tilde(j), 0.0, Scale);
    Scale    -=  vOmega(j);
  }

  vOmega(M - 1) = Scale;

  return vOmega;
}


//[[Rcpp::export]]
arma::vec SimplexUnMapping(arma::vec vOmega, int M){

  int j;

  arma::vec vOmega_tilde(M - 1);

  double Scale = 1.0;

  for (j = 0; j < M - 1; j++) {

    vOmega_tilde(j) = UnMapC(vOmega(j), 0.0, Scale);
    Scale    -=  vOmega(j);

  }

  return vOmega_tilde;
}
