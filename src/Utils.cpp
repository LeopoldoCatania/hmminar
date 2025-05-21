#include <RcppArmadillo.h>
#include "SafeFunctions.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double MixtDensityScale(arma::vec vOmega_log, arma::vec vD_log){

  int iJ = vOmega_log.size();

  arma::vec wp_log = vOmega_log + vD_log;

  double dK = max(wp_log );

  arma::vec wp_log_scaled = wp_log - dK;

  double dLK = 0;
  int j;
  for (j = 0; j < iJ; j++) {
    dLK += exp(wp_log_scaled(j));
  }

  double dLLK = dK + log(dLK);

  dLLK = LogDensityCheck(dLLK);

  return dLLK;
}

//[[Rcpp::export]]
double LogSumExp(arma::vec vLog) {

  double dMax = max(vLog);

  double dLogSumExp = dMax + log(accu(exp(vLog - dMax)));

  return dLogSumExp;

}

int WhichMax(arma::vec vX){

  int iK = vX.size();
  int k;

  int iMax = 0;
  double dMax = vX(0);

  for (k = 1; k < iK; k++) {
    if (vX(k) > dMax) {
      iMax = k;
      dMax = vX(k);
    }
  }

  return iMax;
}

double factorial2(double dX) {

  double dF = Rcpp::internal::factorial(dX);

  return dF;

}

double lfactorial2(double dX) {

  double dlF = Rcpp::internal::lfactorial(dX);

  return dlF;

}
double abs3(double x){
  double abs_x = x;
  if (abs_x < 0) {
    abs_x = -abs_x;
  }
  return abs_x;
}

int rando_index(arma::vec p) {
  int N = p.size();

  double u = R::runif(0, 1.0);
  int    i = 0;
  double s = p(0);

  while((u > s) && (i < N - 1)) {
    i++;
    s += p(i);
  }

  return i;
}

int rando_index2(arma::vec p) {
  int N = p.size();

  double u = R::runif(0, 1.0);
  int    i = 0;
  double s = p(0);

  while((u > s) && (i < N - 1)) {
    i++;
    s += p(i);
  }

  return i;
}

//[[Rcpp::export]]
arma::vec getDelta(arma::mat gamma, int m){
  arma::mat I = eye(m, m);
  arma::mat Umat = ones(m, m);

  arma::vec Uvec(m);Uvec.fill(1);

  arma::mat foo = (I - gamma + Umat).t();

  arma::mat delta = (foo).i() * Uvec;

  // boundaries
  delta = DeltaCheck(delta, m);

  return delta;
}

arma::mat AccessListMat(List list, std::string element_name){
  SEXP foo = wrap(as<NumericMatrix>(list[element_name]));
  arma::mat mat_out = Rcpp::as<arma::mat>(foo);
  return mat_out;
}

//[[Rcpp::export]]
double LW(double x, double tolerance = 1e-10, int maxit = 50)
{
  double ans = x;
  ans = pow(x, 0.5)/2;
  double cutpt = 3.0;

  double L1 = 0.0;
  double L2 = 0.0;
  double wzinit = 0.0;
  double exp1 = 0.0;
  double exp2 = 0.0;
  double delta = 0.0;

  if (x > cutpt) {
    L1 = log(x);
    L2 = log(L1);
    wzinit = L1 - L2 + (L2 + (L2 * (-2.0 + L2)/(2.0) + (L2 *
      (6.0 + L2 * (-9.0 + L2 * 2.0))/(6.0) + L2 * (-12.0 + L2 * (36.0 +
      L2 * (-22.0 + L2 * 3.0)))/(12.0 * L1))/L1)/L1)/L1;
    ans = wzinit;
  }

  int ii = 1;

  while(ii <= maxit && delta > tolerance) {

    exp1 = exp(ans);
    exp2 = ans * exp1;
    delta = (exp2 - x)/(exp2 + exp1 - ((ans + 2.0) * (exp2 -  x)/(2.0 * (ans + 1.0))));
    ans = ans - delta;

    ii += 1;

  }

  if(ii == maxit) {
    printf("Lambert function did not converge");
  }

  return ans;
}

//[[Rcpp::export]]
arma::mat OmegaZB_Indices(int iJ, int iK, int iL) {

  arma::mat mIndices = zeros(iJ * iK * iL, 3);

  int j,k,l;
  int iC = 0;

  for(l = 0; l < iL; l++) {
    for(j = 0; j < iJ; j++) {
      for(k = 0; k < iK; k++) {

        mIndices(iC, 0) = j;
        mIndices(iC, 1) = k;
        mIndices(iC, 2) = l;

        iC += 1;

      }
    }
  }

  return mIndices;

}

double LogOfSum(arma::vec vLog) {

  double c    = max(vLog);
  double dLogSum = c + log(sum(exp(vLog - c)));

  return dLogSum;

}

//[[Rcpp::export]]
arma::vec BoundProbs(arma::vec vP, double probbound) {

  int iJ = vP.size();
  for (int j = 0; j < iJ; j++) {
    if(vP(j) < probbound) {
      vP(j) = probbound;
    }
  }

  vP = vP/accu(vP);

  return vP;
}

//[[Rcpp::export]]
arma::mat BoundGamma(arma::mat mGamma, double probbound){

  int iJ = mGamma.n_cols;
  for (int j = 0; j < iJ; j++) {
    mGamma.row(j) = arma::trans(BoundProbs(mGamma.row(j).t(), probbound));
  }

  return mGamma;
}
