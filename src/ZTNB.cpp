#include <RcppArmadillo.h>
#include "Utils.h"
#include "SafeFunctions.h"
// #include <omp.h>

using namespace arma;
using namespace Rcpp;


//[[Rcpp::export]]
double dZTNB(double dV, double dA, double dB, bool bLog = true) {

  double dPDF = Rf_lgammafn(dA + dV) - Rf_lgammafn(dA) - lfactorial2(dV) + dV * log(dB) - log(pow(1.0 - dB, -dA) - 1.0);

  if (!bLog) {
    dPDF = exp(dPDF);
  }

  return dPDF;

}

//[[Rcpp::export]]
double dSumZTNB(double dV, double dN, double dA, double dB, bool bLog = true) {

  double dLPDF = 0.0;

  double dPDF1 = -INFINITY;
  double dPDF2 = -INFINITY;

  arma::vec vFoo(2);

  int r;
  int iN = (int)dN;

  double dr;

  double dFoo = 0.0;

  double dConst = -lfactorial2(dV) + dV * log(dB) - dN * log(pow(1.0 - dB, -dA) - 1.0);

  bool   bSum1 = false;

  double dL;

  if ((iN == 0) & (dV > 0)) {
    dLPDF = -INFINITY;
  } else if ((iN == 0) & (dV == 0)) {
    dLPDF = -INFINITY;
  } else if ((iN == 0) & (dV == 0)) {
    dLPDF = 0;
  } else if (dV < dN) {
    dLPDF = -INFINITY;
  } else if (iN == 1) {

    dLPDF = Rf_lgammafn(dA + dV) - Rf_lgammafn(dA) - lfactorial2(dV) + dV * log(dB) - log(pow(1.0 - dB, -dA) - 1.0);

  } else {

    if((iN - 1) % 2 == 0) {
      bSum1 = true;
    }

    for (r = 1; r <= iN; r++) {

      dr = r * 1.0;

      if (r == 1) {
        dFoo = log(dN);
      } else if (r == iN) {
        dFoo = 0.0;
      } else {
        dFoo += (log(dN + 1.0 - dr) - log(dr));
      }

      vFoo(1) = dFoo + Rf_lgammafn(dV + dr * dA) - Rf_lgammafn(dr * dA) + dConst;

      if(bSum1) {
        vFoo(0) = dPDF1;
        dPDF1 = LogSumExp(vFoo);
      } else {
        vFoo(0) = dPDF2;
        dPDF2 = LogSumExp(vFoo);
      }

      bSum1 = !bSum1;

    }

    dL = dPDF1 - dPDF2;

    dLPDF = dPDF1 + log(1.0 - exp(-dL));

    if (!bLog) {
      dLPDF = exp(dLPDF);
    }
  }

  return dLPDF;

}

//[[Rcpp::export]]
double dLLK_SumZTNB(arma::vec vV, arma::vec vN, double dA, double dB){

  double dLLK = 0.0;

  int iT = vV.size();

  int t;
  double dFoo;

  for (t = 0; t < iT; t++) {
    dFoo = dSumZTNB(vV(t), vN(t), dA, dB, true);
    if (R_IsNaN(dFoo)) {
      dFoo = -9999.0;
    }
      dLLK += dFoo;
  }

  return dLLK;

}
