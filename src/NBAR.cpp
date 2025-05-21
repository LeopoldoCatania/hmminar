#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double dNB(double dY, double dDelta, double dBeta, bool bLog = true) {

  NumericVector vFoo(1);
  vFoo[0] = dY;
  NumericVector vBaz = lfactorial(vFoo);
  double dConst = vBaz[0];

  double dLPMF = dY * log(dBeta) + Rf_lgammafn(dDelta + dY) - Rf_lgammafn(dDelta) - (dDelta + dY) * log(dBeta + 1.0) - dConst;

  if(!bLog) {
    dLPMF = exp(dLPMF);
  }
  return dLPMF;

}

//[[Rcpp::export]]
double pNB(double dY, double dDelta, double dBeta, bool bLog = true) {

  double dCMF = 0.0;
  int iY = (int)dY;

  if(dY > -1e-5) {
    for (int y = 0; y <= iY; y++) {
      dCMF += dNB((double)y, dDelta, dBeta, false);
    }
  }
  return dCMF;
}


// Negative Binomial Autoregression
//[[Rcpp::export]]
List NBAR(arma::vec vY, double dBeta, double dDelta0, double dOmega0) {

  int iT = vY.size();
  arma::vec vLLK(iT);
  arma::vec vDelta(iT + 1);
  int t;

  vDelta(0) = dDelta0;
  vLLK(0) = dNB(vY(0), vDelta(0), dBeta, true);
  double dLLK = vLLK(0);

  for (t = 1; t < iT+1; t++) {

    vDelta(t) = dOmega0 + vY(t-1);

    if(t < iT) {
      vLLK(t) = dNB(vY(t), vDelta(t), dBeta, true);
      dLLK += vLLK(t);
    }
  }

  List lOut;
  lOut["vLLK"] = vLLK;
  lOut["dLLK"] = dLLK;
  lOut["vDelta"] = vDelta;

  return lOut;
}


// Periodic Negative Binomial Autoregression
//[[Rcpp::export]]
List PNBAR(arma::vec vY, double dBeta, double dDelta0, double dOmega0,
           double dOmega1, double dOmega2, int iS) {

  int iT = vY.size();
  arma::vec vLLK(iT);
  arma::vec vLogOmega_S(iS);
  arma::vec vDelta(iT + 1);
  int t;

  double dS = 1.0 * iS;

  for (t = 0; t < iS; t++) {
    vLogOmega_S(t) = dOmega0 + dOmega1 * cos(2.0 * M_PI * (t*1.0 + 1.0)/dS - dOmega2 * M_PI);
  }

  vDelta(0) = dDelta0;
  vLLK(0) = dNB(vY(0), vDelta(0), dBeta, true);
  double dLLK = vLLK(0);

  int iC = 0;
  for (t = 1; t < iT+1; t++) {

    vDelta(t) = vY(t-1) + exp(vLogOmega_S(iC));

    if(iC < iS-1) {
      iC += 1;
    }else {
      iC = 0;
    }
    if(t < iT) {
      vLLK(t) = dNB(vY(t), vDelta(t), dBeta, true);
      dLLK += vLLK(t);
    }
  }

  List lOut;
  lOut["vLLK"] = vLLK;
  lOut["dLLK"] = dLLK;
  lOut["vDelta"] = vDelta;
  lOut["vLogOmega_S"] = vLogOmega_S;

  return lOut;
}

