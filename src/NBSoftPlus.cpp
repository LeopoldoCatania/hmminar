#include <RcppArmadillo.h>
#include "NBAR.h"

using namespace arma;
using namespace Rcpp;

double SoftPlus(double dX, double dC) {
  double dOut = dC * log(1 + exp(dX/dC));
  return dOut;
}

//[[Rcpp::export]]
double dNBin(double dY, double dN, double dP, bool bLog = true) {

  NumericVector vFoo(1);
  vFoo[0] = dY;
  NumericVector vBaz = lfactorial(vFoo);
  double dConst = vBaz[0];

  double dLPMF = Rf_lgammafn(dN + dY) - Rf_lgammafn(dN) - dConst + dN * log(dP) + dY * log(1.0 - dP);

  if(!bLog) {
    dLPMF = exp(dLPMF);
  }
  return dLPMF;

}

//[[Rcpp::export]]
double pNBin(double dY, double dN, double dP, bool bLog = true) {

  double dCMF = 0.0;
  int iY = (int)dY;

  if(dY > -1e-5) {
    for (int y = 0; y <= iY; y++) {
      dCMF += dNBin((double)y, dN, dP, false);
    }
  }
  return dCMF;
}



//[[Rcpp::export]]
List NBSoftplusSeasonalINGARCH_filter(arma::vec vY, double dM0,
                                      double dAlpha0,
                                      double dAlpha1, double dBeta1,
                                      arma::vec vDelta_S,
                                      int iS, double dC,
                                      double dP) {

  int iT = vY.size();
  arma::vec vM(iT + 1);
  arma::vec vDeltaPlusM(iT + 1);
  arma::vec vN(iT + 1);
  arma::vec vLLK(iT);
  vM(0) = dM0;
  double dLLK = 0.0;

  int t;

  //conditional mean
  vDeltaPlusM(0) = vM(0) + vDelta_S(0);
  //number of success
  vN(0) = vDeltaPlusM(0) * dP/(1.0 - dP);

  vLLK(0) = dNBin(vY(0), vN(0), dP, true);
  dLLK += vLLK(0);
  double dFoo = 0.0;

  int iC = 1;
  for(t = 1; t < iT + 1; t++) {

    dFoo = dAlpha0 + dAlpha1 * (1.0 * vY(t - 1)) + dBeta1 * vM(t - 1);
    vM(t) = SoftPlus(dFoo, dC);

    vDeltaPlusM(t) = vM(t) + vDelta_S(iC);
    vN(t) = vDeltaPlusM(t) * dP/(1.0 - dP);

    if(iC < iS-1) {
      iC += 1;
    }else {
      iC = 0;
    }
    if(t < iT) {
      vLLK(t) = dNBin(vY(t), vN(t), dP, true);
      dLLK += vLLK(t);
    }

  }

  List lOut;
  lOut["vLLK"] = vLLK;
  lOut["dLLK"] = dLLK;
  lOut["vDeltaPlusM"] = vDeltaPlusM;
  lOut["vDelta_S"] = vDelta_S;
  lOut["vM"] = vM;
  lOut["vN"] = vN;

  return lOut;

}


//[[Rcpp::export]]
List NBSoftplusSeasonalINGARCH_DummySeason_filter(arma::vec vY, double dM0,
                                      double dAlpha0,
                                      double dAlpha1, double dBeta1,
                                      arma::vec vSeason, arma::vec vBeta_s,
                                      double dC,
                                      double dP) {

  int iT = vY.size();
  arma::vec vM(iT);
  arma::vec vDeltaPlusM(iT);
  arma::vec vN(iT);
  arma::vec vLLK(iT);
  vM(0) = dM0;
  double dLLK = 0.0;

  int t;

  //conditional mean
  vDeltaPlusM(0) = vM(0) + vBeta_s(vSeason(0));
  //number of success
  vN(0) = vDeltaPlusM(0) * dP/(1.0 - dP);

  vLLK(0) = dNBin(vY(0), vN(0), dP, true);
  dLLK += vLLK(0);
  double dFoo = 0.0;

  for(t = 1; t < iT; t++) {

    dFoo = dAlpha0 + dAlpha1 * (1.0 * vY(t - 1)) + dBeta1 * vM(t - 1);
    vM(t) = SoftPlus(dFoo, dC);

    vDeltaPlusM(t) = vM(t) + vBeta_s(vSeason(t));
    vN(t) = vDeltaPlusM(t) * dP/(1.0 - dP);

    if(t < iT) {
      vLLK(t) = dNBin(vY(t), vN(t), dP, true);
      dLLK += vLLK(t);
    }
  }

  List lOut;
  lOut["vLLK"] = vLLK;
  lOut["dLLK"] = dLLK;
  lOut["vDeltaPlusM"] = vDeltaPlusM;
  lOut["vBeta_s"] = vBeta_s;
  lOut["vM"] = vM;
  lOut["vN"] = vN;
  lOut["vSeason"] = vSeason;

  return lOut;

}
