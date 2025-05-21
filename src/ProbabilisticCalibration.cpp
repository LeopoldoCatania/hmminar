#include <RcppArmadillo.h>
#include "PoissonBinomial.h"

using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
double PIT_MixINAR(double dU, double dY, 
                   double dX, double dAlpha, arma::vec vLambda, arma::vec vOmega, arma::vec vLogFactorial_K) {

  double dPIT = 0.0;

  if (dY == 0) {

    return dPIT;

  }

  double dP_y = pMixPB(dY, dX, dAlpha, vLambda, vOmega, vLogFactorial_K, false);

  double dP_y_m1 = pMixPB(dY - 1.0, dX, dAlpha, vLambda, vOmega, vLogFactorial_K, false);

  if (dU >= dP_y) {

    dPIT = 1.0;
    return dPIT;

  }

  if (dU <= dP_y_m1) {

    dPIT = 0.0;
    return dPIT;

  }

  dPIT = (dU - dP_y_m1) / (dP_y - dP_y_m1);

  return dPIT;

}

//[[Rcpp::export]]
double AveragePIT_MixINAR(double dU, arma::vec vY, 
                          double dAlpha, arma::vec vLambda, arma::vec vOmega, arma::vec vBeta, 
                          arma::mat mD, arma::vec vLogFactorial_K) {
  
  int iT = vY.size();
  int iD = mD.n_cols;
  int t;
  int d;
  
  mD = mD * 1.0;
  
  int iSeason = 0;
  
  double dOut = 0.0;
  
  for(t = 1; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        iSeason = d;
      }
    }
    dOut += PIT_MixINAR(dU, vY(t), vY(t - 1), dAlpha, vLambda * vBeta(iSeason), vOmega, vLogFactorial_K);
  }
  
  dOut = dOut/(iT * 1.0 - 1.0);
  
  return dOut;
}

//[[Rcpp::export]]
List PIT_Bin_MixINAR(int iB, arma::vec vY, 
                  double dAlpha, arma::vec vLambda, arma::vec vOmega, arma::vec vBeta, 
                  arma::mat mD, arma::vec vLogFactorial_K) {
  
  arma::vec vPIT(iB);
  
  double dF_b_iB = 0.0;
  double dF_bm1_iB = 0.0;
  
  int b;
  double db;
  double dB = iB * 1.0;
  
  for (b = 1; b <= iB; b++) {
    
    db = (double)b;
    
    if (b == 1) {
      dF_bm1_iB = AveragePIT_MixINAR((db - 1.0)/dB, vY, dAlpha, vLambda, vOmega, vBeta, mD, vLogFactorial_K);
    } else {
      dF_bm1_iB = dF_b_iB;
    }
    
    dF_b_iB   = AveragePIT_MixINAR(db/dB, vY, dAlpha, vLambda, vOmega, vBeta, mD, vLogFactorial_K);
   
    vPIT(b - 1) = dF_b_iB - dF_bm1_iB;
  }
    
  
  int iDegreeFreedom = iB - (vBeta.size() + vLambda.size() + vOmega.size() - 1.0);
  
  double dTest = accu(pow(vPIT - 1.0/dB, 2.0)/(1.0/dB));
  double dPVal = -999;
  
  if (iDegreeFreedom > 0) {
    dPVal = 1.0 - Rf_pchisq(dTest, iDegreeFreedom * 1.0, 1, 0);
  }
  
  List lOut;
  
  lOut["vPIT"] = vPIT;
  lOut["dTest"] = dTest;
  lOut["dPVal"] = dPVal;
  
  return lOut;
}

