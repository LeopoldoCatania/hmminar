#include <RcppArmadillo.h>
#include "Reparameterization.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List Moments_C(arma::mat PredictedProb, 
             int iJ, int iL,  
             arma::vec vMu, arma::vec vNu, 
             arma::vec vOmega,
             arma::mat mGamma,
             arma::vec vBeta,
             arma::mat mD) {
  
  int iQ = iL * iJ;
  int iD = mD.n_cols;
  
  int iT = PredictedProb.n_rows - 1;
  
  arma::mat mGamma_enlarged = Gamma_enlarged(mGamma, vOmega);
  arma::vec vDelta          = getDelta(mGamma_enlarged, iQ);
  
  arma::mat mIndices = Indices(iJ, iL);
  
  arma::vec vMean(iT);
  arma::vec vVar(iT);
  arma::vec vE2(iT);
  
  int t;
  int q;
  int d;
  
  vMean.zeros();
  vVar.zeros();
  vE2.zeros();
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  vSeason.zeros();
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }
  
  for (t = 0; t < iT; t++) {
    for (q = 0; q < iQ; q++) {
      vMean(t) += (PredictedProb(t, q) * vNu((int)mIndices(q, 1)) * vMu((int)mIndices(q, 0)) * vBeta(vSeason(t)));
      vE2(t)   += (PredictedProb(t, q) * vMu((int)mIndices(q, 0)) * vBeta(vSeason(t)) * vNu((int)mIndices(q, 1)) * 
        (vNu((int)mIndices(q, 1)) * (1.0 + vMu((int)mIndices(q, 0)) * vBeta(vSeason(t))) + 1.0));
    }
    vVar(t) = vE2(t) - pow(vMean(t), 2.0);
  }
  
  List lM;
  
  lM["vVar"]   = vVar;
  lM["vMean"]  = vMean;
  
  return lM;
  
}