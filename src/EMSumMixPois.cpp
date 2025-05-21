#include <RcppArmadillo.h>
#include "SumMixPois.h"
#include "Utils.h"
#include "Neyman.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List EM_SumMixPois(arma::vec vY, int iJ, 
                         double dMu, arma::vec vNu, 
                         arma::vec vOmega,
                         double dTol = 1e-6, int iMaxiter = 4000) {
  
  // definitions
  int iT = vY.size();
  int t;
  int j;
  
  arma::vec LLKSeries(iMaxiter);
  arma::mat mLLK(iJ, iT);
  arma::vec vLLK(iT);
  arma::mat mK(iJ, iT);
  arma::mat mZ(iJ, iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  double dMu_Next       = dMu;
  arma::vec vNu_Next    = vNu;
  arma::vec vOmega_Next = vOmega;
  
  for (t = 0; t < iT; t++) {
    dLLK += dSumMixPois(vY(t), vOmega, vNu, dMu, true); 
  }
  
  LLKSeries(0) = dLLK;
  
  int iter = 1;
  
  double dLPMF_foo = 0.0;
  
  while (dEps > dTol && iter < iMaxiter) {
    //E step
    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        mLLK(j, t) = dNeyman(vY(t), vNu(j), dMu, true);
        mK(j, t) = (vY(t) + 1.0)/vNu(j) * exp(dNeyman(vY(t) + 1.0, vNu(j), dMu, true) - mLLK(j, t));
      }
      
      dLPMF_foo = MixtDensityScale(log(vOmega), mLLK.col(t));
      vLLK(t) = dLPMF_foo; 
      
      for (j = 0; j < iJ; j++) {
        mZ(j, t) = exp(mLLK(j, t) + log(vOmega(j)) - dLPMF_foo);
      }
    }
    
    //M Step
    dMu_Next = 0.0;
    for (j = 0; j < iJ; j++) {
      dMu_Next       += accu(mK.row(j) % mZ.row(j));
      vNu_Next(j)    = accu(vY % mZ.row(j).t()) / accu(mK.row(j) % mZ.row(j));
      vOmega_Next(j) = accu(mZ.row(j));
    }
    
    dMu_Next = dMu_Next / (iT * 1.0);

    vOmega_Next = vOmega_Next/accu(vOmega_Next);
    
    //update Parameters
    dMu = dMu_Next;
    vNu = vNu_Next;
    vOmega = vOmega_Next;
    
    dLLK = 0.0;
    for (t = 0; t < iT; t++) {
      dLLK += dSumMixPois(vY(t), vOmega, vNu, dMu, true);
    }
   // dLLK = accu(vLLK);
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["vNu"] = vNu;
  lOut["vOmega"] = vOmega;
  lOut["dMu"] = dMu;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  lOut["dEps"] = dEps;
  
  lOut["mZ"] = mZ;
  
  return lOut;
  
}