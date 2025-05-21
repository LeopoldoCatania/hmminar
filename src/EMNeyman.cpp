#include <RcppArmadillo.h>
#include "Neyman.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List EM_Neyman(arma::vec vY, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  int t;
  
  arma::vec LLKSeries(iMaxiter);
  arma::vec vK(iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  //parameter initialization
  double dMean = mean(vY);
  double dVar  = var(vY);

  double dMu;
  double dNu;
  
  if (dVar > dMean) {
    dMu = pow(dMean, 2.0)/(dVar - dMean);
    dNu = dVar/dMean - 1.0;
  } else {
    dMu = dMean;
    dNu = dVar;
  }
  
  double dMu_Next = dMu;
  double dNu_Next = dNu;
  
  for (t = 0; t < iT; t++) {
    dLLK += dNeyman(vY(t), dNu, dMu, true);
  }
  
  LLKSeries(0) = dLLK;
  
  int iter = 1;
  
  while (dEps > dTol && iter < iMaxiter) {
    //E step
    for (t = 0; t < iT; t++) {
      vK(t) = (vY(t) + 1.0)/(dNu) * exp(dNeyman(vY(t) + 1.0,  dNu, dMu, true) - dNeyman(vY(t), dNu, dMu, true));
    }
    
    //M Step
    dMu_Next = mean(vK);
    dNu_Next = dMean/dMu_Next;
    
    //update Parameters
    dMu = dMu_Next;
    dNu = dNu_Next;
    
    dLLK = 0.0;
    for (t = 0; t < iT; t++) {
      dLLK += dNeyman(vY(t), dNu, dMu, true);
    }
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["dNu"] = dNu;
  lOut["dMu"] = dMu;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  lOut["dEps"] = dEps;
  
  return lOut;
  
}