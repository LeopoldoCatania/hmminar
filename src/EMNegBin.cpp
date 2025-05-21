#include <RcppArmadillo.h>
#include "NegBin.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List EM_NegBin(arma::vec vY, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  int t;
  int j;
  
  arma::vec LLKSeries(iMaxiter);
  arma::vec vM(iT);
  arma::vec vZ(iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  //parameter initialization
  double dMean = mean(vY);
  double dVar  = var(vY);

  double dN;
  double dP;
  
  if (dVar > dMean) {
    dP = 1.0 - dMean/dVar;
    dN = pow(dMean, 2.0)/(dVar - dMean);
  } else {
    dP = 0.4;
    dN = 1.5;
  }
  
  double dTheta  = dP/(1.0 + dP);
  double dLambda = -dN * log(1.0 - dTheta); 
  double dAlpha  = -1.0 / log(1.0 - dTheta);
  
  double dTheta_Next  = dTheta;
  double dLambda_Next = dLambda;
  
  for (t = 0; t < iT; t++) {
    dLLK += dNegBin(vY(t), dN, dP, true);
  }
  
  LLKSeries(0) = dLLK;
  
  int iter = 1;
  
  while (dEps > dTol && iter < iMaxiter) {
    //E step
    for (t = 0; t < iT; t++) {
      vM(t) = 0.0;
      if (vY(t) > 0) {
        for (j = 1; j <= (int)vY(t); j++) {
          vM(t) += (1.0/(dAlpha * dLambda + j * 1.0 - 1.0));
        }
        vM(t) = vM(t) * dAlpha * dLambda;
      }
      // vM(t) = dAlpha * dLambda * (Rf_digamma(dAlpha * dLambda + vY(t)) - Rf_digamma(dAlpha * dLambda));
      vZ(t) = -((1.0 - dTheta)/dTheta - dAlpha);
    }
    
    //M Step
    dTheta_Next  = accu(vY - vM)/(accu(vY - vM) + accu(vM % vZ));
    dLambda_Next = mean(vM);
    
    dAlpha  = -1.0 / log(1.0 - dTheta_Next);
    
    //update Parameters
    dP = dTheta_Next/(1.0 - dTheta_Next);
    dN = dAlpha * dLambda_Next;
    
    dTheta  = dTheta_Next;
    dLambda = dLambda_Next;
    
    dLLK = 0.0;
    for (t = 0; t < iT; t++) {
      dLLK += dNegBin(vY(t), dN, dP, true);
    }
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["dP"] = dP;
  lOut["dN"] = dN;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  
  return lOut;
  
}