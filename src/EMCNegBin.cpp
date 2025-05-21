#include <RcppArmadillo.h>
#include "NegBin.h"
#include "CNegBin.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;


double Estep_W_CNegBin(double dY, double dP, double dN, double dEta) {
  
  double dP_y = dCNegBin(dY, dN, dP, dEta, true);
  
  double dW_hat = 0.0;
  
  int w;
  
  double dP_X_given_W = 0.0;
  
  for (w = 1; w <= 100; w++) {
    
    dP_X_given_W = dNegBin(dY, (double)w * dN, dP, true);
    
    dW_hat += exp(log((double)w) + dP_X_given_W + Rf_dpois((double)w, dEta, 1) - dP_y);
  }
  
  return dW_hat;
  
}



//[[Rcpp::export]]
Rcpp::List EM_CNegBin(arma::vec vY, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  int t;
  int j; 
  
  arma::vec LLKSeries(iMaxiter);
  arma::vec vM(iT);
  arma::vec vZ(iT);
  arma::vec vW(iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  //parameter initialization
  double dMean = mean(vY);
  double dVar  = var(vY);

  double dN;
  double dP;
  double dEta = 4;
  
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
  double dEta_Next = dEta;
  
  for (t = 0; t < iT; t++) {
    dLLK += dCNegBin(vY(t), dN, dP, dEta, true);
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
      vW(t) = Estep_W_CNegBin(vY(t), dP, dN, dEta);
    }
    
    //M Step
    dTheta_Next  = accu(vY - vM)/(accu(vY - vM) + accu(vM % vZ));
    dLambda_Next = accu(vM)/accu(vW);
    // dEta_Next    = mean(vW);
    
    dAlpha  = -1.0 / log(1.0 - dTheta_Next);
    
    //update Parameters
    dP = dTheta_Next/(1.0 - dTheta_Next);
    dN = dAlpha * dLambda_Next;
    
    dTheta  = dTheta_Next;
    dLambda = dLambda_Next;
    dEta    = dEta_Next;
    
    dLLK = 0.0;
    for (t = 0; t < iT; t++) {
      dLLK += dCNegBin(vY(t), dN, dP, dEta, true);
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
  lOut["dEta"] = dEta;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  
  return lOut;
  
}