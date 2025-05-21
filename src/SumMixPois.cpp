#include <RcppArmadillo.h>
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

double dMixPois(double dY, arma::vec vOmega, arma::vec vNu, bool bLog = true) {
  
  int iJ = vOmega.size();
  int j;
  
  arma::vec vD_log(iJ);
  
  for (j = 0; j < iJ; j++) {
    vD_log(j) = Rf_dpois(dY, vNu(j), 1);
  }
  
  double dLPMF = MixtDensityScale(log(vOmega), vD_log);
  
  if (!bLog) {
    dLPMF = exp(dLPMF);
  }
  
  return dLPMF;
  
}

double rMixPois(arma::vec vOmega, arma::vec vNu) {
  
  int index = rando_index(vOmega);
  double dY = Rf_rpois(vNu(index));
  
  return dY;
  
}

// //[[Rcpp::export]]
// double dSumMixPois(double dY, arma::vec vOmega, arma::vec vNu, double dMu, bool bLog = true, int iTrunc = 40) {
//   
//   arma::vec vPOI_log(iTrunc);
//   arma::vec vMixPOI_log(iTrunc);
//   
//   int i;
//   
//   for (i = 0; i < iTrunc; i++) {
//     
//     vPOI_log(i) = Rf_dpois((double)i, dMu, 1);
//     
//     if (i == 0) {
//       vMixPOI_log(i)  = dMixPois(dY, vOmega, vNu * 1e-8, true);
//     } else {
//       
//       vMixPOI_log(i)  = dMixPois(dY, vOmega, vNu * (double)i, true);
//     }
//   }
//   
//   double dLPMF = MixtDensityScale(vPOI_log, vMixPOI_log);
//   
//   if(!bLog) {
//     dLPMF = exp(dLPMF);
//   }
//   
//   return dLPMF;
//   
// }

//[[Rcpp::export]]
double dSumMixPois(double dY, arma::vec vOmega, arma::vec vNu, double dMu, bool bLog = true, int iTrunc = 100) {
  
  int iJ = vNu.size();
  
  arma::vec vFoo(iJ);
  arma::vec vBaz(iTrunc - 1);
  double dLPMF = 0.0;
  
  vFoo.zeros();
  
  double dLogFoo = 0.0;
  
  int i;
  int j;
  
  double dK;
  double dlK;
  double dlMu= log(dMu);
  double dLogFactY = lfactorial2(dY);
  
  double dStop = dMu - 25.0;
  
  double dPoi_Prop_P = 10;
  
  for (j = 0; j < iJ; j++) {
    dLogFoo = 0.0;
    dPoi_Prop_P = 10.0;
    i = 1;
    
    while (dPoi_Prop_P > dStop && i < iTrunc) {
      
      dK = i * 1.0;
      
      dlK = log(dK);
      
      dLogFoo += dlK;
      
      dPoi_Prop_P = dK * dlMu  - dLogFoo;
      
      vBaz(i - 1) = dPoi_Prop_P - dK * vNu(j) + dY * dlK;
      
      i += 1;
    }
    
    vFoo(j) = dY * log(vNu(j)) + log(vOmega(j)) + LogSumExp(vBaz.subvec(0, i - 2));
    
  }
  
  dLPMF = - dMu - dLogFactY + LogSumExp(vFoo);
  
  if(!bLog) {
    dLPMF = exp(dLPMF);
  }
  
  return dLPMF;
  
}






//[[Rcpp::export]]
double dLLK_SumMixPois(arma::vec vY, arma::vec vOmega, arma::vec vNu, double dMu) {
  double dLLK = 0.0;
  
  int iT = vY.size();
  
  for (int i = 0; i < iT; i++) {
    dLLK += dSumMixPois(vY(i), vOmega, vNu, dMu);
  }
  
  return dLLK;
  
}

//[[Rcpp::export]]
double rSumMixPois(arma::vec vOmega, arma::vec vNu, double dMu) {
  
  double dK = Rf_rpois(dMu);
  
  double dY = 0.0;
  
  if (dK > 0) {
    dY = rMixPois(vOmega, vNu * dK);
  }
  
  return dY;
  
}
