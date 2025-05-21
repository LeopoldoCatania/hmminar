#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

// //[[Rcpp::export]]
// double dNeyman_Int(double dY, double dNu, double dMu) {
//   
//   double dLambda = dMu * exp(-dNu);
//   
//   double PDF = 0.0;
//   double dFoo = 0.0;
//   int r;
//   
//   if (dY < 1) {
//     PDF = exp(-dMu + dLambda);
//   } else {
//     
//     for(r = 0; r <= ((int)dY - 1); r++){
//       dFoo += (pow(dNu, (double)r)/Rcpp::internal::factorial((double)r) * dNeyman_Int(dY - (double)r - 1, dNu, dMu));
//     }
//     
//     PDF = dNu * dLambda * dFoo/ (dY);
//   }
//   
//   return PDF;
//   
// }
// //[[Rcpp::export]]
// double dNeyman(double dY, double dNu, double dMu, bool bLog = true, int iTrunc = 100) {
//   
//   int k;
//   double dK;
//   
//   double dLPMF = 0.0;
//   
//   for (k = 0; k < iTrunc; k++) {
//     
//     dK = (double)k;
//     
//     dLPMF += (exp(-dK * dNu) * pow(dK * dNu, dY) * exp(-dMu) * pow(dMu, dK)/(factorial2(dY) * factorial2(dK)));
//   }
//   
//   if(bLog) {
//     dLPMF = log(dLPMF);
//   }
//   
//   return dLPMF;
// }

//[[Rcpp::export]]
double dNeyman2(double dY, double dNu, double dMu, bool bLog = true, int iTrunc = 100) {
  
  int i;
  
  arma::vec vPOI1_log(iTrunc);
  arma::vec vPOI2_log(iTrunc);
  
  for (i = 0; i < iTrunc; i++) {
    
    vPOI1_log(i) = Rf_dpois((double)i, dMu, 1);
    
    if (i == 0) {
      vPOI2_log(i)  = Rf_dpois(dY, dNu * 1e-8, true);
    } else {
      
      vPOI2_log(i)  = Rf_dpois(dY, dNu * (double)i, true);
      
    }
  }
  
  double dLPMF = MixtDensityScale(vPOI1_log, vPOI2_log);
  
  if(!bLog) {
    dLPMF = exp(dLPMF);
  }
  
  return dLPMF;
}

//[[Rcpp::export]]
double dNeyman(double dY, double dNu, double dMu, bool bLog = true, int iTrunc = 100) {
  
  arma::vec vBaz(iTrunc - 1);
  double dLPMF = 0.0;
  
  double dLogFoo = 0.0;
  
  int i;
  double dK;
  double dlK;
  double dlMu= log(dMu);
  double dLogFactY = lfactorial2(dY);
  
  double dStop = dMu - 25.0;
  
  double dPoi_Prop_P = 10;
  
    dLogFoo = 0.0;
    dPoi_Prop_P = 10.0;
    i = 1;
    
    while (dPoi_Prop_P > dStop && i < iTrunc) {
      
      dK = i * 1.0;
      
      dlK = log(dK);
      
      dLogFoo += dlK;
      
      dPoi_Prop_P = dK * dlMu  - dLogFoo;
      
      vBaz(i - 1) = dPoi_Prop_P - dK * dNu + dY * dlK;
      
      i += 1;
    }
    
  dLPMF = - dMu - dLogFactY + dY * log(dNu) + LogSumExp(vBaz.subvec(0, i - 2));
  
  if(!bLog) {
    dLPMF = exp(dLPMF);
  }
  
  return dLPMF;
  
}

//[[Rcpp::export]]
double rNeyman(double dNu, double dMu) {
  
  double dK = Rf_rpois(dMu);
  
  double dY = 0.0;
  
  if (dK > 0) {
    dY = Rf_rpois(dNu * dK);
  }
  return dY;
}
