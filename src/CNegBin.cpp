#include <RcppArmadillo.h>
#include "NegBin.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double dCNegBin(double dY, double dN, double dP, double dLambda, bool bLog = true, int iTrunc = 50) {

  int i;

  arma::vec vPOI_log(iTrunc);
  arma::vec vNB_log(iTrunc);

  for (i = 0; i < iTrunc; i++) {

    vPOI_log(i) = Rf_dpois((double)i, dLambda, 1);

    if (i == 0) {
      vNB_log(i)  = dNegBin(dY, dN * 1e-8, dP, true);
    } else {

    vNB_log(i)  = dNegBin(dY, dN * (double)i, dP, true);

    }
  }

  double dLPMF = MixtDensityScale(vPOI_log, vNB_log);

  if(!bLog) {
    dLPMF = exp(dLPMF);
  }

  return dLPMF;

}

//[[Rcpp::export]]
double dLLK_cNegBin(arma::vec vY, double dN, double dP, double dLambda) {
  double dLLK = 0.0;
  
  int iT = vY.size();
  
  for (int i = 0; i < iT; i++) {
    dLLK += dCNegBin(vY(i), dN, dP, dLambda);
  }
  
  return dLLK;
  
}

// //[[Rcpp::export]]
// double dCNegBin(double dY, double dN, double dP, double dLambda, bool bLog = true) {
//   
//   double dLPMF = 0.0;
//   arma::vec vFoo((int)dY + 1);
//   
//   int iY;
//   double dX;
//   
//   if (dY < 1) {
//     dLPMF = dLambda + pow(1.0 + dP, -1.0 * dN) - 1.0;
//   } else {
//     
//     dX = dY - 1.0;
//     
//     for (iY = 0; iY <= dY; iY++) {
//       vFoo(iY) = Rf_lchoose(dN + dX - iY*1.0, dN) + (dX - iY*1.0) * log(dP/(1.0 + dP)) * 
//         dCNegBin((double)iY, dN, dP, dLambda, true);
//     }
//     
//     dLPMF = LogSumExp(vFoo);
//     
//     dLPMF += (log(dLambda) + log(dN) + log(dP) - log(dY) - (dN + 1) * log(1.0 + dP));
//     
//   }
//   
//   if(!bLog) {
//     dLPMF = exp(dLPMF);
//   }
//   
//   return dLPMF;
//   
// }

