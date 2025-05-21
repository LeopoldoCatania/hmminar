#include <RcppArmadillo.h>
#include "math.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
int StirlingNumber(int iN, int iM) {
  
  int iS;
  
  if (iM > iN) {
    iS = 0;
  } else if ((iN == 0) && (iM == 0)) {
    iS = 1;
  } else if ((iN == 0) && (iM > 0)) {
    iS = 0;
  } else if ((iM == 0) && (iN > 0)) {
    iS = 0;
  } else {
    iS = StirlingNumber(iN - 1, iM - 1) - (iN - 1) * StirlingNumber(iN - 1, iM);
  }
  
  return iS;
  
}


double AbsStirlingNumber(double dN, double dM) {
  
  int iS = StirlingNumber((int)dN, (int)dM);
  
  double dAbsS = (double)iS;
  
  if (dAbsS < 0) {
    dAbsS = -1.0 * dAbsS;
  }
  
  return dAbsS;
  
}

//[[Rcpp::export]]
double dStirling(double dY, double dM, double dTheta, bool bLog = true) {
  
  NumericVector vY_foo(2);
  vY_foo(0) = dY;
  vY_foo(1) = dM;
  
  double dAlpha = -1.0 / log(1.0 - dTheta);
  
  NumericVector foo = lfactorial(vY_foo);
  
  double dAbsStirling = AbsStirlingNumber(dY, dM);
  
  double dPMF = foo(1) + dY * log(dTheta) + log(dAbsStirling) + dM * log(dAlpha) - foo(0);
  
  if (!bLog) {
    
    dPMF = exp(dPMF);
    
  }
  
  return dPMF;
  
}

