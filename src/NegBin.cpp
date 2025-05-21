#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double dNegBin(double dY, double dN, double dP, bool bLog = true) {
  
  NumericVector vY_foo(1);
  vY_foo(0) = dY;
  
  NumericVector foo = lfactorial(vY_foo);
  
  
  double dPMF = Rf_lgammafn(dN + dY) + dY * log(dP/(1.0 + dP)) + dN * log(1.0 - dP/(1.0 + dP)) - 
    foo(0) - Rf_lgammafn(dN);
  
  if (!bLog) {
    dPMF = exp(dPMF);
  }
  
  return dPMF;
  
}


/*  rlogarithmic() is an implementation with automatic selection of
 *  the LS and LK algorithms of:
 *
 *  Kemp, A. W. (1981), Efficient Generation of Logarithmically
 *  Distributed Pseudo-Random Variables, Journal of the Royal
 *  Statistical Society, Series C. Vol. 30, p. 249-253.
 *  URL http://www.jstor.org/stable/2346348
 *
 *  The algorithms are also discussed in chapter 10 of Devroye (1986).
 */

//[[Rcpp::export]]
double rlogarithmic(double p)
{
  if (p < 0 || p > 1) return R_NaN;
  
  /* limiting case as p approaches zero is point mass at one. */
  if (p == 0) return 1.0;
  
  /* Automatic selection between the LS and LK algorithms */
  if (p < 0.95)
  {
    double s = -p/log1p(-p);
    double x = 1.0;
    double u = Rf_runif(0.0, 1.0);//  unif_rand();
    
    while (u > s)
    {
      u -= s;
      x += 1.0;
      s *= p * (x - 1.0)/x;
    }
    
    return(x);
  }
  
  /* else (p >= 0.95) */
  {
    double r = log1p(-p);
    double v = unif_rand();
    
    if (v >= p)       return 1.0;
    
    double u = Rf_runif(0.0, 1.0);//unif_rand();
    double q = -expm1(r * u);
    
    if (v <= (q * q)) return(round(1.0 + log(v)/log(q)));
    if (v <= q)       return(1.0); /* case q^2 < v <= q */
  return(2.0);		       /* case v > q */
  }
}

//[[Rcpp::export]]
double rNegBin(double dN, double dP) {
  
  double dTheta  = dP/(dP + 1.0);
  double dLambda = -dN * log(1.0 - dTheta);
  
  double dM = Rf_rpois(dLambda);
  
  double dY = 0.0;
  
  if (dM > 0) {
    for (int m = 0; m < (int)dM; m++) {
      dY += rlogarithmic(dTheta);
    }
  }
  
  return dY;
  
}


