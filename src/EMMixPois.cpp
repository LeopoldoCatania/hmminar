#include <RcppArmadillo.h>
#include "Utils.h"
#include "SafeFunctions.h"

using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
double dPOI(double dY, double dLambda, int iLog) {
  
  double dPoi = Rf_dpois(dY, dLambda, 1);
  
  if (iLog == 0) {
    dPoi = exp(dPoi);
    
    dPoi = DensityCheck(dPoi);
    
  }
  
  if (iLog != 0) {
    
    dPoi = LogDensityCheck(dPoi);
    
  }
  
  return dPoi;
  
}

Rcpp::List StartingValueEM_MM_Pois(arma::vec vY, int iK)
{
  int iT = vY.size();
  int k;
  int t;
  
  arma::vec vLambda(iK);
  double dMu = accu(vY)/(iT * 1.0);
  
  double start = 0.8;
  double end   = 1.2;
  double by    = (end - start)/(iK * 1.0);
  double foo = start - by;
  
  arma::vec vP(iK);
  vP.fill(1.0/(iK * 1.0));
  
  for(k = 0; k < iK; k++)
  {
    vLambda(k)       = dMu * foo;
    foo += by;
  }
  
  // initialize weights
  arma::mat mLLK(iK, iT);
  arma::mat mW(iK, iT);
  arma::vec vLLK(iT);
  
  for (t = 0; t < iT; t++)
  {
    for (k = 0; k < iK; k++)
    {
      mLLK(k, t) = dPOI(vY(t), vLambda(k), 1);
    }
    vLLK(t) = MixtDensityScale(log(vP), mLLK.col(t));
    
    for (k = 0; k < iK; k++)
    {
      mW(k, t) = exp(log(vP(k)) + mLLK(k, t) - vLLK(t));
    }
  }
  for (k = 0; k < iK; k++)
  {
    vP(k) = accu(mW.row(k));
  }
  
  vP = vP/(iT * 1.0);
  
  List out;
  out["vLambda"] = vLambda;
  out["vP"]      = vP;
  
  return out;
}


//[[Rcpp::export]]
List EM_Mixture_Pois(arma::vec vY, int iJ, int maxIter = 1e3, double tol = 1e-8) {
  
  List lStarting    = StartingValueEM_MM_Pois(vY, iJ);
  arma::vec vLambda = lStarting["vLambda"];
  arma::vec vP      = lStarting["vP"];
  
  arma::vec vLambda_Next = vLambda;
  arma::vec vP_Next = vP;
  
  int iT = vY.size();
  
  int iter = 0;
  int j;
  int t;
  
  arma::vec vLLKSeries(maxIter + 1);
  
  double eps = 1.0;
  
  arma::mat mLLK(iJ, iT);
  arma::mat mW(iJ, iT);
  arma::vec vLLK(iT);
  
  for (t = 0; t < iT; t++)
  {
    for (j = 0; j < iJ; j++)
    {
      mLLK(j, t) = dPOI(vY(t), vLambda(j), 1);
    }
    vLLK(t) = MixtDensityScale(log(vP), mLLK.col(t));
    for (j = 0; j < iJ; j++)
    {
      mW(j, t) = exp(log(vP(j)) + mLLK(j, t) - vLLK(t));
    }
  }
  vLLKSeries(0) = accu(vLLK);
  while(eps > tol && iter<maxIter)
  {
    mLLK.zeros();
    mW.zeros();
    vLLK.zeros();
    vLambda_Next.zeros();
    for (t = 0; t < iT; t++)
    {
      for (j = 0; j < iJ; j++)
      {
        mLLK(j, t) = dPOI(vY(t), vLambda(j), 1);
      }
      vLLK(t) = MixtDensityScale(log(vP), mLLK.col(t));
      for (j = 0; j < iJ; j++)
      {
        mW(j, t) = exp(log(vP(j)) + mLLK(j, t) - vLLK(t));
        
        vLambda_Next(j) += mW(j, t) * vY(t);
      }
    }
    for (j = 0; j < iJ; j++)
    {
      vP_Next(j) = accu(mW.row(j));
      vLambda_Next(j) = vLambda_Next(j)/vP_Next(j);
    }
    
    vP_Next = vP_Next/(iT * 1.0);
    //Store the llk
    vLLKSeries(iter) = accu(vLLK);
    iter += 1;
    if(iter>10)  {
      eps = abs3((vLLKSeries(iter - 1) - vLLKSeries(iter-2))/(vLLKSeries(iter-2) + 1.0));
    }
    
    //Update Parameters
    vLambda    = vLambda_Next;
    vP     = vP_Next;
  }
  
  // Decoding
  arma::vec vDecoding(iT);
  for (t = 0; t < iT; t++)
  {
    vDecoding(t) = WhichMax(log(vP) + mLLK.col(t) - vLLK(t)); 
  }
  
  List EMOut;
  EMOut["mW"]    = mW;
  vLLKSeries = vLLKSeries.subvec(0, iter - 2);
  EMOut["vLLKSeries"] = vLLKSeries;
  EMOut["mLLK"]       = mLLK;
  EMOut["vDecoding"]  = vDecoding;
  EMOut["eps"]        = eps;
  EMOut["iter"]       = iter;
  EMOut["vLambda"]    = vLambda;
  EMOut["vY"]         = vY;
  EMOut["vP"]         = vP;
  
  return EMOut;
}

