#include <RcppArmadillo.h>
#include "PoissonBinomial.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;
//[[Rcpp::export]]
List EMInar(arma::vec vY, arma::vec vLogFactorial_K, double dTol = 1e-4, int iMaxiter = 1000) {
  
  int iT = vY.size();
  int t;
  
  double dMean = accu(vY)/(iT * 1.0);
  
  double dAlpha = accu((vY.subvec(0, iT - 2) - dMean) % (vY.subvec(1, iT - 1) - dMean))/accu(pow(vY - dMean, 2.0));
  
  if (dAlpha < 0.0) {
    dAlpha = 0.1;
  }
  
  double dLambda = dMean * (1.0 - dAlpha);
  
  double dAlpha_Next = dAlpha;
  double dLambda_Next = dLambda;
  
  arma::vec LLKSeries(iMaxiter);
  arma::vec vEps(iT);
  arma::vec vLLK(iT);
  
  vEps.zeros();
  vLLK.zeros();
  
  double dSum_t = accu(vY) - vY(0);
  double dSum_tm1 = dSum_t - vY(iT - 1) + vY(0);
  
  double dLLK = 0.0;
  
  for (t = 1; t < iT; t++) {
    vLLK(t) = dPB(vY(t), vY(t - 1), dAlpha, dLambda, vLogFactorial_K, true);
    dLLK   += vLLK(t);
  }
  
  LLKSeries(0) = dLLK;
  
  int iter = 1;
  double dEps = 10;
  
  while (dEps > dTol && iter < iMaxiter) {
    
    //E step
    for (t = 1; t < iT; t++) {
      vEps(t) = exp(log(dLambda) + dPB(vY(t) - 1.0, vY(t - 1), dAlpha, dLambda, vLogFactorial_K, true) - vLLK(t));
    }
    
    //M Step
    dAlpha_Next  = (dSum_t - accu(vEps))/dSum_tm1;
    dLambda_Next = accu(vEps)/(iT * 1.0 - 1.0);
    
    //update Parameters
    dLambda = dLambda_Next;
    dAlpha  = dAlpha_Next;
    
    dLLK = 0.0;
    for (t = 1; t < iT; t++) {
      vLLK(t) = dPB(vY(t), vY(t - 1), dAlpha, dLambda, vLogFactorial_K, true);
      dLLK   += vLLK(t);
    }
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["dLambda"] = dLambda;
  lOut["dAlpha"] = dAlpha;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  
  return lOut;
  
}

//[[Rcpp::export]]
List EMMixInar(arma::vec vY, int iJ, arma::vec vOmega, double dAlpha, arma::vec vLambda,
               arma::vec vLogFactorial_K, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  int t;
  int j;
  
  arma::vec LLKSeries(iMaxiter);
  arma::mat mQ = zeros(iJ, iT);
  arma::mat mL = zeros(iJ, iT);
  arma::mat mLLK = zeros(iJ, iT);
  arma::vec vLLK(iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  double dAlpha_Next  = dAlpha;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  
  double dSum_t = accu(vY) - vY(0);
  double dSum_tm1 = dSum_t - vY(iT - 1) + vY(0);
  
  for (t = 1; t < iT; t++) {
    for (j = 0; j < iJ; j++) {
      mLLK(j, t) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j), vLogFactorial_K, true);
    }
    vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
    dLLK   += vLLK(t);
  }
  
  LLKSeries(0) = dLLK;
  
  int iter = 1;
  arma::vec vFoo(iJ);
  
  while (dEps > dTol && iter < iMaxiter) {
    //E step
    for (t = 1; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        mQ(j, t) = exp(log(vOmega(j)) + mLLK(j, t) - vLLK(t));
        mL(j, t) = exp(log(vLambda(j)) + dPB(vY(t) - 1.0, vY(t - 1), dAlpha, vLambda(j), vLogFactorial_K, true) - 
          mLLK(j, t));
      }
    }
    
    //M Step
    dAlpha_Next  = (dSum_t - accu(mQ % mL))/dSum_tm1;
    for (j = 0; j < iJ; j++) {
      vOmega_Next(j)  = accu(mQ.row(j))/(iT * 1.0 - 1.0);
      vLambda_Next(j) = accu(mQ.row(j) % mL.row(j))/accu(mQ.row(j));
    }
    
    //update Parameters
    vOmega  = vOmega_Next;
    vLambda = vLambda_Next;
    dAlpha  = dAlpha_Next;
    
    dLLK = 0.0;
    for (t = 1; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        mLLK(j, t) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j), vLogFactorial_K, true);
      }
      vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
      dLLK += vLLK(t);
    }
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["dAlpha"] = dAlpha;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  
  return lOut;
  
}

//[[Rcpp::export]]
List EMMixInar_Seasonal(arma::vec vY, int iJ, arma::vec vOmega, double dAlpha, arma::vec vLambda, arma::vec vBeta,
                        arma::vec vLogFactorial_K, arma::mat mD, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  double dY0 = vY(0);
  vY = vY.subvec(1, iT - 1);
  iT = iT - 1;
  
  int iD = mD.n_cols;
  int t;
  int j;
  int d;
  
  arma::vec LLKSeries(iMaxiter);
  arma::mat mQ = zeros(iJ, iT);
  arma::mat mL = zeros(iJ, iT);
  arma::mat mLLK = zeros(iJ, iT);
  arma::vec vLLK(iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }
  
  bool bSeason = false;
  if (iD > 1) {
    bSeason = true;
  }
  if (iD == 1) {
    if (accu(mD.col(0)) < iT) {
      bSeason = true;
    }
  }
  
  
  double dAlpha_Next  = dAlpha;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vBeta_Next = vBeta;
  
  arma::vec vBeta_Num_foo(iD);
  arma::vec vBeta_Den_foo(iD);
  
  arma::vec vLambda_Num_foo(iJ);
  arma::vec vLambda_Den_foo(iJ);
  
  int iter = 0;
  
  double dSum_t = accu(vY) - vY(0);
  double dSum_tm1 = dSum_t - vY(iT - 1) + vY(0);
  
  while (dEps > dTol && iter < iMaxiter) {
    
    dLLK = 0.0;
    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          mLLK(j, t) = dPB(vY(t), dY0, dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);
        } else {
          mLLK(j, t) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);
        }
     
      }
      vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
      dLLK   += vLLK(t);
    }
    
    LLKSeries(iter) = dLLK;
    
    //E step
    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          mL(j, t) = exp(log(vLambda(j) * vBeta(vSeason(t))) + dPB(vY(t) - 1.0, dY0, dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true) - 
            mLLK(j, t)); 
        } else {
          mL(j, t) = exp(log(vLambda(j) * vBeta(vSeason(t))) + dPB(vY(t) - 1.0, vY(t - 1), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true) - 
            mLLK(j, t)); 
        }
        mQ(j, t) = exp(log(vOmega(j)) + mLLK(j, t) - vLLK(t));
      }
    }
    
    //M Step
    dAlpha_Next  = (dSum_t - accu(mQ % mL))/dSum_tm1;
    
    vLambda_Num_foo.zeros();
    vLambda_Den_foo.zeros();
    
    for (j = 0; j < iJ; j++) {
      vOmega_Next(j)  = accu(mQ.row(j))/(iT * 1.0 - 1.0);
      
      for (t = 0; t < iT; t++) {
        vLambda_Num_foo(j) += (mQ(j, t) * mL(j, t));
        vLambda_Den_foo(j) += (mQ(j, t) * vBeta(vSeason(t)));
      }
      vLambda_Next(j) = vLambda_Num_foo(j) / vLambda_Den_foo(j);
    }
    
    
    vBeta_Num_foo.zeros();
    vBeta_Den_foo.zeros();
    
    if (bSeason) {
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          vBeta_Num_foo(vSeason(t)) += (mQ(j, t) * mL(j, t));
          vBeta_Den_foo(vSeason(t)) += (mQ(j, t) * vLambda(j));
        }
      }
      
      for (d = 0; d < iD; d++) {
        vBeta_Next(d) = vBeta_Num_foo(d)/vBeta_Den_foo(d);
      }
    }
    
    //update Parameters
    vOmega  = vOmega_Next;
    vLambda = vLambda_Next;
    dAlpha  = dAlpha_Next;
    vBeta   = vBeta_Next;
    
    if(iter > 10) {
      
      dEps = abs3(((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1)));
      
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["dAlpha"] = dAlpha;
  lOut["vBeta"] = vBeta;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  lOut["mLLK"] = mLLK;
  lOut["mL"] = mL;
  lOut["mQ"] = mQ;
  
  return lOut;
  
}

//[[Rcpp::export]]
List EMMixInar_Seasonal2(arma::vec vY, int iJ, arma::vec vOmega, double dAlpha, arma::vec vLambda, arma::vec vBeta,
                        arma::vec vLogFactorial_K, arma::mat mD, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  int iD = mD.n_cols;
  int t;
  int j;
  int d;
  
  arma::vec LLKSeries(iMaxiter);
  arma::mat mQ = zeros(iJ, iT);
  arma::mat mL = zeros(iJ, iT);
  arma::mat mLLK = zeros(iJ, iT);
  arma::vec vLLK(iT);
  
  double dLLK = 0.0;
  double dEps = 10;
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }
  
  bool bSeason = false;
  if (iD > 1) {
    bSeason = true;
  }
  if (iD == 1) {
    if (accu(mD.col(0)) < iT) {
      bSeason = true;
    }
  }
  
  
  double dAlpha_Next  = dAlpha;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vBeta_Next = vBeta;
  
  arma::vec vBeta_Num_foo(iD);
  arma::vec vBeta_Den_foo(iD);
  
  arma::vec vLambda_Num_foo(iJ);
  arma::vec vLambda_Den_foo(iJ);
  
  double dSum_t = accu(vY) - vY(0);
  double dSum_tm1 = dSum_t - vY(iT - 1) + vY(0);
  
  int iter = 0;
  
  for (t = 1; t < iT; t++) {
    for (j = 0; j < iJ; j++) {
      mLLK(j, t) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);
    }
    vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
    dLLK   += vLLK(t);
  }
  
  LLKSeries(0) = dLLK;
  
  while (dEps > dTol && iter < iMaxiter) {
    
    //E step
    for (t = 1; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        mQ(j, t) = exp(log(vOmega(j)) + mLLK(j, t) - vLLK(t));
        mL(j, t) = exp(log(vLambda(j) * vBeta(vSeason(t))) + dPB(vY(t) - 1.0, vY(t - 1), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true) - 
          mLLK(j, t));
      }
    }
    
    //M Step
    dAlpha_Next  = (dSum_t - accu(mQ % mL))/dSum_tm1;
    
    vLambda_Num_foo.zeros();
    vLambda_Den_foo.zeros();
    
    for (j = 0; j < iJ; j++) {
      vOmega_Next(j)  = accu(mQ.row(j))/(iT * 1.0 - 1.0);
      
      for (t = 1; t < iT; t++) {
        // vLambda_Num_foo(j) += (mQ(j, t) * mL(j, t) * vBeta(vSeason(t)));
        vLambda_Num_foo(j) += (mQ(j, t) * mL(j, t));
        vLambda_Den_foo(j) += (mQ(j, t) * vBeta(vSeason(t)));
      }
      
      // vLambda_Next(j) = accu(mQ.row(j) % mL.row(j))/accu(mQ.row(j));
      vLambda_Next(j) = vLambda_Num_foo(j) / vLambda_Den_foo(j);
      
    }
    
    
    vBeta_Num_foo.zeros();
    vBeta_Den_foo.zeros();
    
    if (bSeason) {
      for (t = 1; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          vBeta_Num_foo(vSeason(t)) += (mQ(j, t) * mL(j, t));
          vBeta_Den_foo(vSeason(t)) += (mQ(j, t) * vLambda(j));
        }
      }
      
      for (d = 0; d < iD; d++) {
        vBeta_Next(d) = vBeta_Num_foo(d)/vBeta_Den_foo(d);
      }
    }
    
    //update Parameters
    vOmega  = vOmega_Next;
    vLambda = vLambda_Next;
    dAlpha  = dAlpha_Next;
    vBeta   = vBeta_Next;
    
    dLLK = 0.0;
    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          mLLK(j, t) = dPB(vY(t), vY(t), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);
        } else {
          mLLK(j, t) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);  
        }
        
      }
      vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
      dLLK += vLLK(t);
    }
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      
      dEps = abs3(((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1)));
      
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["dAlpha"] = dAlpha;
  lOut["vBeta"] = vBeta;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  lOut["mLLK"] = mLLK;
  lOut["mL"] = mL;
  lOut["mQ"] = mQ;
  
  return lOut;
  
}


//[[Rcpp::export]]
List EMMixInarP_Seasonal(arma::vec vY, int iJ, int iP, arma::vec vOmega, double dAlpha, arma::vec vLambda, arma::vec vBeta,
                         arma::mat mD, double dTol = 1e-4, int iMaxiter = 1000) {
  
  // definitions
  int iT = vY.size();
  int iD = mD.n_cols;
  int t;
  int j;
  int d;
  
  arma::vec LLKSeries(iMaxiter);
  arma::mat mQ = zeros(iJ, iT);
  arma::mat mL = zeros(iJ, iT);
  arma::mat mLLK = zeros(iJ, iT);
  arma::vec vLLK(iT);
  arma::vec vY_Sumlag(iT);
  vY_Sumlag.zeros();
  
  double dLLK = 0.0;
  double dEps = 10;
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
    // zero based we start from t == p
    if (t >= iP) {
      vY_Sumlag(t) = accu(vY.subvec(t - iP, t - 1));
    }
  }
  
  int iMax = (int) max(vY_Sumlag);
  arma::vec vLogFactorial_K(iMax + 1);
  
  vLogFactorial_K(0) = 0;
  
  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }
  
  double dAlpha_Next  = dAlpha;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vBeta_Next = vBeta;
  
  arma::vec vBeta_Num_foo(iD);
  arma::vec vBeta_Den_foo(iD);
  
  arma::vec vLambda_Num_foo(iJ);
  arma::vec vLambda_Den_foo(iJ);
  
  double dSum_t = accu(vY) - accu(vY.subvec(0, iP - 1));
  double dSum_tm1_p = accu(vY_Sumlag);
  
  for (t = iP; t < iT; t++) {
    for (j = 0; j < iJ; j++) {
      mLLK(j, t) = dPB(vY(t), vY_Sumlag(t), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);
    }
    vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
    dLLK   += vLLK(t);
  }
  
  LLKSeries(0) = dLLK;
  
  int iter = 1;
  
  while (dEps > dTol && iter < iMaxiter) {
    //E step
    for (t = iP; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        mQ(j, t) = exp(log(vOmega(j)) + mLLK(j, t) - vLLK(t));
        mL(j, t) = exp(log(vLambda(j) * vBeta(vSeason(t))) + dPB(vY(t) - 1.0, vY_Sumlag(t), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true) - 
          mLLK(j, t));
      }
    }
    
    //M Step
    dAlpha_Next  = (dSum_t - accu(mQ % mL))/dSum_tm1_p;
    
    vLambda_Num_foo.zeros();
    vLambda_Den_foo.zeros();
    
    for (j = 0; j < iJ; j++) {
      vOmega_Next(j)  = accu(mQ.row(j))/(iT * 1.0 - 1.0);
      
      for (t = iP; t < iT; t++) {
        vLambda_Num_foo(j) += (mQ(j, t) * mL(j, t));
        vLambda_Den_foo(j) += (mQ(j, t) * vBeta(vSeason(t)));
      }
      
      vLambda_Next(j) = vLambda_Num_foo(j) / vLambda_Den_foo(j);
      
    }
    
    vBeta_Num_foo.zeros();
    vBeta_Den_foo.zeros();
    
    for (t = iP; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        vBeta_Num_foo(vSeason(t)) += (mQ(j, t) * mL(j, t));
        vBeta_Den_foo(vSeason(t)) += (mQ(j, t) * vLambda(j));
      }
    }
    
    for (d = 0; d < iD; d++) {
      vBeta_Next(d) = vBeta_Num_foo(d)/vBeta_Den_foo(d);
    }
    
    //update Parameters
    vOmega  = vOmega_Next;
    vLambda = vLambda_Next;
    dAlpha  = dAlpha_Next;
    vBeta   = vBeta_Next;
    
    dLLK = 0.0;
    for (t = iP; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        mLLK(j, t) = dPB(vY(t), vY_Sumlag(t), dAlpha, vLambda(j) * vBeta(vSeason(t)), vLogFactorial_K, true);
      }
      vLLK(t) = MixtDensityScale(log(vOmega), mLLK.col(t));
      dLLK += vLLK(t);
    }
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  List lOut;
  
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["dAlpha"] = dAlpha;
  lOut["vBeta"] = vBeta;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  
  return lOut;
  
}
