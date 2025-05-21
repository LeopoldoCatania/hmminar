#include <RcppArmadillo.h>
#include "PoissonBinomial.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;


double const dTwo_Pi = 2.0 * M_PI;


//[[Rcpp::export]]
arma::vec Sim_MixInar(int iT, arma::vec vOmega, arma::vec vLambda, double dAlpha) {
  
  int t;
  int iFoo;
  
  arma::vec vY(iT);
  
  iFoo = rando_index(vOmega);
  vY(0) = Rf_rpois(vLambda(iFoo)/(1.0 - dAlpha));
  
  for (t = 1; t < iT; t++){
    iFoo = rando_index(vOmega);
    
    vY(t) = Rf_rbinom(vY(t - 1), dAlpha) +  Rf_rpois(vLambda(iFoo));
    
  }
  
  return vY;
  
}

//[[Rcpp::export]]
List Sim_MixInar_Seasonal(int iT, int iS, arma::vec vOmega, arma::vec vLambda, double dAlpha, arma::vec vEta) {
  
  int t;
  int s;
  int iFoo;
  
  arma::vec vPsi(iS);
  arma::vec vY(iT);
  
  // compute periodic pattern
  for (s = 0; s < iS; s++) {
    vPsi(s) = exp(vEta(0) * cos(dTwo_Pi * ((s + 1.0)/(iS * 1.0)) - vEta(1) * M_PI));
  }
  
  iFoo = rando_index(vOmega);
  vY(0) = Rf_rpois(vLambda(iFoo)/(1.0 - dAlpha) * vPsi(0));
  
  s = 1;

  for (t = 1; t < iT; t++){
    iFoo = rando_index(vOmega);
    
    vY(t) = Rf_rbinom(vY(t - 1), dAlpha) +  Rf_rpois(vLambda(iFoo) * vPsi(s));
    
    if (s == iS - 1) {
      s = 0;
    } else {
      s += 1;
    }
    
  }
  
  List lOut;
  
  lOut["vY"] = vY;
  lOut["vPsi"] = vPsi;
  
  return lOut;
  
}

//[[Rcpp::export]]
arma::vec Sim_MixInar_Seasonal_Dummy(int iT, arma::vec vOmega, arma::vec vLambda, double dAlpha, arma::vec vBeta,
                                arma::mat mD) {
  
  int t;
  int d;
  int iFoo;
  int iD = mD.n_cols;
  
  arma::vec vY(iT);
  arma::vec vSeason(iT);
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }
  
  
  iFoo = rando_index(vOmega);
  vY(0) = Rf_rpois(vLambda(iFoo)/(1.0 - dAlpha) * accu(vBeta) / (iD * 1.0));
  
  for (t = 1; t < iT; t++){
    
    iFoo = rando_index(vOmega);
    vY(t) = Rf_rbinom(vY(t - 1), dAlpha) +  Rf_rpois(vLambda(iFoo) * vBeta(vSeason(t)));
    
  }
  
  return vY;
  
}

//[[Rcpp::export]]
List Filtering_MixInar_C(arma::vec vY, int iJ, arma::vec vOmega, double dAlpha, arma::vec vLambda, arma::vec vBeta,
                       arma::mat mD) {
  
  // definitions
  int iT = vY.size();
  int iD = mD.n_cols;
  int t;
  int j;
  int d;
  
  arma::vec vMu(iT - 1);
  arma::vec vSigma2(iT - 1);
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  for (d = 0; d < iD; d++) {
    if (mD(0, d) == 1) {
      vSeason(0) = d;
    }
  }
  
  double dEps_mean = 0.0;
  double dEps_var = 0.0;
  
  for (t = 1; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
    
    dEps_mean = 0.0;
    dEps_var  = 0.0;
    
    for (j = 0; j < iJ; j++) {
      dEps_mean += (vOmega(j) * vLambda(j) * vBeta(vSeason(t)));
      dEps_var  += (vOmega(j) * vLambda(j) * vBeta(vSeason(t)) * (1.0 + vLambda(j) * vBeta(vSeason(t))));
    }
    
    dEps_var = dEps_var - pow(dEps_mean, 2.0);
    
    vMu(t - 1) = dAlpha * vY(t - 1) + dEps_mean;
    vSigma2(t - 1) = dAlpha * (1.0 - dAlpha) * vY(t - 1) + dEps_var;
    
  }
  
  List lOut;
  
  lOut["vMu"] = vMu;
  lOut["vSigma2"] = vSigma2;
  
  return lOut;
  
}


//[[Rcpp::export]]
List FilteringMedian_MixInar_C(arma::vec vY, int iJ, arma::vec vOmega, double dAlpha, arma::vec vLambda, arma::vec vBeta,
                         arma::mat mD, int iB = 5000) {
  
  // definitions
  int iT = vY.size();
  int iD = mD.n_cols;
  int t;
  int j;
  int d;
  int b;
  
  arma::vec vMedian(iT - 1);
  arma::vec vSigma2(iT - 1);
  arma::vec vMu(iT - 1);
  arma::vec vSim(iB);
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  for (d = 0; d < iD; d++) {
    if (mD(0, d) == 1) {
      vSeason(0) = d;
    }
  }
  
  double dEps_mean = 0.0;
  double dEps_var = 0.0;
  
  for (t = 1; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
    
    dEps_mean = 0.0;
    dEps_var  = 0.0;
    
    for (j = 0; j < iJ; j++) {
      dEps_mean += (vOmega(j) * vLambda(j) * vBeta(vSeason(t)));
      dEps_var  += (vOmega(j) * vLambda(j) * vBeta(vSeason(t)) * (1.0 + vLambda(j) * vBeta(vSeason(t))));
    }
    
    dEps_var = dEps_var - pow(dEps_mean, 2.0);
    
    for (b = 0; b < iB; b++) {
      vSim(b) = rMixPB(vY(t - 1), dAlpha, vLambda * vBeta(vSeason(t)), vOmega);
    }
    
    vMu(t - 1) = dAlpha * vY(t - 1) + dEps_mean;
    vSigma2(t - 1) = dAlpha * (1.0 - dAlpha) * vY(t - 1) + dEps_var;
    vMedian(t - 1) = median(vSim); 
  }
  
  List lOut;
  
  lOut["vMu"] = vMu;
  lOut["vSigma2"] = vSigma2;
  lOut["vMedian"] = vMedian;
  
  return lOut;
  
}




