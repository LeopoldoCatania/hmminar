#include <RcppArmadillo.h>
#include "Utils.h"
#include "Neyman.h"

using namespace arma;
using namespace Rcpp; 

//[[Rcpp::export]]
Rcpp::List MSSumMixPois_Sim(int iT, 
                           arma::vec vMu, arma::vec vNu, 
                           arma::vec vOmega,
                           arma::mat mGamma) {
  
  
  arma::vec vY(iT);
  int iL = mGamma.n_cols;
  
  arma::vec vDelta = getDelta(mGamma, iL);
  
  int t;
  
  arma::vec vState(iT);
  arma::vec vMixComponent(iT);
  
  vState(0) = rando_index(vDelta);
  vMixComponent(0) = rando_index(vOmega);
  
  vY(0) = rNeyman(vNu(vMixComponent(0)), vMu(vState(0)));
  
  for (t = 1; t < iT; t++) {
    vState(t)        = rando_index(mGamma.row(vState(t - 1)).t());
    vMixComponent(t) = rando_index(vOmega);
    vY(t) = rNeyman(vNu(vMixComponent(t)), vMu(vState(t)));
  }
  
  List lOut;
  
  lOut["vState"]      = vState;
  lOut["vMixComponent"]          = vMixComponent;
  lOut["vY"] = vY;
  return lOut;
}

//[[Rcpp::export]]
Rcpp::List MSSumMixPois_Sim_Dummy(int iT, 
                            arma::vec vMu, arma::vec vNu, 
                            arma::vec vOmega,
                            arma::mat mGamma, arma::vec vBeta,
                            arma::mat mD) {
  
  
  arma::vec vY(iT);
  int iL = mGamma.n_cols;
  int d;
  int iD = mD.n_cols;
  int t;
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }
  
  arma::vec vDelta = getDelta(mGamma, iL);
  
  arma::vec vState(iT);
  arma::vec vMixComponent(iT);
  
  vState(0) = rando_index(vDelta);
  vMixComponent(0) = rando_index(vOmega);
  
  vY(0) = rNeyman(vNu(vMixComponent(0)), vMu(vState(0)) * vBeta(vSeason(0)));
  
  for (t = 1; t < iT; t++) {
    vState(t)        = rando_index(mGamma.row(vState(t - 1)).t());
    vMixComponent(t) = rando_index(vOmega);
    vY(t) = rNeyman(vNu(vMixComponent(t)), vMu(vState(t)) * vBeta(vSeason(t)));
  }
  
  List lOut;
  
  lOut["vState"]      = vState;
  lOut["vMixComponent"]          = vMixComponent;
  lOut["vY"] = vY;
  return lOut;
}



