#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

const double dLowerProb = 1e-7;
const double dUpperProb = 1 - dLowerProb;

const double dZeroBound = 1e-10;
const double dLowerBound = -1e5;
// const double dUpperBound = 1e5;

const double dUpperBoundSigma = 115.0;

const double dUpperBoundBeta = 3.0;

//[[Rcpp::export]]
arma::mat GammaCheck(arma::mat mGamma, int iJ) {

  int i;
  int j;

  for (i = 0; i < iJ; i++) {
    for (j = 0; j < iJ; j++) {
      if (mGamma(i, j) < dLowerProb) {
        mGamma(i, j) = dLowerProb;
      }
      if (mGamma(i, j) > dUpperProb) {
        mGamma(i, j) = dUpperProb;
      }
    }
  }

  for (i = 0; i < iJ; i++) {
    mGamma.row(i) = mGamma.row(i) / accu(mGamma.row(i));
  }

  return mGamma;

}

//[[Rcpp::export]]
arma::mat OmegaCheck(arma::mat mOmega, int iK, int iJ) {

  int j;
  int k;

  for (j = 0; j < iJ; j++) {
    for (k = 0; k < iK; k++) {
      if (mOmega(k, j) < dLowerProb) {
        mOmega(k, j) = dLowerProb;
      }
      if (mOmega(k, j) > dUpperProb) {
        mOmega(k, j) = dUpperProb;
      }
    }
    mOmega.col(j) = mOmega.col(j)/accu(mOmega.col(j));
  }

  return mOmega;
}

arma::vec DeltaCheck(arma::vec vDelta, int iK) {

  int i;
  for (i = 0; i < iK; i++) {

    if (vDelta(i) < dLowerProb) {
      vDelta(i) = dLowerProb;
    }
    if (vDelta(i) > dUpperProb) {
      vDelta(i) = dUpperProb;
    }

  }

  vDelta = vDelta / accu(vDelta);

  return vDelta;
}

double LogDensityCheck(double dLogDensity) {

  if (std::isnan(dLogDensity)) {
    dLogDensity = dLowerBound;
  }

  if (dLogDensity < dLowerBound) {

    dLogDensity = dLowerBound;

  }

  return dLogDensity;

}

double DensityCheck(double dDensity) {

  if (std::isnan(dDensity)) {
    dDensity = dZeroBound;
  }

  if (dDensity < dZeroBound) {
    dDensity = dZeroBound;
  }

  return dDensity;

}

arma::mat SigmaCheck(arma::mat mSigma, int iN, int iK) {


  int n;
  int k;

  for (n = 0; n < iN; n++) {
    for (k = 0; k < iK; k++) {
      if (mSigma(n, k) < dZeroBound) {
        mSigma(n, k) = dZeroBound;
      }
      if (mSigma(n, k) > dUpperBoundSigma) {
        mSigma(n, k) = dUpperBoundSigma;
      }
    }
  }

  return mSigma;
}

arma::mat mKappaCheck(arma::mat mKappa, int iN, int iL) {


  int n;
  int l;

  for (n = 0; n < iN; n++) {
    for (l = 0; l < iL; l++) {
      if (mKappa(n, l) < dLowerProb) {
        mKappa(n, l) = dLowerProb;
      }
      if (mKappa(n, l) > dUpperProb) {
        mKappa(n, l) = dUpperProb;
      }
    }
  }

  return mKappa;

}

arma::mat BetaCheck(arma::mat mBeta, int iN, int iD) {

  int n;
  int d;

  for (n = 0; n < iN; n++) {
    for (d = 0; d < iD; d++) {
      if (mBeta(n, d) < dZeroBound) {
        mBeta(n, d) = dZeroBound;
      }
      if (mBeta(n, d) > dUpperBoundBeta) {
        mBeta(n, d) = dUpperBoundBeta;
      }
    }
  }

  return mBeta;
}
