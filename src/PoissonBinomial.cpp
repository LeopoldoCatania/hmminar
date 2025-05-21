#include <RcppArmadillo.h>
#include "Utils.h"
//#include <omp.h>

using namespace Rcpp;
using namespace arma;

const double dZeroBound = 1e-15;
double const dTwo_Pi = 2.0 * M_PI;

////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// FOR INAR  ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
double dPB2(double dY, double dX, double dAlpha, double dLambda,
            arma::vec vLogFactorial_K, bool bLog = true) {

  int k = 0;
  int iY = (int)dY;
  int iX = (int)dX;

  int iMax = iY;
  if (iY > iX) {
    iMax = iX;
  }

  double dLPMF = -INFINITY;
  double dFoo = 0.0;
  arma::vec vFoo(2);

  double dLogAlpha = log(dAlpha);
  double dLog1mAlpha = log(1.0 - dAlpha);
  double dLogLambda = log(dLambda);

  int iStart = 0.0;

  double dk;

  if (dY >= 0) {

    for (k = iStart; k <= iMax; k++) {

      dk = (double)k;

      // dFoo = - dLambda + dk * dLogLambda - vLogFactorial_K(k) + vLogFactorial_K(iX)
      //   - vLogFactorial_K(iY - k) - vLogFactorial_K(iX - iY + k) + (dY - dk) * dLogAlpha +
      //     (dX - dY + dk) * dLog1mAlpha;

      dFoo = vLogFactorial_K(iX) - vLogFactorial_K(k) - vLogFactorial_K(iX - k) + dk * dLogAlpha +
        (dX - dk) * dLog1mAlpha - dLambda + (dY - dk) * dLogLambda - vLogFactorial_K(iY - k);

      vFoo(0) = dLPMF;
      vFoo(1) = dFoo;

      dLPMF = LogSumExp(vFoo);

    }
  }

  if(!bLog) {
    dLPMF = exp(dLPMF);
  }

  return dLPMF;

}

//[[Rcpp::export]]
double dPB(double dY, double dX, double dAlpha, double dLambda,
           arma::vec vLogFactorial_K, bool bLog = true) {

  if(dAlpha < dZeroBound) {
    dAlpha = dZeroBound;
  }
  if (dAlpha > 1.0 - dZeroBound) {
    dAlpha = 1.0 - dZeroBound;
  }
  if(dLambda < dZeroBound) {
    dLambda = dZeroBound;
  }

  int k = 0;
  int iY = (int)dY;
  int iX = (int)dX;

  int iMax = iY;

  double dLPMF = -INFINITY;
  double dFoo = 0.0;
  arma::vec vFoo(2);

  double dLogAlpha = log(dAlpha);
  double dLog1mAlpha = log(1.0 - dAlpha);
  double dLogLambda = log(dLambda);

  int iStart = iY - iX;
  if (iStart < 0) {
    iStart = 0;
  }

  double dk;

  if (dY >= 0) {

    if (dX >= 1) {

      for (k = iStart; k <= iMax; k++) {

        dk = (double)k;

        dFoo = - dLambda + dk * dLogLambda - vLogFactorial_K(k) + vLogFactorial_K(iX)
          - vLogFactorial_K(iY - k) - vLogFactorial_K(iX - iY + k) + (dY - dk) * dLogAlpha +
            (dX - dY + dk) * dLog1mAlpha;

        vFoo(0) = dLPMF;
        vFoo(1) = dFoo;

        dLPMF = LogSumExp(vFoo);

      }

    } else {

      if(dY>0) {
        dLPMF = -dLambda + dY * dLogLambda - vLogFactorial_K((int)dY);
      } else {
        dLPMF = -dLambda;
      }
    }
  }
  if(!bLog) {
    dLPMF = exp(dLPMF);
  }


  return dLPMF;

}


//[[Rcpp::export]]
double dMixPB(double dY, double dX, double dAlpha, arma::vec vLambda,
              arma::vec vOmega,
              arma::vec vLogFactorial_K, bool bLog = true){

  int iJ = vOmega.size();
  int j;

  arma::vec vLP(iJ);

  for (j = 0; j < iJ; j++) {
    vLP(j) = log(vOmega(j)) + dPB(dY, dX, dAlpha, vLambda(j), vLogFactorial_K, true);
  }

  double dP = LogSumExp(vLP);

  if (!bLog) {
    dP = exp(dP);
  }

  return dP;
}

//[[Rcpp::export]]
double pPB(double dY, double dX, double dAlpha, double dLambda,
           arma::vec vLogFactorial_K, bool bLog = true) {

  int k = 0;
  int iY = (int)dY;
  int iX = (int)dX;

  int iEnd = iY;

  if (iX < iY) {
    iEnd = iX;
  }

  double dLCMF = -INFINITY;
  double dFoo = 0.0;
  arma::vec vFoo(2);

  double dLogAlpha = log(dAlpha);
  double d1mLogAlpha = log(1.0 - dAlpha);

  double dk;
  double dPois = 0.0;

  if (dY >= 0) {

    for (k = 0; k <= iEnd; k++) {

      dk = (double)k;

      dPois = Rf_ppois(dY - dk, dLambda, 1, 1);

      dFoo = dPois + vLogFactorial_K(iX) - vLogFactorial_K(k) - vLogFactorial_K(iX - k) + dk * dLogAlpha + (dX - dk) * d1mLogAlpha;

      vFoo(0) = dLCMF;
      vFoo(1) = dFoo;

      dLCMF = LogSumExp(vFoo);

    }
  }

  if(!bLog) {
    dLCMF = exp(dLCMF);
  }

  return dLCMF;

}


//[[Rcpp::export]]
double pMixPB(double dY, double dX, double dAlpha, arma::vec vLambda,
              arma::vec vOmega,
              arma::vec vLogFactorial_K, bool bLog = true){

  int iJ = vOmega.size();
  int j;

  arma::vec vLP(iJ);

  for (j = 0; j < iJ; j++) {
    vLP(j) = log(vOmega(j)) + pPB(dY, dX, dAlpha, vLambda(j), vLogFactorial_K, true);
  }

  double dP = LogSumExp(vLP);

  if (!bLog) {
    dP = exp(dP);
  }

  return dP;
}

double rPB(double dX, double dAlpha, double dLambda) {

  double dY = Rf_rbinom(dX, dAlpha) + Rf_rpois(dLambda);

  return dY;

}

double rMixPB(double dX, double dAlpha, arma::vec vLambda, arma::vec vOmega) {

  int indicator = rando_index(vOmega);
  double dY = rPB(dX, dAlpha, vLambda(indicator));
  return dY;

}

//[[Rcpp::export]]
double dLLK_MixPB(arma::vec vY, arma::vec vLambda, arma::vec vOmega,
                  arma::vec vLogFactorial_K, double dAlpha, int iT, int iJ) {

  double dLLK = 0.0;

  int t;
  int j;

  arma::vec vD(iJ);

  for (t = 1; t < iT; t++) {

    for (j = 0; j < iJ; j++) {
      vD(j) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j),
         vLogFactorial_K, true);
    }
    dLLK += MixtDensityScale(log(vOmega), vD);
  }

  return dLLK;

}

//[[Rcpp::export]]
List dLLK_MixPB_seasonal(arma::vec vY, arma::vec vLambda, arma::vec vOmega, arma::vec vEta,
                         arma::vec vLogFactorial_K, double dAlpha, int iT, int iJ, int iS) {

  double dLLK = 0.0;
  arma::vec vD(iJ);
  arma::vec vPsi(iS);

  int t;
  int s;
  int j;

  // compute periodic pattern
  for (s = 0; s < iS; s++) {
    vPsi(s) = exp(vEta(0) * cos(dTwo_Pi * ((s + 1.0)/(iS * 1.0)) - vEta(1) * M_PI));
  }

  s = 1;

  for (t = 1; t < iT; t++) {

    for (j = 0; j < iJ; j++) {
      vD(j) = dPB(vY(t), vY(t - 1), dAlpha, vLambda(j) * vPsi(s),
         vLogFactorial_K, true);
    }

    dLLK += MixtDensityScale(log(vOmega), vD);

    if (s == iS - 1) {
      s = 0;
    } else {
      s += 1;
    }


  }

  List lOut;

  lOut["dLLK"] = dLLK;
  lOut["vPsi"] = vPsi;

  return lOut;

}






////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// THREE COMPONENTS MIXTURE ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
double dSumThreeBinom(double dY, arma::vec vK, arma::vec vPi, bool bLog = true) {

  double dK = sum(vK);

  int j;
  int l;
  double dj;
  double dl;

  double dLPMF = -INFINITY;

  double dFoo = 0.0;
  arma::vec vFoo(2);

  if (dY <= dK) {

    for (j = 0; j <= dY; j++) {
      dj = (double)j;
      for (l = 0; l <= dY - j; l++) {
        dl = (double)l;

        // if ((dY - dj - dl <= vK(0)) & (dl <= vK(1)) & (dj <= vK(2))) {
        if ((dY - dj - dl <= vK(0)) && (dl <= vK(1)) && (dj <= vK(2))) {

          dFoo = Rf_dbinom(dY - dj - dl, vK(0), vPi(0), 1) +
            Rf_dbinom(dl, vK(1), vPi(1), 1) + Rf_dbinom(dj, vK(2), vPi(2), 1);

          vFoo(0) = dLPMF;
          vFoo(1) = dFoo;

          dLPMF = LogSumExp(vFoo);

        }

      }
    }

  }

  if (!bLog) {
    dLPMF = exp(dLPMF);
  }

  return dLPMF;

}

//[[Rcpp::export]]
double dSumThreeBinom_Eff(double dY, arma::vec vK, arma::vec vPi, arma::vec vLogFactorial_K, bool bLog = true) {

  double dK = sum(vK);

  int j;
  int l;
  double dj;
  double dl;
  int iY = (int)dY;

  int iK1 = (int)vK(0);
  int iK2 = (int)vK(1);
  int iK3 = (int)vK(2);

  double dLPMF = -INFINITY;

  double dFoo = 0.0;
  arma::vec vFoo(2);
  arma::vec vPi_log = log(vPi);
  arma::vec v1mPi_log = log(1.0 - vPi);

  if (dY <= dK) {

    double dlogFactorial_K1 = vLogFactorial_K(iK1);
    double dlogFactorial_K2 = vLogFactorial_K(iK2);
    double dlogFactorial_K3 = vLogFactorial_K(iK3);

    double dSumLogFac =  dlogFactorial_K1 + dlogFactorial_K2 + dlogFactorial_K3;

    for (j = 0; j <= dY; j++) {
      dj = (double)j;
      for (l = 0; l <= dY - j; l++) {
        dl = (double)l;

        // if ((dY - dj - dl <= vK(0)) & (dl <= vK(1)) & (dj <= vK(2))) {
        if ((dY - dj - dl <= vK(0)) && (dl <= vK(1)) && (dj <= vK(2))) {

          // dFoo = Rf_dbinom(dY - dj - dl, vK(0), vPi(0), 1) +
          //   Rf_dbinom(dl, vK(1), vPi(1), 1) + Rf_dbinom(dj, vK(2), vPi(2), 1);

          dFoo = dSumLogFac -
            vLogFactorial_K(iY - j - l) - vLogFactorial_K(iK1 - iY + j + l) -
            vLogFactorial_K(l) - vLogFactorial_K(iK2 - l) - vLogFactorial_K(j) -
            vLogFactorial_K(iK3 - j) + (dY - dj - dl) * vPi_log(0) +
            (vK(0) - dY + dj + dl) * v1mPi_log(0) + dl * vPi_log(1) + (vK(1) - dl) * v1mPi_log(1) +
            dj * vPi_log(2) + (vK(2) - dj) * v1mPi_log(2);


          vFoo(0) = dLPMF;
          vFoo(1) = dFoo;

          dLPMF = LogSumExp(vFoo);

        }

      }
    }

  }

  if (!bLog) {
    dLPMF = exp(dLPMF);
  }

  return dLPMF;

}

//[[Rcpp::export]]
double dSumThreeBinom_Approx(double dY, arma::vec vK, arma::vec vPi, bool bLog = true) {

  double dMu = accu(vK % vPi);
  double dSD = pow(accu(vK % (1.0 - vPi) % vPi), 0.5);

  double dLPMF = Rf_dnorm4(dY, dMu, dSD, 1);

  if (!bLog) {
    dLPMF = exp(dLPMF);
  }

  return dLPMF;

}

//[[Rcpp::export]]
double dSumMixThreeBinom(double dY, double dK, arma::vec vPi, arma::vec vOmega, arma::vec vLogFactorial_K, bool bLog = true) {

  int i;
  int h;

  double di;
  double dh;

  double dLPMF = 0.0;

  arma::vec vK(3);

  arma::vec vLog_Omega = log(vOmega);
  arma::vec vGamma_Foo(3);

  int iK = (int)dK;

  int iMax = iK + 1;

  double dConst = vLogFactorial_K(iMax - 1);
  double dMu = 0.0;
  double dSigma = 0.0;
  arma::vec vPQ = vPi % (1.0 - vPi);
  double dConstGauss = -0.5 * log(2.0 * M_PI);

  for(i = 0; i<= dK; i++) {
    di = (double)i;
    for(h = 0; h <= dK - di; h++) {
      dh = (double)h;

      vK(0) = dK - di - dh;
      vK(1) = dh;
      vK(2) = di;

      dMu    = vK(0) * vPi(0) + vK(1) * vPi(1) + vK(2) * vPi(2);
      dSigma = pow(vK(0) * vPQ(0) + vK(1) * vPQ(1) + vK(2) * vPQ(2), 0.5);

      dLPMF +=  exp(dConst +
        vK(0) * vLog_Omega(0) - vLogFactorial_K(iK - i - h) +
        vK(1) * vLog_Omega(1) - vLogFactorial_K(h) +
        vK(2) * vLog_Omega(2) - vLogFactorial_K(i) +
        dConstGauss - log(dSigma) - 0.5 * pow((dY - dMu)/dSigma, 2.0));


    }
  }
  if (bLog) {
    dLPMF = log(dLPMF);
  }

  return dLPMF;
}

//[[Rcpp::export]]
double dSumMixThreeBinom_Approx(double dY, double dK, arma::vec vPi, arma::vec vOmega, arma::vec vLogFactorial_K, bool bLog = true) {

  int i;
  int h;

  double di;
  double dh;

  double dLPMF = 0.0;

  arma::vec vK(3);

  arma::vec vLog_Omega = log(vOmega);
  arma::vec vGamma_Foo(3);

  int iK = (int)dK;

  int iMax = iK + 1;

  double dConst = vLogFactorial_K(iMax - 1);
  double dMu = 0.0;
  double dSigma = 0.0;
  arma::vec vPQ = vPi % (1.0 - vPi);
  double dConstGauss = -0.5 * log(2.0 * M_PI);

  double dAlpha = 0.9;

  double dMax_i = (dK + pow(vOmega(1)/vOmega(2), 1.0/dAlpha) +  pow(vOmega(0)/vOmega(2), 1.0/dAlpha) - 2.0)/
    (1.0 + pow(vOmega(1)/vOmega(2), 1.0/dAlpha) +  pow(vOmega(0)/vOmega(2), 1.0/dAlpha));

  double dMax_h =  pow(vOmega(1)/vOmega(2), 1.0/dAlpha) * (dMax_i - 1.0) + 1.0;


  int Max_i = (int)dMax_i;
  int Max_h = (int)dMax_h;

  int iLower_i = Max_i - 20;
  int iLower_h = Max_h - 20;

  int iUpper_i = Max_i + 20;
  int iUpper_h = Max_h + 20;

  if (iLower_i < 0) {
    iLower_i = 0;
  }
  if (iUpper_i > iK) {
    iUpper_i = iK;
  }

  if (iLower_h < 0) {
    iLower_h = 0;
  }

  for(i = iLower_i; i<= iUpper_i; i++) {
    di = (double)i;
    for(h = iLower_h; h <= dK - di; h++) {

      if (h <= iUpper_h) {

        dh = (double)h;

        vK(0) = dK - di - dh;
        vK(1) = dh;
        vK(2) = di;

        dMu    = vK(0) * vPi(0) + vK(1) * vPi(1) + vK(2) * vPi(2);
        dSigma = pow(vK(0) * vPQ(0) + vK(1) * vPQ(1) + vK(2) * vPQ(2), 0.5);

        dLPMF +=  exp(dConst +
          vK(0) * vLog_Omega(0) - vLogFactorial_K(iK - i - h) +
          vK(1) * vLog_Omega(1) - vLogFactorial_K(h) +
          vK(2) * vLog_Omega(2) - vLogFactorial_K(i) +
          dConstGauss - log(dSigma) - 0.5 * pow((dY - dMu)/dSigma, 2.0));

      }
    }
  }
  if (bLog) {
    dLPMF = log(dLPMF);
  }

  return dLPMF;
}

//[[Rcpp::export]]
arma::mat dSumMixThreeBinom_foo(double dY, double dK, arma::vec vPi, arma::vec vOmega, arma::vec vLogFactorial_K, bool bLog = true) {

  int i;
  int h;

  double di;
  double dh;

  double dLPMF = 0.0;

  arma::vec vK(3);

  arma::vec vLog_Omega = log(vOmega);
  arma::vec vGamma_Foo(3);

  int iK = (int)dK;

  int iMax = iK + 1;

  double dConst = vLogFactorial_K(iMax - 1);
  double dMu = 0.0;
  double dSigma = 0.0;
  arma::vec vPQ = vPi % (1.0 - vPi);
  double dConstGauss = -0.5 * log(2.0 * M_PI);

  arma::mat mK_Foo = zeros(iMax, iMax);

  for(i = 0; i<= dK; i++) {
    di = (double)i;
    for(h = 0; h <= dK - di; h++) {
      dh = (double)h;

      vK(0) = dK - di - dh;
      vK(1) = dh;
      vK(2) = di;

      dMu    = vK(0) * vPi(0) + vK(1) * vPi(1) + vK(2) * vPi(2);
      dSigma = pow(vK(0) * vPQ(0) + vK(1) * vPQ(1) + vK(2) * vPQ(2), 0.5);

      mK_Foo(i, h) = exp(dConst +
        vK(0) * vLog_Omega(0) - vLogFactorial_K(iK - i - h) +
        vK(1) * vLog_Omega(1) - vLogFactorial_K(h) +
        vK(2) * vLog_Omega(2) - vLogFactorial_K(i) +
        dConstGauss - log(dSigma) - 0.5 * pow((dY - dMu)/dSigma, 2.0));

      dLPMF += mK_Foo(i, h) ;

    }
  }
  if (bLog) {
    dLPMF = log(dLPMF);
  }

  return mK_Foo;
}

//[[Rcpp::export]]
double dPoissonSumMixThreeBinom(double dY, double dLambda, arma::vec vPi, arma::vec vOmega, bool bLog = true) {

  int iY = (int)dY;

  int iMax = iY + 100;

  int k;

  double dLPMF = -INFINITY;

  double dFoo = 0.0;
  arma::vec vFoo(2);
  arma::vec vK(3);

  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K.zeros();

  vLogFactorial_K(0) = 0;

  double dLogLambda = log(dLambda);

  for (k = 0; k < iMax; k++) {

    if (k > 0) {
      vLogFactorial_K(k) = log(k) + vLogFactorial_K(k - 1);
    }

    if (k >= iY) {

      dFoo = -dLambda + (double)k * dLogLambda - vLogFactorial_K(k) +
        dSumMixThreeBinom(dY, (double)k, vPi, vOmega, vLogFactorial_K, true);

      vFoo(0) = dLPMF;
      vFoo(1) = dFoo;

      dLPMF = LogSumExp(vFoo);

    }

  }

  if (!bLog) {
    dLPMF = exp(dLPMF);
  }

  return dLPMF;

}

////[[Rcpp::plugins(openmp)]]
//[[Rcpp::export]]
double dLLKPoissonSumMixThreeBinom(arma::vec vY, double dLambda, arma::vec vPi,
                                   arma::vec vOmega, arma::vec vLogFactorial_K, int iT, int iMax, int iCores = 7) {

  // omp_set_num_threads(iCores);

  int k;
  int t;

  double dLLK = 0.0;

  double dLPMF;

  arma::vec vK(3);

  double dLogLambda = log(dLambda);
  double dY;

  int iLower = (int) Rf_qpois(0.001, dLambda, 1, 0) ;
  int iUpper = (int) Rf_qpois(0.001, dLambda, 0, 0) ;

  int iUpper_foo = iUpper;

  for (t = 0; t < iT; t++) {
    // dLPMF = -INFINITY;
    dLPMF = 0.0;

    dY = vY(t);

    iUpper_foo = iUpper;

    if (iUpper_foo < dY) {
      iUpper_foo = (int)dY + 20;
    }
    if (iUpper_foo > iMax) {
      iUpper_foo = iMax;
    }

    // #pragma omp parallel for reduction(+: dLPMF)
    for (k = iLower; k < iUpper_foo; k++) {
      if (k >= dY) {
        dLPMF += exp(-dLambda + (double)k * dLogLambda - vLogFactorial_K(k) +
          dSumMixThreeBinom(dY, (double)k, vPi, vOmega, vLogFactorial_K, true));

      }
    }

    dLLK += log(dLPMF);

  }

  return dLLK;

}

////[[Rcpp::plugins(openmp)]]
//[[Rcpp::export]]
double dLLKPoissonSumMixThreeBinom_fast(arma::vec vY, double dLambda, arma::vec vPi,
                                        arma::vec vOmega, arma::vec vLogFactorial_K, int iT, int iMax, int iCores = 7) {

  // omp_set_num_threads(iCores);

  int k;
  int t;

  double dLLK = 0.0;

  double dLPMF;

  arma::vec vK(3);

  double dLogLambda = log(dLambda);
  double dY;

  int iLower = (int) Rf_qpois(0.001, dLambda, 1, 0) ;
  int iUpper = (int) Rf_qpois(0.001, dLambda, 0, 0) ;

  int iUpper_foo = iUpper;

  for (t = 0; t < iT; t++) {
    // dLPMF = -INFINITY;
    dLPMF = 0.0;

    dY = vY(t);

    iUpper_foo = iUpper;

    if (iUpper_foo < dY) {
      iUpper_foo = (int)dY + 20;
    }
    if (iUpper_foo > iMax) {
      iUpper_foo = iMax;
    }

    // #pragma omp parallel for reduction(+: dLPMF)
    for (k = iLower; k < iUpper_foo; k++) {
      if (k >= dY) {
        dLPMF += exp(-dLambda + (double)k * dLogLambda - vLogFactorial_K(k) +
          dSumMixThreeBinom_Approx(dY, (double)k, vPi, vOmega, vLogFactorial_K, true));

      }
    }

    dLLK += log(dLPMF);

  }

  return dLLK;

}

//[[Rcpp::export]]
double dLLKSumMixThreeBinom(arma::vec vY, double dK, arma::vec vPi,
                            arma::vec vOmega, arma::vec vLogFactorial_K, int iT) {

  int t;

  double dLLK = 0.0;

  double dY;

  for (t = 0; t < iT; t++) {

    dY = vY(t);

    if (dK >= dY) {

      dLLK += dSumMixThreeBinom(dY, dK, vPi, vOmega, vLogFactorial_K, true);

    }

  }

  return dLLK;

}

//[[Rcpp::export]]
arma::vec rSumMixThreeBinom(int iT, double dK, arma::vec vPi, arma::vec vOmega) {

  int t;
  int k;

  arma::vec vY(iT);

  vY.zeros();

  int foo;
  arma::vec vBaz(2);

  if(dK > 0) {

    for (t = 0; t < iT; t++) {
      for (k = 1; k <= dK; k++) {

        foo = rando_index(vOmega);
        vY(t) += Rf_rbinom(1, vPi(foo));

      }
    }

  }

  return vY;

}

//[[Rcpp::export]]
arma::vec rPoissonSumMixThreeBinom(int iT, double dLambda, arma::vec vPi, arma::vec vOmega) {

  int t;
  int k;

  arma::vec vY(iT);

  vY.zeros();

  int foo;
  arma::vec vBaz(2);
  double dK = 0.0;

  for (t = 0; t < iT; t++) {
    dK = Rf_rpois(dLambda);

    if(dK > 0) {
      for (k = 1; k <= dK; k++) {

        foo = rando_index(vOmega);
        vY(t) += Rf_rbinom(1, vPi(foo));

      }
    }

  }

  return vY;

}

////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// TWO COMPONENTS MIXTURE ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
double f_i(double di, double dY, double dK, double dOmega, arma::vec vPi, arma::vec vQ, double dConst) {

  double dMu = vPi(0) * di + vPi(1) * (dK - di);
  double dSigma2 = vQ(0) * di + vQ(1) * (dK - di);

  double dFoo = vPi(0) * vQ(0) - vPi(1) * vQ(1);

  double dOut = log((dK - di)/(di + 1)) + log(dOmega/(1.0 - dOmega)) - 0.5 * (dFoo/dSigma2 + (2.0 * (dY - dMu) * dSigma2 + (dY - dMu) * (dY - dMu) * dFoo)/(dSigma2 * dSigma2));

  return dOut;

}

//[[Rcpp::export]]
double find_iMax(double dY, double dK, double dOmega, arma::vec vPi, arma::vec vQ, double dConst, double dLower, double dUpper,
                 int maxiter = 1e4, double eps = 1e-4) {

  double a=dLower;
  double b=dUpper;

  double x=dLower;
  double x1=dUpper;
  int iter = 1;
  double fa,fx;

  //check
  fa = f_i(a, dY, dK, dOmega, vPi, vQ, dConst);
  fx = f_i(x1, dY, dK, dOmega, vPi, vQ, dConst);

  if(fa * fx > 0){
    return 0.0;
  }

  do
  {
    fa = f_i(a, dY, dK, dOmega, vPi, vQ, dConst);
    fx = f_i(x, dY, dK, dOmega, vPi, vQ, dConst);

    if (fa * fx < 0){
      b = x;
    }else{
      a = x;
    }

    x1 = (a + b)/2.0;
    iter++;

    if (abs3(x1 - x) < eps)
    {
      return x1;
    }
    x = x1;
  }while(iter<maxiter);

  // Rprintf("Bisection Warning: Maximum numeber of iteration reached");
  return 0.0;
}

//[[Rcpp::export]]
double dSumMixTwoBinom(double dY, double dK, arma::vec vPi, double dOmega, arma::vec vLogFactorial_K, bool bLog = true) {

  int iK = (int)dK;
  int i;

  double dConst = vLogFactorial_K(iK) - 0.5 * log(2.0 * M_PI);

  double dLLK = 0.0;

  double dMu = 0.0;
  double dSigma2 = 0.0;

  double di;

  arma::vec vQ = vPi % (1.0 - vPi);
  double dOmega_Log   = log(dOmega);
  double d1mOmega_Log = log(1.0 - dOmega);

  for (i = 0; i <= iK; i++) {

    di = (double)i;

    dMu = vPi(0) * di + vPi(1) * (dK - di);
    dSigma2 = vQ(0) * di + vQ(1) * (dK - di);

    dLLK += exp(dConst - vLogFactorial_K(i) - vLogFactorial_K(iK - i) + di * dOmega_Log + (dK - di) * d1mOmega_Log - 0.5 * (log(dSigma2) + (dY - dMu) * (dY - dMu)/dSigma2));

  }

  if(bLog) {
    dLLK = log(dLLK);
  }

  return dLLK;

}

//[[Rcpp::export]]
double dSumMixTwoBinom_Eff(double dY, double dK, arma::vec vPi, double dOmega, arma::vec vLogFactorial_K, bool bLog = true, double dTol_Log = -17) {

  int iK = (int)dK;
  int i = 0;

  double dConst = vLogFactorial_K(iK) - 0.5 * log(2.0 * M_PI);

  double dLLK = 0.0;

  double dMu = 0.0;
  double dSigma2 = 0.0;

  double di;

  arma::vec vQ = vPi % (1.0 - vPi);
  double dOmega_Log   = log(dOmega);
  double d1mOmega_Log = log(1.0 - dOmega);

  bool bRun = true;
  bool bPassThreshold = false;
  double dFoo = 1.0;

  double di_max = find_iMax(dY, dK, dOmega, vPi, vQ, dConst, 0.9, dK + 1);

  int iTol = 50;

  i = (int)di_max - iTol;
  if (i < 0) {
    i = 0;
  }

  int Max_i = (int)di_max + iTol;

  if (Max_i > iK) {
    Max_i = iK;
  }

  while((i <= Max_i) & bRun) {

    di = (double)i;

    dMu = vPi(0) * di + vPi(1) * (dK - di);
    dSigma2 = vQ(0) * di + vQ(1) * (dK - di);

    dFoo = dConst - vLogFactorial_K(i) - vLogFactorial_K(iK - i) + di * dOmega_Log + (dK - di) * d1mOmega_Log - 0.5 * (log(dSigma2) + (dY - dMu) * (dY - dMu)/dSigma2);

    dLLK += exp(dFoo);

    i += 1;

    if ((dFoo > dTol_Log) & !bPassThreshold) {
      bPassThreshold = true;
    }

    if ((dFoo < dTol_Log) & bPassThreshold) {
      bRun = false;
    }

  }

  if(bLog) {
    dLLK = log(dLLK);
  }

  return dLLK;

}

//[[Rcpp::export]]
arma::vec dSumMixTwoBinom_foo(double dY, double dK, arma::vec vPi, double dOmega, arma::vec vLogFactorial_K, bool bLog = true) {

  int iK = (int)dK;
  int i;

  double dConst = vLogFactorial_K(iK) - 0.5 * log(2.0 * M_PI);

  double dMu = 0.0;
  double dSigma2 = 0.0;

  double di;

  arma::vec vQ = vPi % (1.0 - vPi);
  double dOmega_Log   = log(dOmega);
  double d1mOmega_Log = log(1.0 - dOmega);

  arma::vec vFoo(iK + 1);

  for (i = 0; i <= iK; i++) {

    di = (double)i;

    dMu = vPi(0) * di + vPi(1) * (dK - di);
    dSigma2 = vQ(0) * di + vQ(1) * (dK - di);

    vFoo(i) = exp(dConst - vLogFactorial_K(i) - vLogFactorial_K(iK - i) + di * dOmega_Log + (dK - di) * d1mOmega_Log - 0.5 * (log(dSigma2) + (dY - dMu) * (dY - dMu)/dSigma2));

  }

  return vFoo;

}

////[[Rcpp::plugins(openmp)]]
//[[Rcpp::export]]
double dLLKPoissonSumMixTwoBinom(arma::vec vY, double dLambda, arma::vec vPi,
                                 double dOmega, arma::vec vLogFactorial_K, int iT, int iMax, int iCores = 7) {

  // omp_set_num_threads(iCores);

  int k;
  int t;

  double dLLK = 0.0;

  double dLPMF;

  double dLogLambda = log(dLambda);
  double dY;

  int iLower = (int) Rf_qpois(0.001, dLambda, 1, 0) ;
  int iUpper = (int) Rf_qpois(0.001, dLambda, 0, 0) ;

  int iUpper_foo = iUpper;

  for (t = 0; t < iT; t++) {
    // dLPMF = -INFINITY;
    dLPMF = 0.0;

    dY = vY(t);

    iUpper_foo = iUpper;

    if (iUpper_foo < dY) {
      iUpper_foo = (int)dY + 20;
    }
    if (iUpper_foo > iMax) {
      iUpper_foo = iMax;
    }

    // #pragma omp parallel for reduction(+: dLPMF)
    for (k = iLower; k < iUpper_foo; k++) {
      if (k >= dY) {
        dLPMF += exp(-dLambda + (double)k * dLogLambda - vLogFactorial_K(k) +
          dSumMixTwoBinom_Eff(dY, (double)k, vPi, dOmega, vLogFactorial_K, true));

      }
    }

    dLLK += log(dLPMF);

  }

  return dLLK;

}

//[[Rcpp::export]]
arma::vec rSumMixTwoBinom(int iT, double dK, arma::vec vPi, double dOmega) {

  int t;
  int k;

  arma::vec vY(iT);

  vY.zeros();

  int foo;
  arma::vec vBaz(2);

  arma::vec vOmega(2);

  vOmega(0) = dOmega;
  vOmega(1) = 1.0 - dOmega;

  if(dK > 0) {

    for (t = 0; t < iT; t++) {
      for (k = 1; k <= dK; k++) {

        foo = rando_index(vOmega);
        vY(t) += Rf_rbinom(1, vPi(foo));

      }
    }

  }

  return vY;

}

//[[Rcpp::export]]
arma::vec rNBSumMixTwoBinom(int iT, double dA, double dB, arma::vec vPi, double dOmega) {

  int t;
  int k;

  arma::vec vY(iT);

  vY.zeros();

  int foo;
  arma::vec vBaz(2);
  double dK = 0.0;

  arma::vec vOmega(2);

  vOmega(0) = dOmega;
  vOmega(1) = 1.0 - dOmega;

  for (t = 0; t < iT; t++) {
    dK = Rf_rnbinom(dA, dB);

    if(dK > 0) {
      for (k = 1; k <= dK; k++) {

        foo = rando_index(vOmega);
        vY(t) += Rf_rbinom(1, vPi(foo));

      }
    }

  }

  return vY;

}

//[[Rcpp::export]]
arma::vec rPoissonSumMixTwoBinom_Seasonal(int iT, int iS, arma::vec vEta, arma::vec vPi, double dOmega) {

  int t;
  int k;
  int s;

  arma::vec vY(iT);

  vY.zeros();

  int foo;
  arma::vec vBaz(2);
  double dK = 0.0;

  arma::vec vOmega(2);

  vOmega(0) = dOmega;
  vOmega(1) = 1.0 - dOmega;

  arma::vec vLogLambda(iS);
  arma::vec vLambda(iS);

  // compute periodic pattern
  for (s = 0; s < iS; s++) {
    vLogLambda(s) = vEta(0) + vEta(1) * cos(dTwo_Pi * ((s + 1.0)/(iS * 1.0)) - vEta(2) * M_PI);
    vLambda(s) = exp(vLogLambda(s));
  }

  s = 0;

  double dLambda;

  for (t = 0; t < iT; t++) {

    dLambda = vLambda(s);

    dK = Rf_rpois(dLambda);

    if(dK > 0) {
      for (k = 1; k <= dK; k++) {

        foo = rando_index(vOmega);
        vY(t) += Rf_rbinom(1, vPi(foo));

      }
    }

    if (s == iS - 1) {
      s = 0;
    } else {
      s += 1;
    }

  }

  return vY;

}

////[[Rcpp::plugins(openmp)]]
//[[Rcpp::export]]
Rcpp::List dLLKPoissonSumMixTwoBinom_Seasonal_pre(arma::vec vY, arma::vec vPi,
                                                  double dOmega, arma::vec vEta, arma::vec vLogFactorial_K,
                                                  int iT, int iS, int iMax) {

  int k;
  int t;
  int s;

  double dLLK = 0.0;

  double dLPMF;

  double dY;

  arma::vec vLogLambda(iS);
  arma::vec vLambda(iS);
  arma::vec vLower(iS);
  arma::vec vUpper(iS);

  int iLower[iS];
  int iUpper[iS];

  double dThreshold = 0.001;

  // compute periodic pattern
  for (s = 0; s < iS; s++) {
    vLogLambda(s) = vEta(0) + vEta(1) * cos(dTwo_Pi * ((s + 1.0)/(iS * 1.0)) - vEta(2) * M_PI);
    vLambda(s) = exp(vLogLambda(s));
    iLower[s] = (int) Rf_qpois(dThreshold, vLambda(s), 1, 0);
    iUpper[s] = (int) Rf_qpois(dThreshold, vLambda(s), 0, 0);
  }

  s = 0;
  int iUpper_foo;
  int iLower_foo;
  double dLogLambda;
  double dLambda;

  for (t = 0; t < iT; t++) {

    // dLPMF = -INFINITY;
    dLPMF = 0.0;

    dY = vY(t);

    iUpper_foo = iUpper[s];

    if (iUpper_foo < dY) {
      iUpper_foo = (int)dY + 20;
    }
    if (iUpper_foo > iMax) {
      iUpper_foo = iMax;
    }

    iLower_foo = iLower[s];
    dLogLambda = vLogLambda(s);
    dLambda    = vLambda(s);

    // #pragma omp parallel for reduction(+: dLPMF)
    for (k = iLower_foo; k < iUpper_foo; k++) {
      if (k >= dY) {
        dLPMF += exp(-dLambda + (double)k * dLogLambda - vLogFactorial_K(k) +
          dSumMixTwoBinom_Eff(dY, (double)k, vPi, dOmega, vLogFactorial_K, true));

      }
    }

    dLLK += log(dLPMF);

    if (s == iS - 1) {
      s = 0;
    } else {
      s += 1;
    }

  }

  List lOut;

  lOut["dLLK"] = dLLK;
  lOut["vLogLambda"] = vLogLambda;

  return lOut;

}

//// [[Rcpp::plugins(openmp)]]
//[[Rcpp::export]]
Rcpp::List dLLKPoissonSumMixTwoBinom_Seasonal(arma::vec vY, arma::vec vPi,
                                              double dOmega, arma::vec vEta, arma::vec vLogFactorial_K,
                                              int iT, int iS, int iMax) {

  int k;
  int t;
  int s;

  double dLLK = 0.0;

  double dY;

  arma::vec vLogLambda(iS);
  arma::vec vLambda(iS);
  arma::vec vLower(iS);
  arma::vec vUpper(iS);

  // int iLower[iS];
  int iUpper[iS];

  double dThreshold = 0.001;

  // compute periodic pattern
  for (s = 0; s < iS; s++) {
    vLogLambda(s) = vEta(0) + vEta(1) * cos(dTwo_Pi * ((s + 1.0)/(iS * 1.0)) - vEta(2) * M_PI);
    vLambda(s) = exp(vLogLambda(s));
    // iLower[s] = (int) Rf_qpois(dThreshold, vLambda(s), 1, 0);
    iUpper[s] = (int) Rf_qpois(dThreshold, vLambda(s), 0, 0);
  }

  s = 0;
  int iUpper_foo;
  int iLower_foo;
  double dLogLambda;
  double dLambda;

  arma::vec vLPMF(iMax);
  arma::vec vLLK(iT);

  for (t = 0; t < iT; t++) {

    dY = vY(t);

    iUpper_foo = iUpper[s];

    if (iUpper_foo < dY) {
      iUpper_foo = (int)dY + 20;
    }

    if (iUpper_foo > iMax) {
      iUpper_foo = iMax;
    }

    dLogLambda = vLogLambda(s);
    dLambda    = vLambda(s);

    if (iLower_foo < dY) {
      iLower_foo = (int)dY;
    }

    // #pragma omp parallel for
    for (k = iLower_foo; k < iUpper_foo; k++) {
      vLPMF(k) = -dLambda + (double)k * dLogLambda - vLogFactorial_K(k) +
        dSumMixTwoBinom_Eff(dY, (double)k, vPi, dOmega, vLogFactorial_K, true);

    }

    vLLK(t) = LogSumExp(vLPMF.subvec(iLower_foo, iUpper_foo - 1));

    dLLK += vLLK(t);

    if (s == iS - 1) {
      s = 0;
    } else {
      s += 1;
    }

  }

  List lOut;

  lOut["dLLK"] = dLLK;
  lOut["vLLK"] = vLLK;
  lOut["vLogLambda"] = vLogLambda;


  return lOut;

}



////[[Rcpp::plugins(openmp)]]
//[[Rcpp::export]]
double dLLKNBSumMixTwoBinom(arma::vec vY, double dOmega, arma::vec vPi,
                            double dA, double dB, arma::vec vLogFactorial_K, int iT, int iMax) {

  int k;
  int t;

  double dLLK = 0.0;

  double dLPMF;

  double dY;

  int iLower = (int) Rf_qnbinom(0.001, dA, dB, 1, 0) ;
  int iUpper = (int) Rf_qnbinom(0.001, dA, dB, 0, 0) ;

  int iUpper_foo = iUpper;
  if (iLower == 0) {
    iLower = 1;
  }

  for (t = 0; t < iT; t++) {
    // dLPMF = -INFINITY;
    dLPMF = 0.0;

    dY = vY(t);

    iUpper_foo = iUpper;

    if (iUpper_foo < dY) {
      iUpper_foo = (int)dY + 20;
    }
    if (iUpper_foo > iMax) {
      iUpper_foo = iMax;
    }

    // #pragma omp parallel for reduction(+: dLPMF)
    for (k = iLower; k < iUpper_foo; k++) {
      if (k >= dY) {

        dLPMF += exp(Rf_dnbinom((double)k, dA, dB, 1) +
          dSumMixTwoBinom_Eff(dY, (double)k, vPi, dOmega, vLogFactorial_K, true));;

      }
    }

    if(dLPMF < 1e-10) {
      dLPMF = 1e-10;
    }

    dLLK += log(dLPMF);

  }

  return dLLK;

}


