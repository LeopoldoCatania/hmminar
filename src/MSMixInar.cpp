#include <RcppArmadillo.h>
#include "Reparameterization.h"
#include "FFBS.h"
#include "PoissonBinomial.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List EM_MSMixInar(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::vec vLambda,
                  arma::vec vAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                  arma::mat mD, arma::vec vLogFactorial_K, double dTol = 1e-4, int iMaxiter = 1000) {



  int iQ = iJ * iM;
  int iT = vY.size();
  double dY0 = vY(0);
  vY = vY.subvec(1, iT - 1);
  iT = iT - 1;

  int iD = mD.n_cols;

  int j;
  int j1;
  int j2;
  int l;
  int m;
  int m1;
  int m2;
  int t;
  int d;
  int q;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep
  arma::mat mu_all(iQ, iT);

  arma::mat mGamma_Next = mGamma;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vAlpha_Next = vAlpha;
  arma::vec vDelta_Next = vDelta;
  arma::vec vBeta_Next = vBeta;

  arma::vec vAlpha_Num_Foo(iJ);
  arma::vec vAlpha_Den_Foo(iJ);

  arma::vec vLambda_Num_Foo(iM);
  arma::vec vLambda_Den_Foo(iM);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  //initializations
  // double av [iJ][iJ][iM][iM][iT];
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iM];
      for (m1 = 0; m1 < iM; m1++) {
        av[j1][j2][m1] = new double*[iM];
        for (m2 = 0; m2 < iM; m2++) {
          av[j1][j2][m1][m2] = new double[iT];
        }
      }
    }
  }

  // double au [iJ][iM][iT];
  double*** au = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      au[j1][m1] = new double[iT];
    }
  }

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aLLK[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aLLK[j1][m1] = new double[iT];
    }
  }

  // double aeps [iJ][iM][iT];
  double*** aeps = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aeps[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aeps[j1][m1] = new double[iT];
    }
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization
  arma::vec LLKSeries(iMaxiter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iter = 0;
  double eps = 1.0;

  arma::mat mGamma_L = zeros(iQ, iQ);
  arma::vec vDelta_L(iQ);

  while (eps > dTol && iter < iMaxiter) {

    mGamma_L = Gamma_enlarged(mGamma, vOmega);
    vDelta_L = vDelta_enlarged(vDelta, vOmega);

    // compute densities
    for (t = 0; t < iT; t++) {
      for (m = 0; m < iM; m++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[j][m][t] = dPB(vY(t), dY0, vAlpha(j), vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[j][m][t] = dPB(vY(t), vY(t - 1), vAlpha(j), vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }
        }
      }
      // compute densities for the filter
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
      }
    }

    //// FFBS
    // fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);
    fb = FFBS(exp(mLLK.t()), vDelta_L, mGamma_L, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogSumExp(mlalpha.col(iT - 1));

    ////////////////// E STEP ////////////////////////////////
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);
        mu_all(q, t) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
      }
      // populate aeps;
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          if (t == 0) {
            aeps[j][m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) +
              dPB(vY(t) - 1, dY0, vAlpha(j), vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[j][m][t]);
          } else {
            aeps[j][m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) +
              dPB(vY(t) - 1, vY(t - 1), vAlpha(j), vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[j][m][t]);
          }
        }
      }
      if (t > 0) {
        // populate av
        for (j = 0; j < iQ; j++) {
          for (l = 0; l < iQ; l++) {
            av[(int)mIndices(j, 0)]
            [(int)mIndices(l, 0)]
            [(int)mIndices(j, 1)]
            [(int)mIndices(l, 1)][t] = mGamma_L(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
          }
        }
      }
    }


    ////////////////// M STEP ////////////////////////////////

    //ALPHA
    vAlpha_Num_Foo.zeros();
    vAlpha_Den_Foo.zeros();
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        for (t = 0; t < iT; t++) {
          vAlpha_Num_Foo(j) += (au[j][m][t] * (vY(t) - aeps[j][m][t]));
          if (t == 0) {
            vAlpha_Den_Foo(j) += (au[j][m][t] * dY0);
          } else {
            vAlpha_Den_Foo(j) += (au[j][m][t] * vY(t - 1));
          }
        }
      }
      vAlpha_Next(j) = vAlpha_Num_Foo(j)/vAlpha_Den_Foo(j);
    }

    //DELTA
    vDelta_Next.zeros();
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        vDelta_Next(j) += au[j][m][0];
      }
    }

    vDelta_Next = vDelta_Next/accu(vDelta_Next);
    //
    //OMEGA
    vOmega_Next.zeros();

    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        vOmega_Next(m) += au[j][m][0];
      }
    }

    for (m1 = 0; m1 < iM; m1++) {
      for (t = 0; t < iT; t++) {
        for (j1 = 0; j1 < iJ; j1++) {
          for (j2 = 0; j2 < iJ; j2++) {
            for (m2 = 0; m2 < iM; m2++) {
              vOmega_Next(m1) += av[j1][j2][m1][m2][t];
            }
          }
        }
      }
    }

    vOmega_Next = vOmega_Next/accu(vOmega_Next);

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (m = 0; m < iM; m++) {
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          vLambda_Num_Foo(m) += au[j][m][t] * aeps[j][m][t];
          vLambda_Den_Foo(m) += au[j][m][t] * vBeta(vSeason(t));
        }
      }
      vLambda_Next(m) = vLambda_Num_Foo(m)/vLambda_Den_Foo(m);
    }

    //BETA
    vBeta_Num_Foo.zeros();
    vBeta_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vBeta_Num_Foo(vSeason(t)) += au[j][m][t] * aeps[j][m][t];
          vBeta_Den_Foo(vSeason(t)) += au[j][m][t] * vLambda(m);
        }
      }
    }

    vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    //GAMMA
    mGamma_Next.zeros();
    for (t = 1; t < iT; t++) {
      for (j1 = 0; j1 < iJ; j1++) {
        for (j2 = 0; j2 < iJ; j2++) {
          for (m1 = 0; m1 < iM; m1++) {
            for (m2 = 0; m2 < iM; m2++) {
              mGamma_Next(j1, j2) += av[j2][j1][m1][m2][t];
            }
          }
        }
      }
    }

    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/accu(mGamma_Next.row(j));
    }

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    vAlpha = vAlpha_Next;
    vDelta = vDelta_Next;
    vOmega = vOmega_Next;
    vLambda = vLambda_Next;
    mGamma = mGamma_Next;
    if(iD > 1) {
    vBeta = vBeta_Next;
    }

  }

  delete[] av;
  delete[] au;
  delete[] aLLK;
  delete[] aeps;

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  List lOut;

  lOut["vAlpha"] = vAlpha;
  lOut["vDelta"] = vDelta;
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["vBeta"] = vBeta;
  lOut["mGamma"] = mGamma;
  lOut["mGamma_L"] = mGamma_L;
  lOut["LLKSeries"] = LLKSeries;
  lOut["iJ"] = iJ;
  lOut["iM"] = iM;
  lOut["mLLK"] = mLLK;
  lOut["mu_all"] = mu_all;
  lOut["dLLK"] = dLLK;

  return lOut;

}

//[[Rcpp::export]]
List EM_MSMixInar_OmegaChain(arma::vec vY, int iJ, int iM, arma::mat mOmega, arma::vec vLambda,
                  double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                  arma::mat mD, arma::vec vLogFactorial_K, double dTol = 1e-4, int iMaxiter = 1000) {

  int iQ = iJ * iM;
  int iT = vY.size();
  double dY0 = vY(0);
  vY = vY.subvec(1, iT - 1);
  iT = iT - 1;

  int iD = mD.n_cols;

  int j;
  int j1;
  int j2;
  int l;
  int m;
  int m1;
  int m2;
  int t;
  int d;
  int q;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep
  arma::mat mu_all(iQ, iT);

  arma::mat mGamma_Next = mGamma;
  arma::mat mOmega_Next = mOmega;
  arma::vec vLambda_Next = vLambda;
  double dAlpha_Next = dAlpha;
  arma::vec vDelta_Next = vDelta;
  arma::vec vBeta_Next = vBeta;

  double dAlpha_Num_Foo = 0.0;
  double dAlpha_Den_Foo = 0.0;

  arma::vec vLambda_Num_Foo(iM);
  arma::vec vLambda_Den_Foo(iM);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  //initializations
  // double av [iJ][iJ][iM][iM][iT];
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iM];
      for (m1 = 0; m1 < iM; m1++) {
        av[j1][j2][m1] = new double*[iM];
        for (m2 = 0; m2 < iM; m2++) {
          av[j1][j2][m1][m2] = new double[iT];
        }
      }
    }
  }

  // double au [iJ][iM][iT];
  double*** au = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      au[j1][m1] = new double[iT];
    }
  }

  // double au [iJ][iM][iT];
  double** aLLK = new double*[iM];
  for (m1 = 0; m1 < iM; m1++) {
    aLLK[m1] = new double[iT];
  }

  // double aeps [iJ][iM][iT];
  double** aeps = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aeps[m1] = new double[iT];
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization
  arma::vec LLKSeries(iMaxiter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iter = 0;
  double eps = 1.0;

  arma::mat mGamma_L = Gamma_Enlarged_2_OmegaChain(mGamma, mOmega, mIndices);
  arma::vec vDelta_L = getDelta(mGamma_L, iQ);

  while (eps > dTol && iter < iMaxiter) {

    mGamma_L = Gamma_Enlarged_2_OmegaChain(mGamma, mOmega, mIndices);

    // compute densities
    for (t = 0; t < iT; t++) {
      for (m = 0; m < iM; m++) {
          if (t == 0) {
            aLLK[m][t] = dPB(vY(t), dY0, dAlpha, vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[m][t] = dPB(vY(t), vY(t - 1), dAlpha, vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }
      }
      // compute densities for the filter
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 1)][t];
      }
    }

    //// FFBS
    fb = FFBS(exp(mLLK.t()), vDelta_L, mGamma_L, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogSumExp(mlalpha.col(iT - 1));

    ////////////////// E STEP ////////////////////////////////
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);
        mu_all(q, t) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
      }
      // populate aeps;
        for (m = 0; m < iM; m++) {
          if (t == 0) {
            aeps[m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) +
              dPB(vY(t) - 1, dY0, dAlpha, vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[m][t]);
          } else {
            aeps[m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) +
              dPB(vY(t) - 1, vY(t - 1), dAlpha, vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[m][t]);
          }
        }
      if (t > 0) {
        // populate av
        for (j = 0; j < iQ; j++) {
          for (l = 0; l < iQ; l++) {
            av[(int)mIndices(j, 0)]
            [(int)mIndices(l, 0)]
            [(int)mIndices(j, 1)]
            [(int)mIndices(l, 1)][t] = mGamma_L(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
          }
        }
      }
    }

    ////////////////// M STEP ////////////////////////////////

    //ALPHA
    dAlpha_Num_Foo = 0.0;
    dAlpha_Den_Foo = 0.0;
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        for (t = 0; t < iT; t++) {
          dAlpha_Num_Foo += (au[j][m][t] * (vY(t) - aeps[m][t]));
          if (t == 0) {
            dAlpha_Den_Foo += (au[j][m][t] * dY0);
          } else {
            dAlpha_Den_Foo += (au[j][m][t] * vY(t - 1));
          }
        }
      }
    }
    dAlpha_Next = dAlpha_Num_Foo/dAlpha_Den_Foo;

    //DELTA
    vDelta_L.zeros();
    for (q = 0; q < iQ; q++) {
      vDelta_L(q) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][0];
    }

    vDelta_L = vDelta_L/accu(vDelta_L);
    //
    //OMEGA
    mOmega_Next.zeros();

    for (t = 1; t < iT; t++) {
      for (m = 0; m < iM; m++) {
        for (j = 0; j < iJ; j++) {
          mOmega_Next(m, j) += au[j][m][t];
        }
      }
    }
    for (j = 0; j < iJ; j++) {
      mOmega_Next.col(j) = mOmega_Next.col(j)/accu(mOmega_Next.col(j));
    }

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (m = 0; m < iM; m++) {
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          vLambda_Num_Foo(m) += au[j][m][t] * aeps[m][t];
          vLambda_Den_Foo(m) += au[j][m][t] * vBeta(vSeason(t));
        }
      }
      vLambda_Next(m) = vLambda_Num_Foo(m)/vLambda_Den_Foo(m);
    }

    //BETA
    vBeta_Num_Foo.zeros();
    vBeta_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vBeta_Num_Foo(vSeason(t)) += au[j][m][t] * aeps[m][t];
          vBeta_Den_Foo(vSeason(t)) += au[j][m][t] * vLambda(m);
        }
      }
    }

    vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    //GAMMA
    mGamma_Next.zeros();
    for (t = 1; t < iT; t++) {
      for (j1 = 0; j1 < iJ; j1++) {
        for (j2 = 0; j2 < iJ; j2++) {
          for (m1 = 0; m1 < iM; m1++) {
            for (m2 = 0; m2 < iM; m2++) {
              // mGamma_Next(j1, j2) += av[j2][j1][m1][m2][t];
              mGamma_Next(j1, j2) += av[j2][j1][m2][m1][t];
            }
          }
        }
      }
    }

    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/accu(mGamma_Next.row(j));
    }

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    dAlpha = dAlpha_Next;
    mOmega = mOmega_Next;
    vLambda = vLambda_Next;
    mGamma = mGamma_Next;
    if(iD > 1) {
      vBeta = vBeta_Next;
    }

  }

  delete[] av;
  delete[] au;
  delete[] aLLK;
  delete[] aeps;

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  List lOut;

  lOut["dAlpha"] = dAlpha;
  lOut["mOmega"] = mOmega;
  lOut["vLambda"] = vLambda;
  lOut["vBeta"] = vBeta;
  lOut["mGamma"] = mGamma;
  lOut["mGamma_L"] = mGamma_L;
  lOut["LLKSeries"] = LLKSeries;
  lOut["iJ"] = iJ;
  lOut["iM"] = iM;
  lOut["mLLK"] = mLLK;
  lOut["mu_all"] = mu_all;
  lOut["dLLK"] = dLLK;

  return lOut;

}


//[[Rcpp::export]]
List EM_MSMixInar_v2(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::vec vLambda, arma::vec vMu,
                  double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                  arma::mat mD, double dTol = 1e-4, int iMaxiter = 1000) {



  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  int j;
  int j1;
  int j2;
  int l;
  int m;
  int m1;
  int m2;
  int t;
  int d;
  int q;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  arma::mat mGamma_Next = mGamma;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vMu_Next = vMu;
  arma::vec vDelta_Next = vDelta;
  arma::vec vBeta_Next = vBeta;
  double dAlpha_Next = dAlpha;

  double dAlpha_Num_Foo;
  double dAlpha_Den_Foo;

  arma::vec vLambda_Num_Foo(iM);
  arma::vec vLambda_Den_Foo(iM);

  arma::vec vMu_Num_Foo(iJ);
  arma::vec vMu_Den_Foo(iJ);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  //initializations
  // double av [iJ][iJ][iM][iM][iT];
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iM];
      for (m1 = 0; m1 < iM; m1++) {
        av[j1][j2][m1] = new double*[iM];
        for (m2 = 0; m2 < iM; m2++) {
          av[j1][j2][m1][m2] = new double[iT];
        }
      }
    }
  }

  // double au [iJ][iM][iT];
  double*** au = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      au[j1][m1] = new double[iT];
    }
  }

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aLLK[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aLLK[j1][m1] = new double[iT];
    }
  }

  // double aeps [iJ][iM][iT];
  double*** aeps = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aeps[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aeps[j1][m1] = new double[iT];
    }
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization
  arma::vec LLKSeries(iMaxiter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  int iter = 0;
  double eps = 1.0;

  arma::mat mGamma_L = zeros(iQ, iQ);
  arma::vec vDelta_L(iQ);

  while (eps > dTol && iter < iMaxiter) {

    mGamma_L = Gamma_enlarged(mGamma, vOmega);
    vDelta_L = vDelta_enlarged(vDelta, vOmega);

    // compute densities
    for (t = 0; t < iT; t++) {
      for (m = 0; m < iM; m++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[j][m][t] = dPB(vY(t), vY(t), dAlpha, vMu(j) * vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[j][m][t] = dPB(vY(t), vY(t - 1), dAlpha, vMu(j) * vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }
        }
      }
      // compute densities for the filter
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
      }
    }

    //// FFBS
    fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogSumExp(mlalpha.col(iT - 1));

    ////////////////// E STEP ////////////////////////////////
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);
      }
      // populate aeps;
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          if (t == 0) {
            aeps[j][m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) + log(vMu(j)) +
              dPB(vY(t) - 1, vY(t), dAlpha, vMu(j) * vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[j][m][t]);
          } else {
            aeps[j][m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) + log(vMu(j)) +
              dPB(vY(t) - 1, vY(t - 1), dAlpha, vMu(j) * vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[j][m][t]);
          }
        }
      }
      if (t > 0) {
        // populate av
        for (j = 0; j < iQ; j++) {
          for (l = 0; l < iQ; l++) {
            av[(int)mIndices(j, 0)]
            [(int)mIndices(l, 0)]
            [(int)mIndices(j, 1)]
            [(int)mIndices(l, 1)][t] = mGamma_L(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
          }
        }
      }
    }

    ////////////////// M STEP ////////////////////////////////

    //ALPHA
    dAlpha_Num_Foo = 0.0;
    dAlpha_Den_Foo = 0.0;
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        for (t = 1; t < iT; t++) {
          dAlpha_Num_Foo += (au[j][m][t] * (vY(t) - aeps[j][m][t]));
          dAlpha_Den_Foo += (au[j][m][t] * vY(t - 1));
        }
      }
    }
    dAlpha_Next = dAlpha_Num_Foo/dAlpha_Den_Foo;

    //DELTA
    vDelta_Next.zeros();
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        vDelta_Next(j) += au[j][m][0];
      }
    }

    vDelta_Next = vDelta_Next/accu(vDelta_Next);

    //OMEGA
    vOmega_Next.zeros();

    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        vOmega_Next(m) += au[j][m][0];
      }
    }

    for (m1 = 0; m1 < iM; m1++) {
      for (t = 1; t < iT; t++) {
        for (j1 = 0; j1 < iJ; j1++) {
          for (j2 = 0; j2 < iJ; j2++) {
            for (m2 = 0; m2 < iM; m2++) {
              vOmega_Next(m1) += av[j1][j2][m1][m2][t];
            }
          }
        }
      }
    }

    vOmega_Next = vOmega_Next/accu(vOmega_Next);

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (m = 0; m < iM; m++) {
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          vLambda_Num_Foo(m) += (au[j][m][t] * aeps[j][m][t]);
          vLambda_Den_Foo(m) += (au[j][m][t] * vBeta(vSeason(t)) * vMu(j));
        }
      }
      vLambda_Next(m) = vLambda_Num_Foo(m)/vLambda_Den_Foo(m);
    }

    //BETA
    vBeta_Num_Foo.zeros();
    vBeta_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vBeta_Num_Foo(vSeason(t)) += (au[j][m][t] * aeps[j][m][t]);
          vBeta_Den_Foo(vSeason(t)) += (au[j][m][t] * vLambda(m) * vMu(j));
        }
      }
    }

    vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    //GAMMA
    mGamma_Next.zeros();
    for (t = 1; t < iT; t++) {
      for (j1 = 0; j1 < iJ; j1++) {
        for (j2 = 0; j2 < iJ; j2++) {
          for (m1 = 0; m1 < iM; m1++) {
            for (m2 = 0; m2 < iM; m2++) {
              mGamma_Next(j1, j2) += av[j2][j1][m1][m2][t];
            }
          }
        }
      }
    }

    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/accu(mGamma_Next.row(j));
    }


    //Mu
    vMu_Num_Foo.zeros();
    vMu_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vMu_Num_Foo(j) += (au[j][m][t] * aeps[j][m][t]);
          vMu_Den_Foo(j) += (au[j][m][t] * vLambda(m) * vBeta(vSeason(t)));
        }
      }
    }

    vMu_Next = vMu_Num_Foo/vMu_Den_Foo;
    vMu_Next = vMu_Next/accu(vMu_Next);

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    dAlpha = dAlpha_Next;
    vMu    = vMu_Next;
    vDelta = vDelta_Next;
    vOmega = vOmega_Next;
    vLambda = vLambda_Next;
    mGamma = mGamma_Next;
    vBeta = vBeta_Next;

  }

  delete[] av;
  delete[] au;
  delete[] aLLK;
  delete[] aeps;

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  List lOut;

  lOut["dAlpha"] = dAlpha;
  lOut["vMu"] = vMu;
  lOut["vDelta"] = vDelta;
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["vBeta"] = vBeta;
  lOut["mGamma"] = mGamma;
  lOut["LLKSeries"] = LLKSeries;

  lOut["iJ"] = iJ;
  lOut["iM"] = iM;

  return lOut;

}

//[[Rcpp::export]]
List Filtering_MSMixInar_v2_Seasonal_Dummy_Inv(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::vec vLambda, arma::vec vMu,
                                               double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                                               arma::mat mD) {

  int t;
  int m;
  int j;
  int d;
  int q;

  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j = 0; j < iJ; j++) {
    aLLK[j] = new double*[iM];
    for (m = 0; m < iM; m++) {
      aLLK[j][m] = new double[iT];
    }
  }

  arma::mat mGamma_L = Gamma_enlarged(mGamma, vOmega);
  arma::vec vDelta_L = vDelta_enlarged(vDelta, vOmega);

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization

  // compute densities
  for (t = 0; t < iT; t++) {
    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          aLLK[j][m][t] = dPB(vY(t), vY(t), dAlpha, vLambda(m) * vBeta(vSeason(t))/vMu(j), vLogFactorial_K, true);
        } else {
          aLLK[j][m][t] = dPB(vY(t), vY(t - 1), dAlpha, vLambda(m) * vBeta(vSeason(t))/vMu(j), vLogFactorial_K, true);
        }
      }
    }
    // compute densities for the filter
    for (q = 0; q < iQ; q++) {
      mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
    }
  }

  //// FFBS
  fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

  mlalpha = AccessListMat(fb, "mlalpha");
  mlbeta  = AccessListMat(fb, "mlbeta");

  double dLLK = LogSumExp(mlalpha.col(iT - 1));

  arma::mat PredictedProb_L(iT + 1, iQ);
  arma::mat FilteredProb_L(iT, iQ);
  arma::mat SmoothedProb_L(iT, iQ);

  arma::mat PredictedProb = zeros(iT + 1, iJ);
  arma::mat FilteredProb = zeros(iT, iJ);
  arma::mat SmoothedProb = zeros(iT, iJ);

  arma::vec vMean_Eps = zeros(iT + 1);
  arma::vec vE2_Eps = zeros(iT + 1);

  arma::vec vVar = zeros(iT + 1);
  arma::vec vMean = zeros(iT + 1);

  double dC;
  double dLLK_foo;

  PredictedProb_L.row(0) = vDelta_L.t();
  PredictedProb.row(0) = vDelta.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb_L.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb_L.row(t + 1) = FilteredProb_L.row(t) * mGamma_L;
    SmoothedProb_L.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);

    for (q = 0; q < iQ; q++) {
      PredictedProb(t + 1, (int)mIndices(q, 0)) += PredictedProb_L(t + 1, q);
      FilteredProb(t, (int)mIndices(q, 0)) += PredictedProb_L(t, q);
      SmoothedProb(t, (int)mIndices(q, 0)) += SmoothedProb_L(t, q);
      if (t + 1 < iT) {
        vMean_Eps(t + 1) += (PredictedProb_L(t + 1, q) * (vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))/vMu((int)mIndices(q, 0))));
        vE2_Eps(t + 1) += (PredictedProb_L(t + 1, q) * ((vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))/vMu((int)mIndices(q, 0))) *
          ((vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))/vMu((int)mIndices(q, 0))) + 1.0)));
      }
    }
    if (t + 1 < iT) {
      vMean(t + 1) = dAlpha * vY(t) + vMean_Eps(t + 1);
      vVar(t + 1) = dAlpha * (1.0 - dAlpha) * vY(t) + vE2_Eps(t + 1) - pow(vMean_Eps(t + 1), 2.0);
    }
  }

  List lOut;

  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;
  lOut["vMu"] = vMean;
  lOut["vSigma2"] = vVar;

  return lOut;

}


// mLambda is iJ x iM
//[[Rcpp::export]]
List EM_MSMixInar_v3(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::mat mLambda,
                  double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                  arma::mat mD, double dTol = 1e-4, int iMaxiter = 1000) {



  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  int j;
  int j1;
  int j2;
  int l;
  int m;
  int m1;
  int m2;
  int t;
  int d;
  int q;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  arma::mat mGamma_Next = mGamma;
  arma::vec vOmega_Next = vOmega;
  arma::mat mLambda_Next = mLambda;
  double dAlpha_Next = dAlpha;
  arma::vec vDelta_Next = vDelta;
  arma::vec vBeta_Next = vBeta;

  double dAlpha_Num_Foo = 0.0;
  double dAlpha_Den_Foo = 0.0;

  double dLambda_Num_Foo = 0.0;
  double dLambda_Den_Foo = 0.0;

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  //initializations
  // double av [iJ][iJ][iM][iM][iT];
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iM];
      for (m1 = 0; m1 < iM; m1++) {
        av[j1][j2][m1] = new double*[iM];
        for (m2 = 0; m2 < iM; m2++) {
          av[j1][j2][m1][m2] = new double[iT];
        }
      }
    }
  }

  // double au [iJ][iM][iT];
  double*** au = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      au[j1][m1] = new double[iT];
    }
  }

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aLLK[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aLLK[j1][m1] = new double[iT];
    }
  }

  // double aeps [iJ][iM][iT];
  double*** aeps = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aeps[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aeps[j1][m1] = new double[iT];
    }
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization
  arma::vec LLKSeries(iMaxiter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  int iter = 0;
  double eps = 1.0;

  arma::mat mGamma_L = zeros(iQ, iQ);
  arma::vec vDelta_L(iQ);

  while (eps > dTol && iter < iMaxiter) {

    mGamma_L = Gamma_enlarged(mGamma, vOmega);
    vDelta_L = vDelta_enlarged(vDelta, vOmega);

    // compute densities
    for (t = 0; t < iT; t++) {
      for (m = 0; m < iM; m++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[j][m][t] = dPB(vY(t), vY(t), dAlpha, mLambda(j, m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[j][m][t] = dPB(vY(t), vY(t - 1), dAlpha, mLambda(j, m) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }
        }
      }
      // compute densities for the filter
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
      }
    }

    //// FFBS
    fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogSumExp(mlalpha.col(iT - 1));

    ////////////////// E STEP ////////////////////////////////
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);
      }
      // populate aeps;
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          if (t == 0) {
            aeps[j][m][t] = exp(log(mLambda(j, m)) + log(vBeta(vSeason(t))) +
              dPB(vY(t) - 1, vY(t), dAlpha, mLambda(j, m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[j][m][t]);
          } else {
            aeps[j][m][t] = exp(log(mLambda(j, m)) + log(vBeta(vSeason(t))) +
              dPB(vY(t) - 1, vY(t - 1), dAlpha, mLambda(j, m) * vBeta(vSeason(t)), vLogFactorial_K, true) - aLLK[j][m][t]);
          }
        }
      }
      if (t > 0) {
        // populate av
        for (j = 0; j < iQ; j++) {
          for (l = 0; l < iQ; l++) {
            av[(int)mIndices(j, 0)]
            [(int)mIndices(l, 0)]
            [(int)mIndices(j, 1)]
            [(int)mIndices(l, 1)][t] = mGamma_L(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
          }
        }
      }
    }


    ////////////////// M STEP ////////////////////////////////

    //ALPHA
    dAlpha_Num_Foo = 0.0;
    dAlpha_Den_Foo = 0.0;
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        for (t = 1; t < iT; t++) {
          dAlpha_Num_Foo += au[j][m][t] * (vY(t) - aeps[j][m][t]);
          dAlpha_Den_Foo += au[j][m][t] * vY(t - 1);
        }
      }
    }

    dAlpha_Next = dAlpha_Num_Foo/dAlpha_Den_Foo;

    //DELTA
    vDelta_Next.zeros();
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        vDelta_Next(j) += au[j][m][0];
      }
    }

    vDelta_Next = vDelta_Next/accu(vDelta_Next);

    //OMEGA
    vOmega_Next.zeros();

    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        vOmega_Next(m) += au[j][m][0];
      }
    }

    for (m1 = 0; m1 < iM; m1++) {
      for (t = 1; t < iT; t++) {
        for (j1 = 0; j1 < iJ; j1++) {
          for (j2 = 0; j2 < iJ; j2++) {
            for (m2 = 0; m2 < iM; m2++) {
              vOmega_Next(m1) += av[j1][j2][m1][m2][t];
            }
          }
        }
      }
    }

    vOmega_Next = vOmega_Next/accu(vOmega_Next);

    //LAMBDA
    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        dLambda_Num_Foo = 0.0;
        dLambda_Den_Foo = 0.0;
        for (t = 0; t < iT; t++) {
          dLambda_Num_Foo += au[j][m][t] * aeps[j][m][t];
          dLambda_Den_Foo += au[j][m][t] * vBeta(vSeason(t));
        }
       mLambda_Next(j, m) = dLambda_Num_Foo/dLambda_Den_Foo;
      }
    }

    //BETA
    vBeta_Num_Foo.zeros();
    vBeta_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vBeta_Num_Foo(vSeason(t)) += au[j][m][t] * aeps[j][m][t];
          vBeta_Den_Foo(vSeason(t)) += au[j][m][t] * mLambda(j, m);
        }
      }
    }

    vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    //GAMMA
    mGamma_Next.zeros();
    for (t = 1; t < iT; t++) {
      for (j1 = 0; j1 < iJ; j1++) {
        for (j2 = 0; j2 < iJ; j2++) {
          for (m1 = 0; m1 < iM; m1++) {
            for (m2 = 0; m2 < iM; m2++) {
              mGamma_Next(j1, j2) += av[j2][j1][m1][m2][t];
            }
          }
        }
      }
    }

    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/accu(mGamma_Next.row(j));
    }

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    dAlpha = dAlpha_Next;
    vDelta = vDelta_Next;
    vOmega = vOmega_Next;
    mLambda = mLambda_Next;
    mGamma = mGamma_Next;
    vBeta = vBeta_Next;

  }

  delete[] av;
  delete[] au;
  delete[] aLLK;
  delete[] aeps;

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  List lOut;

  lOut["vAlpha"] = dAlpha;
  lOut["vDelta"] = vDelta;
  lOut["vOmega"] = vOmega;
  lOut["mLambda"] = mLambda;
  lOut["vBeta"] = vBeta;
  lOut["mGamma"] = mGamma;
  lOut["LLKSeries"] = LLKSeries;
  lOut["iJ"] = iJ;
  lOut["iM"] = iM;

  return lOut;

}



//[[Rcpp::export]]
List Filtering_MSMixInar_Seasonal_Dummy(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::vec vLambda,
                                           arma::vec vAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                                           arma::mat mD) {

  int t;
  int m;
  int j;
  int d;
  int q;

  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j = 0; j < iJ; j++) {
    aLLK[j] = new double*[iM];
    for (m = 0; m < iM; m++) {
      aLLK[j][m] = new double[iT];
    }
  }

  arma::mat mGamma_L = Gamma_enlarged(mGamma, vOmega);
  arma::vec vDelta_L = vDelta_enlarged(vDelta, vOmega);

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization

  // compute densities
  for (t = 0; t < iT; t++) {
    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          aLLK[j][m][t] = dPB(vY(t), vY(t), vAlpha(j), vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
        } else {
          aLLK[j][m][t] = dPB(vY(t), vY(t - 1), vAlpha(j),  vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
        }
      }
    }
    // compute densities for the filter
    for (q = 0; q < iQ; q++) {
      mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
    }
  }

  //// FFBS
  fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

  mlalpha = AccessListMat(fb, "mlalpha");
  mlbeta  = AccessListMat(fb, "mlbeta");

  double dLLK = LogSumExp(mlalpha.col(iT - 1));

  arma::mat PredictedProb_L(iT + 1, iQ);
  arma::mat FilteredProb_L(iT, iQ);
  arma::mat SmoothedProb_L(iT, iQ);

  arma::mat PredictedProb = zeros(iT + 1, iJ);
  arma::mat FilteredProb = zeros(iT, iJ);
  arma::mat SmoothedProb = zeros(iT, iJ);

  arma::vec vVar = zeros(iT + 1);
  arma::vec vMean = zeros(iT + 1);
  arma::vec vE2 = zeros(iT + 1);

  double dC;
  double dLLK_foo;

  PredictedProb_L.row(0) = vDelta_L.t();
  PredictedProb.row(0) = vDelta.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb_L.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb_L.row(t + 1) = FilteredProb_L.row(t) * mGamma_L;
    SmoothedProb_L.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);

    for (q = 0; q < iQ; q++) {
      PredictedProb(t + 1, (int)mIndices(q, 0)) += PredictedProb_L(t + 1, q);
      FilteredProb(t, (int)mIndices(q, 0)) += PredictedProb_L(t, q);
      SmoothedProb(t, (int)mIndices(q, 0)) += SmoothedProb_L(t, q);
      if (t + 1 < iT) {
        vMean(t + 1) += (PredictedProb_L(t + 1, q) * (vAlpha((int)mIndices(q, 0)) * vY(t) + vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))));
        vE2(t + 1) += (PredictedProb_L(t + 1, q) * (
          //2nd moment B
          vAlpha((int)mIndices(q, 0)) * (1.0 - vAlpha((int)mIndices(q, 0)) * vY(t)) + pow(vAlpha((int)mIndices(q, 0)) * vY(t), 2.0) +
          //2nd moment eps
          (vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) * ((vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) + 1.0) +
          // 2 times prod of first moments
          2.0 * (vAlpha((int)mIndices(q, 0)) * vY(t) * vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) ));
      }
    }
    if (t + 1 < iT) {
      vVar(t + 1) = vE2(t + 1) - pow(vMean(t + 1), 2.0);
    }
  }

  List lOut;

  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;
  lOut["vMu"] = vMean;
  lOut["vSigma2"] = vVar;

  return lOut;

}

//[[Rcpp::export]]
List Filtering_MSMixInar_v3_Seasonal_Dummy(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::mat mLambda,
                                          double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                                          arma::mat mD) {

  int t;
  int m;
  int j;
  int d;
  int q;

  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j = 0; j < iJ; j++) {
    aLLK[j] = new double*[iM];
    for (m = 0; m < iM; m++) {
      aLLK[j][m] = new double[iT];
    }
  }

  arma::mat mGamma_L = Gamma_enlarged(mGamma, vOmega);
  arma::vec vDelta_L = vDelta_enlarged(vDelta, vOmega);

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization

  // compute densities
  for (t = 0; t < iT; t++) {
    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          aLLK[j][m][t] = dPB(vY(t), vY(t), dAlpha, mLambda(j, m) * vBeta(vSeason(t)), vLogFactorial_K, true);
        } else {
          aLLK[j][m][t] = dPB(vY(t), vY(t - 1), dAlpha,  mLambda(j, m) * vBeta(vSeason(t)), vLogFactorial_K, true);
        }
      }
    }
    // compute densities for the filter
    for (q = 0; q < iQ; q++) {
      mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
    }
  }

  //// FFBS
  fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

  mlalpha = AccessListMat(fb, "mlalpha");
  mlbeta  = AccessListMat(fb, "mlbeta");

  double dLLK = LogSumExp(mlalpha.col(iT - 1));

  arma::mat PredictedProb_L(iT + 1, iQ);
  arma::mat FilteredProb_L(iT, iQ);
  arma::mat SmoothedProb_L(iT, iQ);

  arma::mat PredictedProb = zeros(iT + 1, iJ);
  arma::mat FilteredProb = zeros(iT, iJ);
  arma::mat SmoothedProb = zeros(iT, iJ);

  arma::vec vVar = zeros(iT + 1);
  arma::vec vMean = zeros(iT + 1);
  arma::vec vE2 = zeros(iT + 1);

  double dC;
  double dLLK_foo;

  PredictedProb_L.row(0) = vDelta_L.t();
  PredictedProb.row(0) = vDelta.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb_L.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb_L.row(t + 1) = FilteredProb_L.row(t) * mGamma_L;
    SmoothedProb_L.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);

    for (q = 0; q < iQ; q++) {
      PredictedProb(t + 1, (int)mIndices(q, 0)) += PredictedProb_L(t + 1, q);
      FilteredProb(t, (int)mIndices(q, 0)) += PredictedProb_L(t, q);
      SmoothedProb(t, (int)mIndices(q, 0)) += SmoothedProb_L(t, q);
      if (t + 1 < iT) {
        vMean(t + 1) += (PredictedProb_L(t + 1, q) * (dAlpha * vY(t) + mLambda((int)mIndices(q, 0), (int)mIndices(q, 1)) * vBeta(vSeason(t + 1))));
        vE2(t + 1) += (PredictedProb_L(t + 1, q) * (
          //2nd moment B
          dAlpha * (1.0 - dAlpha * vY(t)) + pow(dAlpha * vY(t), 2.0) +
            //2nd moment eps
            (mLambda((int)mIndices(q, 0), (int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) * ((mLambda((int)mIndices(q, 0), (int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) + 1.0) +
            // 2 times prod of first moments
            2.0 * (dAlpha * vY(t) * mLambda((int)mIndices(q, 0), (int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) ));
      }
    }
    if (t + 1 < iT) {
      vVar(t + 1) = vE2(t + 1) - pow(vMean(t + 1), 2.0);
    }
  }

  List lOut;

  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;
  lOut["vMu"] = vMean;
  lOut["vSigma2"] = vVar;

  return lOut;

}


//[[Rcpp::export]]
arma::vec Sim_MSMixInar_Seasonal_Dummy(int iT, arma::vec vOmega, arma::vec vLambda, arma::vec vAlpha, arma::mat mGamma,
                                       arma::vec vDelta, arma::vec vBeta, arma::mat mD) {

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

  arma::vec vState(iT);

  vState(0) = rando_index(vDelta);

  iFoo = rando_index(vOmega);
  vY(0) = Rf_rpois(vLambda(iFoo)/(1.0 - vAlpha(vState(0))) * accu(vBeta) / (iD * 1.0));

  for (t = 1; t < iT; t++){

    vState(t) = rando_index(mGamma.row(vState(t - 1)).t());

    iFoo = rando_index(vOmega);
    vY(t) = Rf_rbinom(vY(t - 1), vAlpha(vState(t))) +  Rf_rpois(vLambda(iFoo) * vBeta(vSeason(t)));

  }

  return vY;

}

//[[Rcpp::export]]
arma::vec Sim_MSMixInar_v2_Seasonal_Dummy(int iT, arma::vec vOmega, arma::vec vLambda, arma::vec vMu, double dAlpha, arma::mat mGamma,
                                       arma::vec vDelta, arma::vec vBeta, arma::mat mD) {

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

  arma::vec vState(iT);

  vState(0) = rando_index(vDelta);

  iFoo = rando_index(vOmega);
  vY(0) = Rf_rpois(vLambda(iFoo)/(1.0 - dAlpha) * accu(vBeta) / (iD * 1.0) * accu(vDelta % vMu));

  for (t = 1; t < iT; t++){

    vState(t) = rando_index(mGamma.row(vState(t - 1)).t());

    iFoo = rando_index(vOmega);
    vY(t) = Rf_rbinom(vY(t - 1), dAlpha) +  Rf_rpois(vMu(vState(t)) * vLambda(iFoo) * vBeta(vSeason(t)));

  }

  return vY;

}

//[[Rcpp::export]]
arma::vec Sim_MSMixInar_v3_Seasonal_Dummy(int iT, arma::vec vOmega, arma::mat mLambda, double dAlpha, arma::mat mGamma,
                                       arma::vec vDelta, arma::vec vBeta, arma::mat mD) {

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

  arma::vec vState(iT);

  vState(0) = rando_index(vDelta);
  iFoo = rando_index(vOmega);

  vY(0) = Rf_rpois(mLambda(vState(0), iFoo)/(1.0 - dAlpha) * accu(vBeta) / (iD * 1.0));

  for (t = 1; t < iT; t++){

    vState(t) = rando_index(mGamma.row(vState(t - 1)).t());

    iFoo = rando_index(vOmega);
    vY(t) = Rf_rbinom(vY(t - 1), dAlpha) +  Rf_rpois(mLambda(vState(t), iFoo) * vBeta(vSeason(t)));

  }

  return vY;

}

//[[Rcpp::export]]
List Filtering_MSMixInar_v2_Seasonal_Dummy(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::vec vLambda, arma::vec vMu,
                                           double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                                           arma::mat mD) {

  int t;
  int m;
  int j;
  int d;
  int q;

  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j = 0; j < iJ; j++) {
    aLLK[j] = new double*[iM];
    for (m = 0; m < iM; m++) {
      aLLK[j][m] = new double[iT];
    }
  }

  arma::mat mGamma_L = Gamma_enlarged(mGamma, vOmega);
  arma::vec vDelta_L = vDelta_enlarged(vDelta, vOmega);

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization

  // compute densities
  for (t = 0; t < iT; t++) {
    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        if (t == 0) {
          aLLK[j][m][t] = dPB(vY(t), vY(t), dAlpha, vMu(j) * vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
        } else {
          aLLK[j][m][t] = dPB(vY(t), vY(t - 1), dAlpha, vMu(j) * vLambda(m) * vBeta(vSeason(t)), vLogFactorial_K, true);
        }
      }
    }
    // compute densities for the filter
    for (q = 0; q < iQ; q++) {
      mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
    }
  }

  //// FFBS
  fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

  mlalpha = AccessListMat(fb, "mlalpha");
  mlbeta  = AccessListMat(fb, "mlbeta");

  double dLLK = LogSumExp(mlalpha.col(iT - 1));

  arma::mat PredictedProb_L(iT + 1, iQ);
  arma::mat FilteredProb_L(iT, iQ);
  arma::mat SmoothedProb_L(iT, iQ);

  arma::mat PredictedProb = zeros(iT + 1, iJ);
  arma::mat FilteredProb = zeros(iT, iJ);
  arma::mat SmoothedProb = zeros(iT, iJ);

  arma::vec vMean_Eps = zeros(iT + 1);
  arma::vec vE2_Eps = zeros(iT + 1);

  arma::vec vVar = zeros(iT + 1);
  arma::vec vMean = zeros(iT + 1);

  double dC;
  double dLLK_foo;

  PredictedProb_L.row(0) = vDelta_L.t();
  PredictedProb.row(0) = vDelta.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb_L.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb_L.row(t + 1) = FilteredProb_L.row(t) * mGamma_L;
    SmoothedProb_L.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);

    for (q = 0; q < iQ; q++) {
      PredictedProb(t + 1, (int)mIndices(q, 0)) += PredictedProb_L(t + 1, q);
      FilteredProb(t, (int)mIndices(q, 0)) += PredictedProb_L(t, q);
      SmoothedProb(t, (int)mIndices(q, 0)) += SmoothedProb_L(t, q);
      if (t + 1 < iT) {
        vMean_Eps(t + 1) += (PredictedProb_L(t + 1, q) * (vMu((int)mIndices(q, 0)) * vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))));
        vE2_Eps(t + 1) += (PredictedProb_L(t + 1, q) * ((vMu((int)mIndices(q, 0)) * vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) *
          ((vMu((int)mIndices(q, 0)) * vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t + 1))) + 1.0)));
      }
    }
    if (t + 1 < iT) {
      vMean(t + 1) = dAlpha * vY(t) + vMean_Eps(t + 1);
      vVar(t + 1) = dAlpha * (1.0 - dAlpha) * vY(t) + vE2_Eps(t + 1) - pow(vMean_Eps(t + 1), 2.0);
    }
  }

  List lOut;

  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;
  lOut["vMu"] = vMean;
  lOut["vSigma2"] = vVar;

  return lOut;

}

//[[Rcpp::export]]
List EM_MSMixInar_v2_Inv(arma::vec vY, int iJ, int iM, arma::vec vOmega, arma::vec vLambda, arma::vec vMu,
                     double dAlpha, arma::mat mGamma, arma::vec vDelta, arma::vec vBeta,
                     arma::mat mD, double dTol = 1e-4, int iMaxiter = 1000) {



  int iQ = iJ * iM;
  int iT = vY.size();
  int iD = mD.n_cols;

  int j;
  int j1;
  int j2;
  int l;
  int m;
  int m1;
  int m2;
  int t;
  int d;
  int q;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  second rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) second rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) second rep

  arma::mat mGamma_Next = mGamma;
  arma::vec vOmega_Next = vOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vMu_Next = vMu;
  arma::vec vDelta_Next = vDelta;
  arma::vec vBeta_Next = vBeta;
  double dAlpha_Next = dAlpha;

  double dAlpha_Num_Foo;
  double dAlpha_Den_Foo;

  arma::vec vLambda_Num_Foo(iM);
  arma::vec vLambda_Den_Foo(iM);

  arma::vec vMu_Num_Foo(iJ);
  arma::vec vMu_Den_Foo(iJ);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  //initializations
  // double av [iJ][iJ][iM][iM][iT];
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iM];
      for (m1 = 0; m1 < iM; m1++) {
        av[j1][j2][m1] = new double*[iM];
        for (m2 = 0; m2 < iM; m2++) {
          av[j1][j2][m1][m2] = new double[iT];
        }
      }
    }
  }

  // double au [iJ][iM][iT];
  double*** au = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      au[j1][m1] = new double[iT];
    }
  }

  // double au [iJ][iM][iT];
  double*** aLLK = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aLLK[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aLLK[j1][m1] = new double[iT];
    }
  }

  // double aeps [iJ][iM][iT];
  double*** aeps = new double**[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    aeps[j1] = new double*[iM];
    for (m1 = 0; m1 < iM; m1++) {
      aeps[j1][m1] = new double[iT];
    }
  }

  arma::mat mIndices = Indices(iJ, iM); //indices in the L reparameterization
  arma::vec LLKSeries(iMaxiter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  // Seasonality
  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  int iMax = (int) max(vY);
  arma::vec vLogFactorial_K(iMax + 1);

  vLogFactorial_K(0) = 0;

  for (t = 1; t < iMax + 1; t++) {
    vLogFactorial_K(t) = vLogFactorial_K(t - 1) + log(t * 1.0);
  }

  int iter = 0;
  double eps = 1.0;

  arma::mat mGamma_L = zeros(iQ, iQ);
  arma::vec vDelta_L(iQ);

  while (eps > dTol && iter < iMaxiter) {

    mGamma_L = Gamma_enlarged(mGamma, vOmega);
    vDelta_L = vDelta_enlarged(vDelta, vOmega);

    // compute densities
    for (t = 0; t < iT; t++) {
      for (m = 0; m < iM; m++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[j][m][t] = dPB(vY(t), vY(t), dAlpha, vLambda(m) * vBeta(vSeason(t)) / vMu(j), vLogFactorial_K, true);
          } else {
            aLLK[j][m][t] = dPB(vY(t), vY(t - 1), dAlpha,  vLambda(m) * vBeta(vSeason(t)) / vMu(j), vLogFactorial_K, true);
          }
        }
      }
      // compute densities for the filter
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 0)][(int)mIndices(q, 1)][t];
      }
    }

    //// FFBS
    fb = FFBS_robust(mLLK.t(), vDelta_L, mGamma_L, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogSumExp(mlalpha.col(iT - 1));

    ////////////////// E STEP ////////////////////////////////
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);
      }
      // populate aeps;
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          if (t == 0) {
            aeps[j][m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) - log(vMu(j)) +
              dPB(vY(t) - 1, vY(t), dAlpha, vLambda(m) * vBeta(vSeason(t)) / vMu(j), vLogFactorial_K, true) - aLLK[j][m][t]);
          } else {
            aeps[j][m][t] = exp(log(vLambda(m)) + log(vBeta(vSeason(t))) - log(vMu(j)) +
              dPB(vY(t) - 1, vY(t - 1), dAlpha, vLambda(m) * vBeta(vSeason(t)) / vMu(j), vLogFactorial_K, true) - aLLK[j][m][t]);
          }
        }
      }
      if (t > 0) {
        // populate av
        for (j = 0; j < iQ; j++) {
          for (l = 0; l < iQ; l++) {
            av[(int)mIndices(j, 0)]
            [(int)mIndices(l, 0)]
            [(int)mIndices(j, 1)]
            [(int)mIndices(l, 1)][t] = mGamma_L(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
          }
        }
      }
    }

    ////////////////// M STEP ////////////////////////////////

    //ALPHA
    dAlpha_Num_Foo = 0.0;
    dAlpha_Den_Foo = 0.0;
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        for (t = 1; t < iT; t++) {
          dAlpha_Num_Foo += (au[j][m][t] * (vY(t) - aeps[j][m][t]));
          dAlpha_Den_Foo += (au[j][m][t] * vY(t - 1));
        }
      }
    }
    dAlpha_Next = dAlpha_Num_Foo/dAlpha_Den_Foo;

    //DELTA
    vDelta_Next.zeros();
    for (j = 0; j < iJ; j++) {
      for (m = 0; m < iM; m++) {
        vDelta_Next(j) += au[j][m][0];
      }
    }

    vDelta_Next = vDelta_Next/accu(vDelta_Next);

    //OMEGA
    vOmega_Next.zeros();

    for (m = 0; m < iM; m++) {
      for (j = 0; j < iJ; j++) {
        vOmega_Next(m) += au[j][m][0];
      }
    }

    for (m1 = 0; m1 < iM; m1++) {
      for (t = 1; t < iT; t++) {
        for (j1 = 0; j1 < iJ; j1++) {
          for (j2 = 0; j2 < iJ; j2++) {
            for (m2 = 0; m2 < iM; m2++) {
              vOmega_Next(m1) += av[j1][j2][m1][m2][t];
            }
          }
        }
      }
    }

    vOmega_Next = vOmega_Next/accu(vOmega_Next);

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (m = 0; m < iM; m++) {
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          vLambda_Num_Foo(m) += (au[j][m][t] * aeps[j][m][t]);
          vLambda_Den_Foo(m) += (au[j][m][t] * vBeta(vSeason(t)) / vMu(j));
        }
      }
      vLambda_Next(m) = vLambda_Num_Foo(m)/vLambda_Den_Foo(m);
    }

    //BETA
    vBeta_Num_Foo.zeros();
    vBeta_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vBeta_Num_Foo(vSeason(t)) += (au[j][m][t] * aeps[j][m][t]);
          vBeta_Den_Foo(vSeason(t)) += (au[j][m][t] * vLambda(m) / vMu(j));
        }
      }
    }

    vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    //GAMMA
    mGamma_Next.zeros();
    for (t = 1; t < iT; t++) {
      for (j1 = 0; j1 < iJ; j1++) {
        for (j2 = 0; j2 < iJ; j2++) {
          for (m1 = 0; m1 < iM; m1++) {
            for (m2 = 0; m2 < iM; m2++) {
              mGamma_Next(j1, j2) += av[j2][j1][m1][m2][t];
            }
          }
        }
      }
    }

    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/accu(mGamma_Next.row(j));
    }


    //Mu
    vMu_Num_Foo.zeros();
    vMu_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (j = 0; j < iJ; j++) {
        for (m = 0; m < iM; m++) {
          vMu_Den_Foo(j) += (au[j][m][t] * aeps[j][m][t]);
          vMu_Num_Foo(j) += (au[j][m][t] * vLambda(m) * vBeta(vSeason(t)));
        }
      }
    }

    vMu_Next = vMu_Num_Foo/vMu_Den_Foo;

    vMu_Next = vMu_Next/accu(vMu_Next);

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    dAlpha = dAlpha_Next;
    vMu    = vMu_Next;
    vDelta = vDelta_Next;
    vOmega = vOmega_Next;
    vLambda = vLambda_Next;
    mGamma = mGamma_Next;
    vBeta = vBeta_Next;

  }

  delete[] av;
  delete[] au;
  delete[] aLLK;
  delete[] aeps;

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  List lOut;

  lOut["dAlpha"] = dAlpha;
  lOut["vMu"] = vMu;
  lOut["vDelta"] = vDelta;
  lOut["vOmega"] = vOmega;
  lOut["vLambda"] = vLambda;
  lOut["vBeta"] = vBeta;
  lOut["mGamma"] = mGamma;
  lOut["LLKSeries"] = LLKSeries;

  lOut["iJ"] = iJ;
  lOut["iM"] = iM;

  return lOut;

}


//[[Rcpp::export]]
arma::vec Sim_MSMixInar_v2_Seasonal_Dummy_Inv(int iT, arma::vec vOmega, arma::vec vLambda, arma::vec vMu, double dAlpha, arma::mat mGamma,
                                          arma::vec vDelta, arma::vec vBeta, arma::mat mD) {

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

  arma::vec vState(iT);

  vState(0) = rando_index(vDelta);

  iFoo = rando_index(vOmega);
  vY(0) = Rf_rpois(vLambda(iFoo)/(1.0 - dAlpha) * (accu(vBeta) / (iD * 1.0)) / accu(vDelta % vMu));

  for (t = 1; t < iT; t++){

    vState(t) = rando_index(mGamma.row(vState(t - 1)).t());

    iFoo = rando_index(vOmega);
    vY(t) = Rf_rbinom(vY(t - 1), dAlpha) +  Rf_rpois(vLambda(iFoo) * vBeta(vSeason(t))/vMu(vState(t)));

  }

  return vY;

}

// iJ is the state space of S^alpha
// iL is the State space of S^eta
// iK is the State space of Z

//[[Rcpp::export]]
List EM_MSMixInarTwoChains_Dummy_C2(arma::vec vY, int iJ, int iK, int iL,
                                   arma::mat mGamma_Eta,
                                   arma::mat mGamma_Alpha,
                                   arma::mat mOmega,
                                   arma::vec vLambda,
                                   arma::vec vAlpha,
                                   arma::vec vBeta,
                                   arma::vec vSeason,
                                   arma::vec vLogFactorial_K,
                                   int maxIter = 1e3, double tol = 1e-4, bool bSeasonal = true
) {

  //////////////////// DEFINITIONS /////////////////////////////////

  int iT = vY.size();
  int iD = vBeta.size();

  double dY0 = vY(0);
  vY = vY.subvec(1, iT - 1);
  vSeason = vSeason.subvec(1, iT - 1);

  iT = iT - 1;

  int iQ = iJ * iL * iK;

  int t;
  int j;
  int k;
  int l;
  int q;

  int j1;
  int j2;
  int l1;
  int l2;
  int k1;
  int k2;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  third rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) third rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) third rep

  arma::cube cv_Alpha = zeros(iJ, iJ, iT); //v_i,j,t P(S_t^alpha, S_{t-1}^alpha | Y_{1:T})
  arma::cube cv_Eta   = zeros(iL, iL, iT); //v_i,j,t P(S_t^eta, S_{t-1}^eta | Y_{1:T})
  // cv_Alpha(j1, j2, t) = P(S_t^alpha = j1, S_{t-1}^alpha = j2 | Y_{1:T})
  // cv_Eta(l1, l2, t) = P(S_t^eta = l1, S_{t-1}^eta = l2 | Y_{1:T})

  arma::mat mu_Eta(iL, iT); // smoothed prob S^Eta
  arma::mat mu_Alpha(iJ, iT); // smoothed prob S^Alpha
  arma::mat mu_Z(iK, iT); // smoothed prob Z
  arma::mat mu_all(iQ, iT); // smoothed prob all
  arma::mat mEta = zeros(iQ, iT);

  arma::cube cZ(iK, iT, iL); // P(Z_t = K, S_t^Eta = j | Y_{1:T})

  //initializations

  // double av [iJ][iJ][iL][iL][iK][iK][iT];
  // av[j1][j2][l1][l2][k1][k2] = P(S_t^alpha = j1, S_{t-1}^alpha = j2,
  //                                S_t^eta = l1, S_{t-1}^eta = l2,
  //                                Z_t = k1, Z_{t-1} = k2| Y_{1:T})
  double******* av = new double******[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double*****[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double****[iL];
      for (l1 = 0; l1 < iL; l1++) {
        av[j1][j2][l1] = new double***[iL];
        for (l2 = 0; l2 < iL; l2++) {
          av[j1][j2][l1][l2] = new double**[iK];
          for (k1 = 0; k1 < iK; k1++) {
            av[j1][j2][l1][l2][k1] = new double*[iK];
            for (k2 = 0; k2 < iK; k2++) {
              av[j1][j2][l1][l2][k1][k2] = new double[iT];
            }
          }
        }
      }
    }
  }

  // double au [iJ][iK][iL][iT];
  // ah[j1][l1][k1] = P(S_t^alpha = j1, S_t^eta = l1, Z_t = k1)
  double**** au = new double***[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double**[iK];
    for (k1 = 0; k1 < iK; k1++) {
      au[j1][k1] = new double*[iL];
      for (l1 = 0; l1 < iL; l1++) {
        au[j1][k1][l1] = new double[iT];
      }
    }
  }

  // double aLLK [iK][iJ][iT];
  // aLLK = P(Y_t \vert S_t^alpha = j1, Z_t = k1, Y_{t-1})
  double*** aLLK = new double**[iK];
  for (k1 = 0; k1 < iK; k1++) {
    aLLK[k1] = new double*[iJ];
    for (j1 = 0; j1 < iJ; j1++) {
      aLLK[k1][j1] = new double[iT];
    }
  }

  //indices in the enlarged reparameterization
  // First index for S_t^alpha
  // Second index for Z_t
  // Third index for S_t^eta
  arma::mat mIndices = OmegaZB_Indices(iJ, iK, iL);

  arma::vec  LLKSeries(maxIter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  int iter = 0;
  double eps = 1.0;

  arma::vec vLambda_Num_Foo(iK);
  arma::vec vLambda_Den_Foo(iK);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  arma::vec vAlpha_Num_Foo(iJ);
  arma::vec vAlpha_Den_Foo(iJ);

  // parameter definition
  arma::mat mGamma_Alpha_Next = mGamma_Alpha;
  arma::mat mGamma_Eta_Next = mGamma_Eta;
  arma::mat mOmega_Next = mOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vBeta_Next = vBeta;
  arma::vec vAlpha_Next = vAlpha;

  arma::mat mGamma_OmegaZB = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);
  arma::vec vDelta_OmegaZB = getDelta(mGamma_OmegaZB, iQ);
  //
  while (eps > tol && iter < maxIter) {

    mGamma_OmegaZB = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);
    /// Compute log densities
    for (t = 0; t < iT; t++) {
      for (k = 0; k < iK; k++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[k][j][t] = dPB(vY(t), dY0, vAlpha(j), vLambda(k) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[k][j][t] = dPB(vY(t), vY(t - 1), vAlpha(j),  vLambda(k) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }
        }
      }
    }
    /// Compute log densities
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 1)][(int)mIndices(q, 0)][t];
      }
    }

    ////////////////// E STEP ////////////////////////////////

    //// FFBS
    // fb = FFBS(exp(mLLK.t()), vDelta_OmegaZB, mGamma_OmegaZB, iQ, iT);
    fb = FFBS_robust(mLLK.t(), vDelta_OmegaZB, mGamma_OmegaZB, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogOfSum(mlalpha.col(iT - 1));

    mu_Alpha.zeros();
    mu_Eta.zeros();
    mu_Z.zeros();
    cZ.zeros();

    //// poulate mlu_OmegaZB
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)]
        [(int)mIndices(q, 2)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);

        mu_all(q, t) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][(int)mIndices(q, 2)][t];

        if (t == 0) {
          mEta(q, t) = exp(log(vLambda((int)mIndices(q, 1))) + log(vBeta(vSeason(t))) +
            dPB(vY(t) - 1.0, dY0, vAlpha((int)mIndices(q, 0)), vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t)), vLogFactorial_K, true) - mLLK(q, t));
        } else {
          mEta(q, t) = exp(log(vLambda((int)mIndices(q, 1))) + log(vBeta(vSeason(t))) +
            dPB(vY(t) - 1.0, vY(t - 1), vAlpha((int)mIndices(q, 0)), vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t)), vLogFactorial_K, true) - mLLK(q, t));
        }

      }
      for (j = 0; j < iJ; j++) {
        for (k = 0; k < iK; k++) {
          for (l = 0; l < iL; l++) {
            mu_Z(k, t)     += au[j][k][l][t];
            mu_Alpha(j, t) += au[j][k][l][t];
            mu_Eta(l, t)   += au[j][k][l][t];
            cZ(k, t, l)    += au[j][k][l][t];
          }
        }
      }
    }

    cv_Alpha.zeros();
    cv_Eta.zeros();
    //
    for (t = 1; t < iT; t++) {
      // populate cv
      for (j = 0; j < iQ; j++) {
        for (l = 0; l < iQ; l++) {

          // av[j1][j2][l1][l2][k1][k2] = P(S_t^alpha = j1, S_{t-1}^alpha = j2,
          //                                S_t^eta = l1, S_{t-1}^eta = l2,
          //                                Z_t = k1, Z_{t-1} = k2| Y_{1:T})

          // mGamma_Enlarged(q1, q2) = mOmega((int)mIndices(q2, 1), (int)mIndices(q2, 2)) * mGamma_Eta((int)mIndices(q1, 2), (int)mIndices(q2, 2)) *
          // mGamma_Alpha((int)mIndices(q1, 0), (int)mIndices(q2, 0));

          av[(int)mIndices(j, 0)]
          [(int)mIndices(l, 0)]
          [(int)mIndices(j, 2)]
          [(int)mIndices(l, 2)]
          [(int)mIndices(j, 1)]
          [(int)mIndices(l, 1)][t] = mGamma_OmegaZB(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
        }
      }
      // populate cv_Omega, vc_B
      for (j1 = 0; j1 < iJ; j1++) {
        for (l1 = 0; l1 < iL; l1++) {
          for (k1 = 0; k1 < iK; k1++) {
            for (j2 = 0; j2 < iJ; j2++) {
              for (l2 = 0; l2 < iL; l2++) {
                for (k2 = 0; k2 < iK; k2++) {

                  // cv_Alpha(j1, j2, t) = P(S_t^alpha = j1, S_{t-1}^alpha = j2 | Y_{1:T})
                  // cv_Eta(l1, l2, t) = P(S_t^eta = l1, S_{t-1}^eta = l2 | Y_{1:T})
                  cv_Alpha(j1, j2, t) += av[j1][j2][l1][l2][k1][k2][t];
                  cv_Eta(l1, l2, t)   += av[j1][j2][l1][l2][k1][k2][t];
                }
              }
            }
          }
        }
      }
    }

    ////////////////// M STEP ////////////////////////////////

    //DELTA
    for (q = 0; q < iQ; q++) {
      vDelta_OmegaZB(q) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][(int)mIndices(q, 2)][0];
    }

    vDelta_OmegaZB = vDelta_OmegaZB/accu(vDelta_OmegaZB);

    // Gamma_omega
    for (j = 0; j < iJ; j++) {
      for (l = 0; l < iJ ; l++) {
        mGamma_Alpha_Next(j, l) = accu(cv_Alpha.tube(l, j));
      }
    }
    // Gamma_B
    for (j = 0; j < iL; j++) {
      for (l = 0; l < iL; l++) {
        mGamma_Eta_Next(j, l) = accu(cv_Eta.tube(l, j));
      }
    }
    // Normalize Gamma_omega
    for (j = 0; j < iJ; j++) {
      mGamma_Alpha_Next.row(j) = mGamma_Alpha_Next.row(j)/sum(mGamma_Alpha_Next.row(j));
    }
    // Normalize Gamma_B
    for (j = 0; j < iL; j++) {
      mGamma_Eta_Next.row(j) = mGamma_Eta_Next.row(j)/sum(mGamma_Eta_Next.row(j));
    }

    // mOmega
    mOmega_Next.zeros();
    // mOmega
    for (t = 1; t < iT; t++) {
      for (k = 0; k < iK; k++) {
        for (l = 0; l < iL; l++) {
          mOmega_Next(k, l) += cZ(k, t, l);
        }
      }
    }

    // normalize omega
    for (l = 0; l < iL; l++) {
      mOmega_Next.col(l) = mOmega_Next.col(l) / sum(mOmega_Next.col(l));
    }

    // mOmega_Next  = OmegaCheck(mOmega_Next, iK, iJ);

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        vLambda_Num_Foo((int)mIndices(q, 1)) +=  (mu_all(q, t) * mEta(q, t));
        vLambda_Den_Foo((int)mIndices(q, 1)) +=  (mu_all(q, t) * vBeta(vSeason(t)));
      }
    }

    vLambda_Next = vLambda_Num_Foo/vLambda_Den_Foo;

    //BETA
    if (bSeasonal) {
      vBeta_Num_Foo.zeros();
      vBeta_Den_Foo.zeros();

      for (t = 0; t < iT; t++) {
        for (q = 0; q < iQ; q++) {
          vBeta_Num_Foo(vSeason(t)) += (mu_all(q, t) * mEta(q, t));
          vBeta_Den_Foo(vSeason(t)) += (mu_all(q, t) * vLambda_Next((int)mIndices(q, 1)));
        }
      }

      vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    }

    //ALPHA
    vAlpha_Num_Foo.zeros();
    vAlpha_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        vAlpha_Num_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * (vY(t) - mEta(q, t)));
        if (t == 0) {
          vAlpha_Den_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * dY0);
        } else {
          vAlpha_Den_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * (vY(t - 1)));
        }
      }
    }

    vAlpha_Next = vAlpha_Num_Foo/vAlpha_Den_Foo;

    for(j = 0; j < iJ; j++) {
      if(vAlpha_Next(j) < 1e-7) {
        vAlpha_Next(j) = 1e-7;
      }
    }

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    // Update Parameters
    mGamma_Eta   = mGamma_Eta_Next;
    mOmega = mOmega_Next;
    mGamma_Alpha = mGamma_Alpha_Next;
    vLambda = vLambda_Next;
    vAlpha = vAlpha_Next;
    if (bSeasonal) {
      vBeta = vBeta_Next;
    }

  }

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  arma::mat PredictedProb(iT + 1, iQ);
  arma::mat FilteredProb(iT, iQ);
  arma::mat SmoothedProb(iT, iQ);

  double dC;
  double dLLK_foo;

  PredictedProb.row(0) = vDelta_OmegaZB.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb.row(t + 1) = FilteredProb.row(t) * mGamma_OmegaZB;
    SmoothedProb.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);
  }

  for (j1 = 0; j1 < iJ; j1++) {
    for (j2 = 0; j2 < iJ; j2++) {
      for (l1 = 0; l1 < iL; l1++) {
        for (l2 = 0; l2 < iL; l2++) {
          for (k1 = 0; k1 < iK; k1++) {
            for (k2 = 0; k2 < iK; k2++) {
              delete av[j1][j2][l1][l2][k1][k2];
            }
          }
        }
      }
    }
  }

  delete[] av;
  for (j1 = 0; j1 < iJ; j1++) {
    for (k1 = 0; k1 < iK; k1++) {
      for (l1 = 0; l1 < iL; l1++) {
        delete au[j1][k1][l1];
      }
    }
  }

  delete[] au;

  for (k1 = 0; k1 < iK; k1++) {
    for (j1 = 0; j1 < iJ; j1++) {
      delete aLLK[k1][j1];
    }
  }

  delete[] aLLK;


  List lOut;

  lOut["mGamma_Eta"] = mGamma_Eta;
  lOut["mGamma_Alpha"] = mGamma_Alpha;
  lOut["mOmega"] = mOmega;
  lOut["vLambda"] = vLambda;
  lOut["vAlpha"] = vAlpha;
  lOut["vBeta"] = vBeta;

  lOut["mLLK"] = mLLK;

  lOut["LLKSeries"] = LLKSeries;
  lOut["eps"] = eps;

  lOut["mu_Alpha"] = mu_Alpha;
  lOut["mu_Z"]     = mu_Z;
  lOut["mu_Eta"]     = mu_Eta;
  lOut["mu_all"]     = mu_all;

  lOut["mIndices"] = mIndices;
  lOut["vSeason"]       = vSeason;

  lOut["mGamma_OmegaZB"] = mGamma_OmegaZB;
  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;

  lOut["dLLK"] = dLLK;

  lOut["mlalpha"] = mlalpha;
  lOut["mlbeta"] = mlbeta;

  lOut["mEta"] = mEta;

  lOut["vY"] = vY;
  lOut["dY0"] = dY0;

  return lOut;

}


// iJ is the state space of S^alpha
// iL is the State space of S^eta
// iK is the State space of Z

//[[Rcpp::export]]
List EM_MSMixInarTwoChains_Dummy_C_AlphaCor(arma::vec vY, int iJ, int iK, int iL,
                                    arma::mat mGamma_Eta,
                                    arma::mat mGamma_Alpha,
                                    arma::mat mOmega,
                                    arma::vec vLambda,
                                    arma::vec vAlpha,
                                    arma::vec vBeta,
                                    double dVarPhi,
                                    arma::vec vSeason,
                                    arma::vec vLogFactorial_K,
                                    arma::vec vO,
                                    int maxIter = 1e3, double tol = 1e-4, bool bSeasonal = true, bool bFilter = false
) {

  //////////////////// DEFINITIONS /////////////////////////////////

  if (bFilter) {
    maxIter = 1;
  }

  int iT = vY.size();
  int iD = vBeta.size();

  double dY0 = vY(0);
  double dO0 = vO(0);
  int iSeason0 = (int)vSeason(0);

  vY = vY.subvec(1, iT - 1);
  vSeason = vSeason.subvec(1, iT - 1);
  vO = vO.subvec(1, iT - 1);

  iT = iT - 1;

  int iQ = iJ * iL * iK;

  int t;
  int j;
  int k;
  int l;
  int q;

  int j1;
  int j2;
  int l1;
  int l2;
  int k1;
  int k2;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  third rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) third rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) third rep

  arma::cube cv_Alpha = zeros(iJ, iJ, iT); //v_i,j,t P(S_t^alpha, S_{t-1}^alpha | Y_{1:T})
  arma::cube cv_Eta   = zeros(iL, iL, iT); //v_i,j,t P(S_t^eta, S_{t-1}^eta | Y_{1:T})
  // cv_Alpha(j1, j2, t) = P(S_t^alpha = j1, S_{t-1}^alpha = j2 | Y_{1:T})
  // cv_Eta(l1, l2, t) = P(S_t^eta = l1, S_{t-1}^eta = l2 | Y_{1:T})

  arma::mat mu_Eta(iL, iT); // smoothed prob S^Eta
  arma::mat mu_Alpha(iJ, iT); // smoothed prob S^Alpha
  arma::mat mu_Z(iK, iT); // smoothed prob Z
  arma::mat mu_all(iQ, iT); // smoothed prob all
  arma::mat mEta = zeros(iQ, iT);

  arma::cube cZ(iK, iT, iL); // P(Z_t = K, S_t^Eta = j | Y_{1:T})

  //initializations

  // double av [iJ][iJ][iL][iL][iK][iK][iT];
  // av[j1][j2][l1][l2][k1][k2] = P(S_t^alpha = j1, S_{t-1}^alpha = j2,
  //                                S_t^eta = l1, S_{t-1}^eta = l2,
  //                                Z_t = k1, Z_{t-1} = k2| Y_{1:T})
  double******* av = new double******[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double*****[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double****[iL];
      for (l1 = 0; l1 < iL; l1++) {
        av[j1][j2][l1] = new double***[iL];
        for (l2 = 0; l2 < iL; l2++) {
          av[j1][j2][l1][l2] = new double**[iK];
          for (k1 = 0; k1 < iK; k1++) {
            av[j1][j2][l1][l2][k1] = new double*[iK];
            for (k2 = 0; k2 < iK; k2++) {
              av[j1][j2][l1][l2][k1][k2] = new double[iT];
            }
          }
        }
      }
    }
  }

  // double au [iJ][iK][iL][iT];
  // ah[j1][l1][k1] = P(S_t^alpha = j1, S_t^eta = l1, Z_t = k1)
  double**** au = new double***[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double**[iK];
    for (k1 = 0; k1 < iK; k1++) {
      au[j1][k1] = new double*[iL];
      for (l1 = 0; l1 < iL; l1++) {
        au[j1][k1][l1] = new double[iT];
      }
    }
  }

  // double aLLK [iK][iJ][iT];
  // aLLK = P(Y_t \vert S_t^alpha = j1, Z_t = k1, Y_{t-1})
  double*** aLLK = new double**[iK];
  for (k1 = 0; k1 < iK; k1++) {
    aLLK[k1] = new double*[iJ];
    for (j1 = 0; j1 < iJ; j1++) {
      aLLK[k1][j1] = new double[iT];
    }
  }

  //indices in the enlarged reparameterization
  // First index for S_t^alpha
  // Second index for Z_t
  // Third index for S_t^eta
  arma::mat mIndices = OmegaZB_Indices(iJ, iK, iL);

  arma::vec  LLKSeries(maxIter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  int iter = 0;
  double eps = 1.0;

  arma::vec vLambda_Num_Foo(iK);
  arma::vec vLambda_Den_Foo(iK);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  arma::vec vAlpha_Num_Foo(iJ);
  arma::vec vAlpha_Den_Foo(iJ);

  double dVarPhi_Num = 0.0;
  double dVarPhi_Den = 0.0;

  // parameter definition
  arma::mat mGamma_Alpha_Next = mGamma_Alpha;
  arma::mat mGamma_Eta_Next = mGamma_Eta;
  arma::mat mOmega_Next = mOmega;
  arma::vec vLambda_Next = vLambda;
  arma::vec vBeta_Next = vBeta;
  arma::vec vAlpha_Next = vAlpha;
  double dVarPhi_Next = dVarPhi;

  arma::mat mGamma_OmegaZB = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);
  arma::vec vDelta_OmegaZB = getDelta(mGamma_OmegaZB, iQ);
  //
  while (eps > tol && iter < maxIter) {

    mGamma_OmegaZB = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);
    /// Compute log densities
    for (t = 0; t < iT; t++) {
      for (k = 0; k < iK; k++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[k][j][t] = dPB(vY(t), dY0, (1 - vO(t)) * vAlpha(j) + vO(t) * dVarPhi, vLambda(k) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[k][j][t] = dPB(vY(t), vY(t - 1),  (1 - vO(t)) * vAlpha(j) + vO(t) * dVarPhi,  vLambda(k) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }

        }
      }
    }
    /// Compute log densities
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 1)][(int)mIndices(q, 0)][t];
      }
    }

    ////////////////// E STEP ////////////////////////////////

    //// FFBS
    // fb = FFBS(exp(mLLK.t()), vDelta_OmegaZB, mGamma_OmegaZB, iQ, iT);
    fb = FFBS_robust(mLLK.t(), vDelta_OmegaZB, mGamma_OmegaZB, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogOfSum(mlalpha.col(iT - 1));

    mu_Alpha.zeros();
    mu_Eta.zeros();
    mu_Z.zeros();
    cZ.zeros();

    //// poulate mlu_OmegaZB
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)]
        [(int)mIndices(q, 2)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);

        mu_all(q, t) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][(int)mIndices(q, 2)][t];

        if (t == 0) {
          mEta(q, t) = exp(log(vLambda((int)mIndices(q, 1))) + log(vBeta(vSeason(t))) +
            dPB(vY(t) - 1.0, dY0, (1 - vO(t)) * vAlpha((int)mIndices(q, 0)) + vO(t) * dVarPhi, vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t)), vLogFactorial_K, true) - mLLK(q, t));
        } else {
          mEta(q, t) = exp(log(vLambda((int)mIndices(q, 1))) + log(vBeta(vSeason(t))) +
            dPB(vY(t) - 1.0, vY(t - 1),  (1 - vO(t)) * vAlpha((int)mIndices(q, 0)) + vO(t) * dVarPhi, vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t)), vLogFactorial_K, true) - mLLK(q, t));
        }

      }
      for (j = 0; j < iJ; j++) {
        for (k = 0; k < iK; k++) {
          for (l = 0; l < iL; l++) {
            mu_Z(k, t)     += au[j][k][l][t];
            mu_Alpha(j, t) += au[j][k][l][t];
            mu_Eta(l, t)   += au[j][k][l][t];
            cZ(k, t, l)    += au[j][k][l][t];
          }
        }
      }
    }

    cv_Alpha.zeros();
    cv_Eta.zeros();
    //
    for (t = 1; t < iT; t++) {
      // populate cv
      for (j = 0; j < iQ; j++) {
        for (l = 0; l < iQ; l++) {

          // av[j1][j2][l1][l2][k1][k2] = P(S_t^alpha = j1, S_{t-1}^alpha = j2,
          //                                S_t^eta = l1, S_{t-1}^eta = l2,
          //                                Z_t = k1, Z_{t-1} = k2| Y_{1:T})

          // mGamma_Enlarged(q1, q2) = mOmega((int)mIndices(q2, 1), (int)mIndices(q2, 2)) * mGamma_Eta((int)mIndices(q1, 2), (int)mIndices(q2, 2)) *
          // mGamma_Alpha((int)mIndices(q1, 0), (int)mIndices(q2, 0));

          av[(int)mIndices(j, 0)]
          [(int)mIndices(l, 0)]
          [(int)mIndices(j, 2)]
          [(int)mIndices(l, 2)]
          [(int)mIndices(j, 1)]
          [(int)mIndices(l, 1)][t] = mGamma_OmegaZB(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
        }
      }
      // populate cv_Omega, vc_B
      for (j1 = 0; j1 < iJ; j1++) {
        for (l1 = 0; l1 < iL; l1++) {
          for (k1 = 0; k1 < iK; k1++) {
            for (j2 = 0; j2 < iJ; j2++) {
              for (l2 = 0; l2 < iL; l2++) {
                for (k2 = 0; k2 < iK; k2++) {

                  // cv_Alpha(j1, j2, t) = P(S_t^alpha = j1, S_{t-1}^alpha = j2 | Y_{1:T})
                  // cv_Eta(l1, l2, t) = P(S_t^eta = l1, S_{t-1}^eta = l2 | Y_{1:T})
                  cv_Alpha(j1, j2, t) += av[j1][j2][l1][l2][k1][k2][t];
                  cv_Eta(l1, l2, t)   += av[j1][j2][l1][l2][k1][k2][t];
                }
              }
            }
          }
        }
      }
    }

    ////////////////// M STEP ////////////////////////////////

    //DELTA
    for (q = 0; q < iQ; q++) {
      vDelta_OmegaZB(q) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][(int)mIndices(q, 2)][0];
    }

    vDelta_OmegaZB = vDelta_OmegaZB/accu(vDelta_OmegaZB);

    // Gamma_omega
    for (j = 0; j < iJ; j++) {
      for (l = 0; l < iJ ; l++) {
        mGamma_Alpha_Next(j, l) = accu(cv_Alpha.tube(l, j));
      }
    }
    // Gamma_B
    for (j = 0; j < iL; j++) {
      for (l = 0; l < iL; l++) {
        mGamma_Eta_Next(j, l) = accu(cv_Eta.tube(l, j));
      }
    }
    // Normalize Gamma_omega
    for (j = 0; j < iJ; j++) {
      mGamma_Alpha_Next.row(j) = mGamma_Alpha_Next.row(j)/sum(mGamma_Alpha_Next.row(j));
    }
    // Normalize Gamma_B
    for (j = 0; j < iL; j++) {
      mGamma_Eta_Next.row(j) = mGamma_Eta_Next.row(j)/sum(mGamma_Eta_Next.row(j));
    }

    // mOmega
    mOmega_Next.zeros();
    // mOmega
    for (t = 1; t < iT; t++) {
      for (k = 0; k < iK; k++) {
        for (l = 0; l < iL; l++) {
          mOmega_Next(k, l) += cZ(k, t, l);
        }
      }
    }

    // normalize omega
    for (l = 0; l < iL; l++) {
      mOmega_Next.col(l) = mOmega_Next.col(l) / sum(mOmega_Next.col(l));
    }

    // mOmega_Next  = OmegaCheck(mOmega_Next, iK, iJ);

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        vLambda_Num_Foo((int)mIndices(q, 1)) +=  (mu_all(q, t) * mEta(q, t));
        vLambda_Den_Foo((int)mIndices(q, 1)) +=  (mu_all(q, t) * vBeta(vSeason(t)));
      }
    }

    vLambda_Next = vLambda_Num_Foo/vLambda_Den_Foo;

    //BETA
    if (bSeasonal) {
      vBeta_Num_Foo.zeros();
      vBeta_Den_Foo.zeros();

      for (t = 0; t < iT; t++) {
        for (q = 0; q < iQ; q++) {
          vBeta_Num_Foo(vSeason(t)) += (mu_all(q, t) * mEta(q, t));
          vBeta_Den_Foo(vSeason(t)) += (mu_all(q, t) * vLambda_Next((int)mIndices(q, 1)));
        }
      }

      vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    }

    //ALPHA
    vAlpha_Num_Foo.zeros();
    vAlpha_Den_Foo.zeros();
    dVarPhi_Num = 0.0;
    dVarPhi_Den = 0.0;

    for (t = 0; t < iT; t++) {
      if (vO(t) > 0)  {
        for (q = 0; q < iQ; q++) {
          dVarPhi_Num += (mu_all(q, t) * (vY(t) - mEta(q, t)));
          if (t == 0) {
            dVarPhi_Den += (mu_all(q, t) * dY0 );
          } else {
            dVarPhi_Den += (mu_all(q, t) * vY(t - 1));
          }
        }
      } else  {
        for (q = 0; q < iQ; q++) {
          vAlpha_Num_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * (vY(t) - mEta(q, t)));
          if (t == 0) {
            vAlpha_Den_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * dY0 );
          } else {
            vAlpha_Den_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * vY(t - 1));
          }
        }
      }

    }

    vAlpha_Next = vAlpha_Num_Foo/vAlpha_Den_Foo;
    dVarPhi_Next = dVarPhi_Num/dVarPhi_Den;

    if (dVarPhi_Next < 1e-7) {
      dVarPhi_Next = dVarPhi;
    }
     for (j = 0; j < iJ; j++) {
       if(vAlpha_Next(j) < 1e-7) {
         vAlpha_Next(j) = vAlpha(j);
       }
     }

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    // Update Parameters
    if (!bFilter) {
      mGamma_Eta   = mGamma_Eta_Next;
      mOmega = mOmega_Next;
      mGamma_Alpha = mGamma_Alpha_Next;
      vLambda = vLambda_Next;
      vAlpha = vAlpha_Next;
      dVarPhi = dVarPhi_Next;
      if (bSeasonal) {
        vBeta = vBeta_Next;
      }
    }

  }

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  arma::mat PredictedProb(iT + 1, iQ);
  arma::mat FilteredProb(iT, iQ);
  arma::mat SmoothedProb(iT, iQ);

  double dC;
  double dLLK_foo;

  PredictedProb.row(0) = vDelta_OmegaZB.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb.row(t + 1) = FilteredProb.row(t) * mGamma_OmegaZB;
    SmoothedProb.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);
  }

  for (j1 = 0; j1 < iJ; j1++) {
    for (j2 = 0; j2 < iJ; j2++) {
      for (l1 = 0; l1 < iL; l1++) {
        for (l2 = 0; l2 < iL; l2++) {
          for (k1 = 0; k1 < iK; k1++) {
            for (k2 = 0; k2 < iK; k2++) {
              delete av[j1][j2][l1][l2][k1][k2];
            }
          }
        }
      }
    }
  }

  delete[] av;
  for (j1 = 0; j1 < iJ; j1++) {
    for (k1 = 0; k1 < iK; k1++) {
      for (l1 = 0; l1 < iL; l1++) {
        delete au[j1][k1][l1];
      }
    }
  }

  delete[] au;

  for (k1 = 0; k1 < iK; k1++) {
    for (j1 = 0; j1 < iJ; j1++) {
      delete aLLK[k1][j1];
    }
  }

  delete[] aLLK;


  List lOut;

  lOut["mGamma_Eta"] = mGamma_Eta;
  lOut["mGamma_Alpha"] = mGamma_Alpha;
  lOut["mOmega"] = mOmega;
  lOut["vLambda"] = vLambda;
  lOut["vAlpha"] = vAlpha;
  lOut["dVarPhi"] = dVarPhi;
  lOut["vBeta"] = vBeta;

  lOut["mLLK"] = mLLK;

  lOut["LLKSeries"] = LLKSeries;
  lOut["eps"] = eps;

  lOut["mu_Alpha"] = mu_Alpha;
  lOut["mu_Z"]     = mu_Z;
  lOut["mu_Eta"]     = mu_Eta;
  lOut["mu_all"]     = mu_all;
  lOut["mEta"] = mEta;

  lOut["mIndices"] = mIndices;
  lOut["vSeason"]       = vSeason;
  lOut["vO"]       = vO;

  lOut["mGamma_OmegaZB"] = mGamma_OmegaZB;
  lOut["vDelta_OmegaZB"] = vDelta_OmegaZB;
  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;

  lOut["dLLK"] = dLLK;

  lOut["mlalpha"] = mlalpha;
  lOut["mlbeta"] = mlbeta;

  lOut["mEta"] = mEta;

  lOut["vY"] = vY;
  lOut["dY0"] = dY0;
  lOut["dO0"] = dO0;
  lOut["iSeason0"] = iSeason0;

  return lOut;

}


/// iL = iK, mOmega = mI

// iJ is the state space of S^alpha
// iL is the State space of S^eta
// iK is the State space of Z

//[[Rcpp::export]]
List EM_MSMixInarTwoChains_Dummy_C_AlphaCor_OmegaIdentity(arma::vec vY, int iJ, int iK, int iL,
                                            arma::mat mGamma_Eta,
                                            arma::mat mGamma_Alpha,
                                            arma::mat mOmega,
                                            arma::vec vLambda,
                                            arma::vec vAlpha,
                                            arma::vec vBeta,
                                            double dVarPhi,
                                            arma::vec vSeason,
                                            arma::vec vLogFactorial_K,
                                            arma::vec vO,
                                            int maxIter = 1e3, double tol = 1e-4, bool bSeasonal = true, bool bFilter = false
) {

  //////////////////// DEFINITIONS /////////////////////////////////

  if (bFilter) {
    maxIter = 1;
  }

  int iT = vY.size();
  int iD = vBeta.size();

  double dY0 = vY(0);
  double dO0 = vO(0);
  int iSeason0 = (int)vSeason(0);

  vY = vY.subvec(1, iT - 1);
  vSeason = vSeason.subvec(1, iT - 1);
  vO = vO.subvec(1, iT - 1);

  iT = iT - 1;

  int iQ = iJ * iL * iK;

  int t;
  int j;
  int k;
  int l;
  int q;

  int j1;
  int j2;
  int l1;
  int l2;
  int k1;
  int k2;

  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob)  third rep
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) third rep
  arma::mat mLLK(iQ, iT);    // log state densities (ZIS) third rep

  arma::cube cv_Alpha = zeros(iJ, iJ, iT); //v_i,j,t P(S_t^alpha, S_{t-1}^alpha | Y_{1:T})
  arma::cube cv_Eta   = zeros(iL, iL, iT); //v_i,j,t P(S_t^eta, S_{t-1}^eta | Y_{1:T})
  // cv_Alpha(j1, j2, t) = P(S_t^alpha = j1, S_{t-1}^alpha = j2 | Y_{1:T})
  // cv_Eta(l1, l2, t) = P(S_t^eta = l1, S_{t-1}^eta = l2 | Y_{1:T})

  arma::mat mu_Eta(iL, iT); // smoothed prob S^Eta
  arma::mat mu_Alpha(iJ, iT); // smoothed prob S^Alpha
  arma::mat mu_Z(iK, iT); // smoothed prob Z
  arma::mat mu_all(iQ, iT); // smoothed prob all
  arma::mat mEta = zeros(iQ, iT);

  arma::cube cZ(iK, iT, iL); // P(Z_t = K, S_t^Eta = j | Y_{1:T})

  //initializations

  // double av [iJ][iJ][iL][iL][iK][iK][iT];
  // av[j1][j2][l1][l2][k1][k2] = P(S_t^alpha = j1, S_{t-1}^alpha = j2,
  //                                S_t^eta = l1, S_{t-1}^eta = l2,
  //                                Z_t = k1, Z_{t-1} = k2| Y_{1:T})
  double******* av = new double******[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double*****[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double****[iL];
      for (l1 = 0; l1 < iL; l1++) {
        av[j1][j2][l1] = new double***[iL];
        for (l2 = 0; l2 < iL; l2++) {
          av[j1][j2][l1][l2] = new double**[iK];
          for (k1 = 0; k1 < iK; k1++) {
            av[j1][j2][l1][l2][k1] = new double*[iK];
            for (k2 = 0; k2 < iK; k2++) {
              av[j1][j2][l1][l2][k1][k2] = new double[iT];
            }
          }
        }
      }
    }
  }

  // double au [iJ][iK][iL][iT];
  // ah[j1][l1][k1] = P(S_t^alpha = j1, S_t^eta = l1, Z_t = k1)
  double**** au = new double***[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    au[j1] = new double**[iK];
    for (k1 = 0; k1 < iK; k1++) {
      au[j1][k1] = new double*[iL];
      for (l1 = 0; l1 < iL; l1++) {
        au[j1][k1][l1] = new double[iT];
      }
    }
  }

  // double aLLK [iK][iJ][iT];
  // aLLK = P(Y_t \vert S_t^alpha = j1, Z_t = k1, Y_{t-1})
  double*** aLLK = new double**[iK];
  for (k1 = 0; k1 < iK; k1++) {
    aLLK[k1] = new double*[iJ];
    for (j1 = 0; j1 < iJ; j1++) {
      aLLK[k1][j1] = new double[iT];
    }
  }

  //indices in the enlarged reparameterization
  // First index for S_t^alpha
  // Second index for Z_t
  // Third index for S_t^eta
  arma::mat mIndices = OmegaZB_Indices(iJ, iK, iL);

  arma::vec  LLKSeries(maxIter + 1); // likelihood at each iteration
  double dLLK = 0.0;

  int iter = 0;
  double eps = 1.0;

  arma::vec vLambda_Num_Foo(iK);
  arma::vec vLambda_Den_Foo(iK);

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  arma::vec vAlpha_Num_Foo(iJ);
  arma::vec vAlpha_Den_Foo(iJ);

  double dVarPhi_Num = 0.0;
  double dVarPhi_Den = 0.0;

  // parameter definition
  arma::mat mGamma_Alpha_Next = mGamma_Alpha;
  arma::mat mGamma_Eta_Next = mGamma_Eta;
  arma::vec vLambda_Next = vLambda;
  arma::vec vBeta_Next = vBeta;
  arma::vec vAlpha_Next = vAlpha;
  double dVarPhi_Next = dVarPhi;

  arma::mat mGamma_OmegaZB = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);
  arma::vec vDelta_OmegaZB = getDelta(mGamma_OmegaZB, iQ);
  //
  while (eps > tol && iter < maxIter) {

    mGamma_OmegaZB = Gamma_Enlarged_2(mGamma_Eta, mOmega, mGamma_Alpha, mIndices);
    /// Compute log densities
    for (t = 0; t < iT; t++) {
      for (k = 0; k < iK; k++) {
        for (j = 0; j < iJ; j++) {
          if (t == 0) {
            aLLK[k][j][t] = dPB(vY(t), dY0, (1 - vO(t)) * vAlpha(j) + vO(t) * dVarPhi, vLambda(k) * vBeta(vSeason(t)), vLogFactorial_K, true);
          } else {
            aLLK[k][j][t] = dPB(vY(t), vY(t - 1),  (1 - vO(t)) * vAlpha(j) + vO(t) * dVarPhi,  vLambda(k) * vBeta(vSeason(t)), vLogFactorial_K, true);
          }

        }
      }
    }
    /// Compute log densities
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        mLLK(q, t) = aLLK[(int)mIndices(q, 1)][(int)mIndices(q, 0)][t];
      }
    }

    ////////////////// E STEP ////////////////////////////////

    //// FFBS
    // fb = FFBS(exp(mLLK.t()), vDelta_OmegaZB, mGamma_OmegaZB, iQ, iT);
    fb = FFBS_robust(mLLK.t(), vDelta_OmegaZB, mGamma_OmegaZB, iQ, iT);

    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");

    dLLK = LogOfSum(mlalpha.col(iT - 1));

    mu_Alpha.zeros();
    mu_Eta.zeros();
    mu_Z.zeros();
    cZ.zeros();

    //// poulate mlu_OmegaZB
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        au[(int)mIndices(q, 0)]
        [(int)mIndices(q, 1)]
        [(int)mIndices(q, 2)][t] = exp(mlalpha(q, t) + mlbeta(q, t) - dLLK);

        mu_all(q, t) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][(int)mIndices(q, 2)][t];

        if (t == 0) {
          mEta(q, t) = exp(log(vLambda((int)mIndices(q, 1))) + log(vBeta(vSeason(t))) +
            dPB(vY(t) - 1.0, dY0, (1 - vO(t)) * vAlpha((int)mIndices(q, 0)) + vO(t) * dVarPhi, vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t)), vLogFactorial_K, true) - mLLK(q, t));
        } else {
          mEta(q, t) = exp(log(vLambda((int)mIndices(q, 1))) + log(vBeta(vSeason(t))) +
            dPB(vY(t) - 1.0, vY(t - 1),  (1 - vO(t)) * vAlpha((int)mIndices(q, 0)) + vO(t) * dVarPhi, vLambda((int)mIndices(q, 1)) * vBeta(vSeason(t)), vLogFactorial_K, true) - mLLK(q, t));
        }

      }
      for (j = 0; j < iJ; j++) {
        for (k = 0; k < iK; k++) {
          for (l = 0; l < iL; l++) {
            mu_Z(k, t)     += au[j][k][l][t];
            mu_Alpha(j, t) += au[j][k][l][t];
            mu_Eta(l, t)   += au[j][k][l][t];
            cZ(k, t, l)    += au[j][k][l][t];
          }
        }
      }
    }

    cv_Alpha.zeros();
    cv_Eta.zeros();
    //
    for (t = 1; t < iT; t++) {
      // populate cv
      for (j = 0; j < iQ; j++) {
        for (l = 0; l < iQ; l++) {

          // av[j1][j2][l1][l2][k1][k2] = P(S_t^alpha = j1, S_{t-1}^alpha = j2,
          //                                S_t^eta = l1, S_{t-1}^eta = l2,
          //                                Z_t = k1, Z_{t-1} = k2| Y_{1:T})

          // mGamma_Enlarged(q1, q2) = mOmega((int)mIndices(q2, 1), (int)mIndices(q2, 2)) * mGamma_Eta((int)mIndices(q1, 2), (int)mIndices(q2, 2)) *
          // mGamma_Alpha((int)mIndices(q1, 0), (int)mIndices(q2, 0));

          av[(int)mIndices(j, 0)]
          [(int)mIndices(l, 0)]
          [(int)mIndices(j, 2)]
          [(int)mIndices(l, 2)]
          [(int)mIndices(j, 1)]
          [(int)mIndices(l, 1)][t] = mGamma_OmegaZB(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
        }
      }
      // populate cv_Omega, vc_B
      for (j1 = 0; j1 < iJ; j1++) {
        for (l1 = 0; l1 < iL; l1++) {
          for (k1 = 0; k1 < iK; k1++) {
            for (j2 = 0; j2 < iJ; j2++) {
              for (l2 = 0; l2 < iL; l2++) {
                for (k2 = 0; k2 < iK; k2++) {

                  // cv_Alpha(j1, j2, t) = P(S_t^alpha = j1, S_{t-1}^alpha = j2 | Y_{1:T})
                  // cv_Eta(l1, l2, t) = P(S_t^eta = l1, S_{t-1}^eta = l2 | Y_{1:T})
                  cv_Alpha(j1, j2, t) += av[j1][j2][l1][l2][k1][k2][t];
                  cv_Eta(l1, l2, t)   += av[j1][j2][l1][l2][k1][k2][t];
                }
              }
            }
          }
        }
      }
    }

    ////////////////// M STEP ////////////////////////////////

    //DELTA
    for (q = 0; q < iQ; q++) {
      vDelta_OmegaZB(q) = au[(int)mIndices(q, 0)][(int)mIndices(q, 1)][(int)mIndices(q, 2)][0];
    }

    vDelta_OmegaZB = vDelta_OmegaZB/accu(vDelta_OmegaZB);

    // Gamma_omega
    for (j = 0; j < iJ; j++) {
      for (l = 0; l < iJ ; l++) {
        mGamma_Alpha_Next(j, l) = accu(cv_Alpha.tube(l, j));
      }
    }
    // Gamma_B
    for (j = 0; j < iL; j++) {
      for (l = 0; l < iL; l++) {
        mGamma_Eta_Next(j, l) = accu(cv_Eta.tube(l, j));
      }
    }
    // Normalize Gamma_omega
    for (j = 0; j < iJ; j++) {
      mGamma_Alpha_Next.row(j) = mGamma_Alpha_Next.row(j)/sum(mGamma_Alpha_Next.row(j));
    }
    // Normalize Gamma_B
    for (j = 0; j < iL; j++) {
      mGamma_Eta_Next.row(j) = mGamma_Eta_Next.row(j)/sum(mGamma_Eta_Next.row(j));
    }

    // mOmega_Next  = OmegaCheck(mOmega_Next, iK, iJ);

    //LAMBDA
    vLambda_Num_Foo.zeros();
    vLambda_Den_Foo.zeros();

    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        vLambda_Num_Foo((int)mIndices(q, 1)) +=  (mu_all(q, t) * mEta(q, t));
        vLambda_Den_Foo((int)mIndices(q, 1)) +=  (mu_all(q, t) * vBeta(vSeason(t)));
      }
    }

    vLambda_Next = vLambda_Num_Foo/vLambda_Den_Foo;

    //BETA
    if (bSeasonal) {
      vBeta_Num_Foo.zeros();
      vBeta_Den_Foo.zeros();

      for (t = 0; t < iT; t++) {
        for (q = 0; q < iQ; q++) {
          vBeta_Num_Foo(vSeason(t)) += (mu_all(q, t) * mEta(q, t));
          vBeta_Den_Foo(vSeason(t)) += (mu_all(q, t) * vLambda_Next((int)mIndices(q, 1)));
        }
      }

      vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    }

    //ALPHA
    vAlpha_Num_Foo.zeros();
    vAlpha_Den_Foo.zeros();
    dVarPhi_Num = 0.0;
    dVarPhi_Den = 0.0;

    for (t = 0; t < iT; t++) {
      if (vO(t) > 0)  {
        for (q = 0; q < iQ; q++) {
          dVarPhi_Num += (mu_all(q, t) * (vY(t) - mEta(q, t)));
          if (t == 0) {
            dVarPhi_Den += (mu_all(q, t) * dY0 );
          } else {
            dVarPhi_Den += (mu_all(q, t) * vY(t - 1));
          }
        }
      } else  {
        for (q = 0; q < iQ; q++) {
          vAlpha_Num_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * (vY(t) - mEta(q, t)));
          if (t == 0) {
            vAlpha_Den_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * dY0 );
          } else {
            vAlpha_Den_Foo((int)mIndices(q, 0)) += (mu_all(q, t) * vY(t - 1));
          }
        }
      }

    }

    vAlpha_Next = vAlpha_Num_Foo/vAlpha_Den_Foo;
    dVarPhi_Next = dVarPhi_Num/dVarPhi_Den;

    if (dVarPhi_Next < 1e-7) {
      dVarPhi_Next = dVarPhi;
    }
    for (j = 0; j < iJ; j++) {
      if(vAlpha_Next(j) < 1e-7) {
        vAlpha_Next(j) = vAlpha(j);
      }
    }

    // Store the llk
    LLKSeries(iter) = dLLK;

    if(iter > 10) {
      eps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }

    iter += 1;

    // Update Parameters
    if (!bFilter) {
      mGamma_Eta   = mGamma_Eta_Next;
      mGamma_Alpha = mGamma_Alpha_Next;
      vLambda = vLambda_Next;
      vAlpha = vAlpha_Next;
      dVarPhi = dVarPhi_Next;
      if (bSeasonal) {
        vBeta = vBeta_Next;
      }
    }

  }

  LLKSeries = LLKSeries.subvec(0, iter - 1);

  arma::mat PredictedProb(iT + 1, iQ);
  arma::mat FilteredProb(iT, iQ);
  arma::mat SmoothedProb(iT, iQ);

  double dC;
  double dLLK_foo;

  PredictedProb.row(0) = vDelta_OmegaZB.t();

  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb.row(t + 1) = FilteredProb.row(t) * mGamma_OmegaZB;
    SmoothedProb.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);
  }

  for (j1 = 0; j1 < iJ; j1++) {
    for (j2 = 0; j2 < iJ; j2++) {
      for (l1 = 0; l1 < iL; l1++) {
        for (l2 = 0; l2 < iL; l2++) {
          for (k1 = 0; k1 < iK; k1++) {
            for (k2 = 0; k2 < iK; k2++) {
              delete av[j1][j2][l1][l2][k1][k2];
            }
          }
        }
      }
    }
  }

  delete[] av;
  for (j1 = 0; j1 < iJ; j1++) {
    for (k1 = 0; k1 < iK; k1++) {
      for (l1 = 0; l1 < iL; l1++) {
        delete au[j1][k1][l1];
      }
    }
  }

  delete[] au;

  for (k1 = 0; k1 < iK; k1++) {
    for (j1 = 0; j1 < iJ; j1++) {
      delete aLLK[k1][j1];
    }
  }

  delete[] aLLK;


  List lOut;

  lOut["mGamma_Eta"] = mGamma_Eta;
  lOut["mGamma_Alpha"] = mGamma_Alpha;
  lOut["mOmega"] = mOmega;
  lOut["vLambda"] = vLambda;
  lOut["vAlpha"] = vAlpha;
  lOut["dVarPhi"] = dVarPhi;
  lOut["vBeta"] = vBeta;

  lOut["mLLK"] = mLLK;

  lOut["LLKSeries"] = LLKSeries;
  lOut["eps"] = eps;

  lOut["mu_Alpha"] = mu_Alpha;
  lOut["mu_Z"]     = mu_Z;
  lOut["mu_Eta"]     = mu_Eta;
  lOut["mu_all"]     = mu_all;
  lOut["mEta"] = mEta;

  lOut["mIndices"] = mIndices;
  lOut["vSeason"]       = vSeason;
  lOut["vO"]       = vO;

  lOut["mGamma_OmegaZB"] = mGamma_OmegaZB;
  lOut["vDelta_OmegaZB"] = vDelta_OmegaZB;
  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothedProb"] = SmoothedProb;

  lOut["dLLK"] = dLLK;

  lOut["mlalpha"] = mlalpha;
  lOut["mlbeta"] = mlbeta;

  lOut["mEta"] = mEta;

  lOut["vY"] = vY;
  lOut["dY0"] = dY0;
  lOut["dO0"] = dO0;
  lOut["iSeason0"] = iSeason0;

  return lOut;

}


//[[Rcpp::export]]
List InitializeBeta_INAR(arma::vec vY, arma::vec vBeta, arma::vec vSeason, arma::vec vLogFactorial_K, int iMaxiter = 1e3, double dTol = 1e-4) {

  int iT = vY.size();
  int t;
  int iD = vBeta.size();

  double dMean = accu(vY)/(iT * 1.0);

  double dAlpha = accu((vY.subvec(0, iT - 2) - dMean) % (vY.subvec(1, iT - 1) - dMean))/accu(pow(vY - dMean, 2.0));

  if (dAlpha < 0.0) {
    dAlpha = 0.1;
  }

  double dLambda = dMean * (1.0 - dAlpha);

  arma::vec vBeta_Next = vBeta;
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
    vLLK(t) = dPB(vY(t), vY(t - 1), dAlpha, dLambda * vBeta(vSeason(t)), vLogFactorial_K, true);
    dLLK   += vLLK(t);
  }

  LLKSeries(0) = dLLK;

  int iter = 1;
  double dEps = 10;

  arma::vec vBeta_Num_Foo(iD);
  arma::vec vBeta_Den_Foo(iD);

  while (dEps > dTol && iter < iMaxiter) {

    //E step
    for (t = 1; t < iT; t++) {
      vEps(t) = exp(log(dLambda * vBeta(vSeason(t))) + dPB(vY(t) - 1.0, vY(t - 1), dAlpha, dLambda * vBeta(vSeason(t)), vLogFactorial_K, true) - vLLK(t));
    }

    //M Step
    vBeta_Num_Foo.zeros();
    vBeta_Den_Foo.zeros();
    //
    dAlpha_Next  = (dSum_t - accu(vEps))/dSum_tm1;
    dLambda_Next = accu(vEps)/(iT * 1.0 - 1.0);

    for (t = 1; t < iT; t++) {
      vBeta_Num_Foo(vSeason(t)) += (vEps(t));
      vBeta_Den_Foo(vSeason(t)) += (dLambda);
    }

    vBeta_Next = vBeta_Num_Foo/vBeta_Den_Foo;

    //update Parameters
    dLambda = dLambda_Next;
    dAlpha  = dAlpha_Next;

    vBeta = vBeta_Next;

    dLLK = 0.0;
    for (t = 1; t < iT; t++) {
      vLLK(t) = dPB(vY(t), vY(t - 1), dAlpha, dLambda * vBeta(vSeason(t)), vLogFactorial_K, true);
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
  lOut["vBeta"] = vBeta;
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;

  return lOut;

}

// mOmega is iL times iE

//[[Rcpp::export]]
List Sim_MSMixInarTwoChains(int iT, arma::mat mOmega, arma::vec vLambda, arma::vec vAlpha,
                                 arma::mat mGamma_Alpha, arma::mat mGamma_Eta, arma::vec vBeta, arma::mat mD) {

  int t;
  int d;
  int iFoo;
  int iA = mGamma_Alpha.n_cols;
  int iE = mGamma_Eta.n_cols;
  int iD = mD.n_cols;

  mD = mD * 1.0;
  arma::vec vSeason(iT);

  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }

  arma::vec vY(iT);
  arma::vec vEta(iT);

  arma::vec vState_Alpha(iT);
  arma::vec vState_Eta(iT);

  arma::vec vDelta_Alpha = getDelta(mGamma_Alpha, iA);
  arma::vec vDelta_Eta = getDelta(mGamma_Eta, iE);

  vState_Alpha(0) = rando_index(vDelta_Alpha);
  vState_Eta(0) = rando_index(vDelta_Eta);

  iFoo = rando_index(mOmega.col(vState_Eta(0)));
  vEta(0) = Rf_rpois(vLambda(iFoo)/(1.0 - vAlpha(vState_Alpha(0))));

  vY(0) = vEta(0);

  for (t = 1; t < iT; t++){

    vState_Alpha(t) = rando_index(mGamma_Alpha.row(vState_Alpha(t - 1)).t());
    vState_Eta(t) = rando_index(mGamma_Eta.row(vState_Eta(t - 1)).t());

    iFoo = rando_index(mOmega.col(vState_Eta(t - 1)));
    vEta(t) = Rf_rpois(vLambda(iFoo) * vBeta(vSeason(t)));
    vY(t) = Rf_rbinom(vY(t - 1), vAlpha(vState_Alpha(t))) + vEta(t);

  }

  List lOut;
  lOut["vState_Alpha"] = vState_Alpha;
  lOut["vState_Eta"] = vState_Eta;
  lOut["vY"] = vY;
  lOut["vEta"] = vEta;

  return lOut;

}

//[[Rcpp::export]]
List Sim_MSMixInarTwoChains_AlphaCor(int iT, arma::mat mOmega, arma::vec vLambda, arma::vec vAlpha,
                            arma::mat mGamma_Alpha, arma::mat mGamma_Eta, arma::vec vBeta,
                            double dVarPhi,
                            arma::vec vO,
                            arma::vec vSeason) {

  int t;
  int iFoo;
  int iA = mGamma_Alpha.n_cols;
  int iE = mGamma_Eta.n_cols;

  arma::vec vY(iT);
  arma::vec vEta(iT);
  arma::vec vB(iT);

  arma::vec vState_Alpha(iT);
  arma::vec vState_Eta(iT);

  arma::vec vDelta_Alpha = getDelta(mGamma_Alpha, iA);
  arma::vec vDelta_Eta = getDelta(mGamma_Eta, iE);

  iFoo = rando_index(mOmega.col(rando_index(vDelta_Eta)));
  double dY0 = Rf_rpois(vLambda(iFoo));

  vState_Alpha(0) = rando_index(vDelta_Alpha);
  vState_Eta(0) = rando_index(vDelta_Eta);

  iFoo = rando_index(mOmega.col(vState_Eta(0)));
  // vEta(0) = Rf_rpois(vLambda(iFoo)/(1.0 - vAlpha(vState_Alpha(0))));
  vEta(0) = Rf_rpois(vLambda(iFoo));

  vB(0) = Rf_rbinom(dY0, (1 - vO(0)) * vAlpha(vState_Alpha(0)) + vO(0) * dVarPhi);
  vY(0) = vB(0) + vEta(0);

  for (t = 1; t < iT; t++){

    vState_Alpha(t) = rando_index(mGamma_Alpha.row(vState_Alpha(t - 1)).t());
    vState_Eta(t) = rando_index(mGamma_Eta.row(vState_Eta(t - 1)).t());

    iFoo = rando_index(mOmega.col(vState_Eta(t - 1)));
    vEta(t) = Rf_rpois(vLambda(iFoo) * vBeta(vSeason(t)));
    vB(t) = Rf_rbinom(vY(t - 1), (1 - vO(t)) * vAlpha(vState_Alpha(t)) + vO(t) * dVarPhi);
    vY(t) = vB(t) + vEta(t);
  }

  List lOut;
  lOut["vState_Alpha"] = vState_Alpha;
  lOut["vState_Eta"] = vState_Eta;
  lOut["vY"] = vY;
  lOut["vB"] = vB;
  lOut["vEta"] = vEta;

  return lOut;

}


//[[Rcpp::export]]
List Sim_MSMixInarOneChain_AlphaCor(int iT, arma::vec vLambda, arma::vec vAlpha,
                                     arma::mat mGamma, arma::vec vBeta,
                                     double dVarPhi,
                                     arma::vec vO,
                                     arma::vec vSeason) {

  int t;
  int iFoo;
  int iQ = mGamma.n_cols;

  arma::vec vY(iT);
  arma::vec vEta(iT);
  arma::vec vB(iT);

  arma::vec vState(iT);

  arma::vec vDelta = getDelta(mGamma, iQ);

  iFoo = rando_index(vDelta);
  double dY0 = Rf_rpois(vLambda(iFoo));

  vState(0) = rando_index(mGamma.row(iFoo).t());

  vEta(0) = Rf_rpois(vLambda(vState(0)));

  vB(0) = Rf_rbinom(dY0, (1 - vO(0)) * vAlpha(vState(0)) + vO(0) * dVarPhi);
  vY(0) = vB(0) + vEta(0);

  for (t = 1; t < iT; t++){

    vState(t) = rando_index(mGamma.row(vState(t - 1)).t());

    vEta(t) = Rf_rpois(vLambda(vState(t)) * vBeta(vSeason(t)));
    vB(t) = Rf_rbinom(vY(t - 1), (1 - vO(t)) * vAlpha(vState(t)) + vO(t) * dVarPhi);
    vY(t) = vB(t) + vEta(t);
  }

  List lOut;
  lOut["vState"] = vState;
  lOut["vY"] = vY;
  lOut["vB"] = vB;
  lOut["vEta"] = vEta;

  return lOut;

}


//[[Rcpp::export]]
List Sim_INAR(int iT, double dAlpha, double dLambda) {

  int t;

  arma::vec vY(iT);
  arma::vec vEta(iT);
  arma::vec vB(iT);

  vEta(0) = Rf_rpois(dLambda/(1.0 - dAlpha));

  vY(0) = vEta(0);
  vB(0) = 0.0;

  for (t = 1; t < iT; t++){

    vEta(t) = Rf_rpois(dLambda);
    vB(t) = Rf_rbinom(vY(t - 1), dAlpha);
    vY(t) = vB(t) + vEta(t);
  }

  List lOut;
  lOut["vY"] = vY;
  lOut["vB"] = vB;
  lOut["vEta"] = vEta;

  return lOut;

}




//[[Rcpp::export]]
double dE1(arma::mat mGamma_OmegaZB, arma::mat mIndices, arma::vec vLambda, arma::vec vAlpha) {

  int iQ = mIndices.n_rows;
  arma::vec vDelta = getDelta(mGamma_OmegaZB, iQ);

  int q;

  double dLambda_bar = 0.0;
  double dAlpha_bar = 0.0;

  for (q = 0; q < iQ; q++) {

    dLambda_bar += (vDelta(q) * vLambda((int)mIndices(q, 1)));
    dAlpha_bar += (vDelta(q) * vAlpha((int)mIndices(q, 0)));

  }

  double dE1 = dLambda_bar/(1.0 - dAlpha_bar);

  return dE1;
}

//[[Rcpp::export]]
double dE2(arma::mat mGamma_OmegaZB, arma::mat mIndices, arma::vec vLambda, arma::vec vAlpha) {

  double dMean = dE1(mGamma_OmegaZB, mIndices, vLambda, vAlpha);

  double dOne = 0.0;
  double dAlpha2_bar = 0.0;
  double dLambda_op_Lambda_bar = 0.0;

  int iQ = mIndices.n_rows;
  arma::vec vDelta = getDelta(mGamma_OmegaZB, iQ);

  int q;

  for (q = 0; q < iQ; q++) {

    dOne += (vDelta(q) * (vAlpha((int)mIndices(q, 0)) * (1.0 - vAlpha((int)mIndices(q, 0))) + 2.0 * vAlpha((int)mIndices(q, 0)) * vLambda((int)mIndices(q, 1))));
    dAlpha2_bar += (vDelta(q) * pow(vAlpha((int)mIndices(q, 0)), 2.0));
    dLambda_op_Lambda_bar += (vDelta(q) * vLambda((int)mIndices(q, 1)) * (1.0 + vLambda((int)mIndices(q, 1))));

  }

  double dE2 = (dMean * dOne + dLambda_op_Lambda_bar)/(1.0 - dAlpha2_bar);

  return dE2;
}

//[[Rcpp::export]]
double dFoo(arma::mat mGamma_OmegaZB, arma::mat mIndices, arma::vec vLambda, arma::vec vAlpha) {

  double dMean = dE1(mGamma_OmegaZB, mIndices, vLambda, vAlpha);

  int iQ = mIndices.n_rows;
  arma::vec vDelta = getDelta(mGamma_OmegaZB, iQ);

  int q;
  int i;

  double dFoo = 0.0;

  for (q = 0; q < iQ; q++) {
    for (i = 0; i < iQ; i++) {

      dFoo += (vLambda((int)mIndices(q, 1))  * vAlpha((int)mIndices(q, 0)) *
        vDelta(i) *  mGamma_OmegaZB(i, q) * (vAlpha((int)mIndices(i, 0))* dMean + vLambda((int)mIndices(i, 1))));
      // vDelta(i) *  mGamma_OmegaZB(q, i) * (vAlpha((int)mIndices(i, 0))* dMean + vLambda((int)mIndices(i, 1))));


    }
  }

  return dFoo;
}


