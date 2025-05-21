#include <RcppArmadillo.h>
#include "SumMixPois.h"
#include "Utils.h"
#include "Neyman.h"
#include "FFBS.h"
#include "SafeFunctions.h"
#include "Reparameterization.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List EM_MSSumMixPois(arma::vec vY, int iJ, int iL,  
                           arma::vec vMu, arma::vec vNu, 
                           arma::vec vOmega,
                           arma::mat mGamma,
                           double dTol = 1e-7, int iMaxiter = 4000) {
  
  int iQ = iL * iJ;
  int iT = vY.size();
  
  int t;
  int j;
  int l;
  int q;
  
  int j1;
  int l1;
  int j2;
  int l2;
  
  arma::cube cLLK(iL, iT, iJ);
  arma::vec LLKSeries(iMaxiter);
  arma::mat mLLK(iQ, iT);
  
  arma::vec vOmega_Next = vOmega;
  arma::mat mGamma_Next = mGamma;
  arma::vec vMu_Next = vMu;
  arma::vec vNu_Next = vNu;
  
  arma::vec vDelta(iQ);
  arma::mat mGamma_enlarged(iQ, iQ);
  
  arma::cube ca(iL, iT, iJ);  //P(S_t, Z_t| Y_{1:T})
  arma::mat  mz(iL, iT);      //P(Z_t| Y_{1:T})
  arma::mat  mu(iJ, iT);      //P(S_t| Y_{1:T})
  arma::cube cv(iJ, iJ, iT);  //P(S_t, S_{t-1}| Y_{1:T})
  
  arma::cube cK(iL, iT, iJ); //E[K_t | Y_t, S_t, Z_t]
  
  // double av [iJ][iJ][iL][iL][iT]; P(S_t, S_{t-1}, Z_t, Z_{t-1} | Y_{1:T})
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iL];
      for (l1 = 0; l1 < iL; l1++) {
        av[j1][j2][l1] = new double*[iL];
        for (l2 = 0; l2 < iL; l2++) {
          av[j1][j2][l1][l2] = new double[iT];
        }
      }
    }
  }
  
  arma::mat mIndices = Indices(iJ, iL);
  
  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob) 
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) 
  
  int iter = 0;
  double dEps = 1.0;
  double dLLK = 0.0;
  double dFoo = 0.0;
  double dBaz = 0.0;
  // 
  while (dEps > dTol && iter < iMaxiter) {
    
    mGamma_enlarged = Gamma_enlarged(mGamma, vOmega);
    vDelta          = getDelta(mGamma_enlarged, iQ);
    
    //calculate densities
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        
        mLLK(q, t) = dNeyman(vY(t), vNu((int)mIndices(q, 1)), vMu((int)mIndices(q, 0)), true);
        
        cLLK((int)mIndices(q, 1), t, (int)mIndices(q, 0)) =  mLLK(q, t) ;
        
      }
    }
    
    //filtering
    fb = FFBS(exp(mLLK.t()), vDelta, mGamma_enlarged, iQ, iT);
    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");
    
    dLLK = LogSumExp(mlalpha.col(iT - 1));
    
    mz.zeros();
    mu.zeros();
    
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        
        ca((int)mIndices(q, 1), t, (int)mIndices(q, 0)) = mlalpha(q, t) + mlbeta(q, t) - dLLK;
        
      }
      
      for (l = 0; l < iL; l++) {
        for(j = 0; j < iJ; j++) {
          
          mz(l, t) += exp(ca(l, t, j));
          mu(j, t) += exp(ca(l, t, j));
          
          dFoo = dNeyman(vY(t) + 1.0, vNu(l), vMu(j), true);
          
          cK(l, t, j) = (vY(t) + 1)/vNu(l) * exp(dFoo - cLLK(l, t, j));
          
        }
      }
      
    }
    cv.zeros();
    //
    for (t = 1; t < iT; t++) {
      // populate cv
      for (j = 0; j < iQ; j++) {
        for (l = 0; l < iQ; l++) {
          av[(int)mIndices(j, 0)]
          [(int)mIndices(l, 0)]
          [(int)mIndices(j, 1)]
          [(int)mIndices(l, 1)][t] = mGamma_enlarged(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
        }
      }
      // populate cv_Omega, vc_B
      for (j1 = 0; j1 < iJ; j1++) {
        for (l1 = 0; l1 < iL; l1++) {
          for (j2 = 0; j2 < iJ; j2++) {
            for (l2 = 0; l2 < iL; l2++) {
              cv(j1, j2, t)     += av[j1][j2][l1][l2][t];
            }
          }
        }
      }
    }
    // Gamma
    for (j = 0; j < iJ; j++) {
      for (l = 0; l < iJ; l++) {
        mGamma_Next(j, l) = accu(cv.tube(l, j));
      }
    }
    // Normalize Gamma
    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/sum(mGamma_Next.row(j));
    }
    // vOmega
    for (l = 0; l < iL; l++) {
      vOmega_Next(l) = accu(mz.row(l).t());
    }
    
    vOmega_Next = vOmega_Next/accu(vOmega_Next);
    
    //Mu
    vMu_Next.zeros();
    for (j = 0; j < iJ; j++) {
      dFoo = 0.0;
      dBaz = 0.0;
      for (t = 0; t < iT; t++) {
        for (l = 0; l < iL; l++) {
          dFoo += (cK(l, t, j) * exp(ca(l, t, j)));
          dBaz += exp(ca(l, t, j));
        }
      }
      vMu_Next(j) = dFoo/dBaz;
    }
    
    //Nu
    vNu_Next.zeros();
    for (l = 0; l < iL; l++) {
      dFoo = 0.0;
      dBaz = 0.0;
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          dFoo += (vY(t) * exp(ca(l, t, j)));
          dBaz += (cK(l, t, j) * exp(ca(l, t, j)));
        }
      }
      vNu_Next(l) = dFoo/dBaz;
    }
    
    vOmega = vOmega_Next;
    mGamma = mGamma_Next;
    vMu    = vMu_Next;
    vNu    = vNu_Next;
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  arma::mat PredictedProb(iT + 1, iQ);
  arma::mat FilteredProb(iT, iQ);
  arma::mat SmoothedProb(iT, iQ);
  
  double dC;
  double dLLK_foo;
  
  PredictedProb.row(0) = vDelta.t();
  
  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb.row(t + 1) = FilteredProb.row(t) * mGamma_enlarged;
    SmoothedProb.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);
  }
  
  delete[] av;
  
  List lOut;
  
  lOut["vNu"] = vNu;
  lOut["vOmega"] = vOmega;
  lOut["vMu"] = vMu;
  lOut["mGamma"] = mGamma;
  
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"] = dLLK;
  lOut["dEps"] = dEps;
  
  lOut["mLLK"] = mLLK;
  lOut["mz"] = mz;
  lOut["mu"] = mu;
  
  lOut["mGamma_enlarged"] = mGamma_enlarged;
  lOut["PredictedProb"]   = PredictedProb;
  lOut["FilteredProb"]    = FilteredProb;
  lOut["SmoothedProb"]    = SmoothedProb;
  
  return lOut;
  
}

//[[Rcpp::export]]
Rcpp::List EM_MSSumMixPois_Dummy(arma::vec vY, int iJ, int iL,  
                           arma::vec vMu, arma::vec vNu, 
                           arma::vec vOmega,
                           arma::mat mGamma,
                           arma::vec vBeta,
                           arma::mat mD,
                           double dTol = 1e-7, int iMaxiter = 4000) {
  
  int iQ = iL * iJ;
  int iT = vY.size();
  int iD = mD.n_cols;
  
  int t;
  int j;
  int l;
  int q;
  int d;
  
  int j1;
  int l1;
  int j2;
  int l2;
  
  arma::cube cLLK(iL, iT, iJ);
  arma::vec LLKSeries(iMaxiter);
  arma::mat mLLK(iQ, iT);
  
  arma::vec vOmega_Next = vOmega;
  arma::mat mGamma_Next = mGamma;
  arma::vec vBeta_Next = vBeta;
  arma::vec vMu_Next = vMu;
  arma::vec vNu_Next = vNu;
  
  arma::vec vDelta(iQ);
  arma::mat mGamma_enlarged(iQ, iQ);
  
  arma::cube ca(iL, iT, iJ);  //P(S_t, Z_t| Y_{1:T})
  arma::mat  mz(iL, iT);      //P(Z_t| Y_{1:T})
  arma::mat  mu(iJ, iT);      //P(S_t| Y_{1:T})
  arma::cube cv(iJ, iJ, iT);  //P(S_t, S_{t-1}| Y_{1:T})
  
  arma::cube cK(iL, iT, iJ); //E[K_t | Y_t, S_t, Z_t]
  
  mD = mD * 1.0;
  arma::vec vSeason(iT);
  
  vSeason.zeros();
  
  for (t = 0; t < iT; t++) {
    for (d = 0; d < iD; d++) {
      if (mD(t, d) == 1) {
        vSeason(t) = d;
      }
    }
  }
  
  // double av [iJ][iJ][iL][iL][iT]; P(S_t, S_{t-1}, Z_t, Z_{t-1} | Y_{1:T})
  double***** av = new double****[iJ];
  for (j1 = 0; j1 < iJ; j1++) {
    av[j1] = new double***[iJ];
    for (j2 = 0; j2 < iJ; j2++) {
      av[j1][j2] = new double**[iL];
      for (l1 = 0; l1 < iL; l1++) {
        av[j1][j2][l1] = new double*[iL];
        for (l2 = 0; l2 < iL; l2++) {
          av[j1][j2][l1][l2] = new double[iT];
        }
      }
    }
  }
  
  arma::mat mIndices = Indices(iJ, iL);
  
  List fb;
  arma::mat mlalpha(iQ, iT); // log alpha (forward prob) 
  arma::mat mlbeta(iQ, iT);  // log beta  (backword prob) 
  
  int iter = 0;
  double dEps = 1.0;
  double dLLK = 0.0;
  double dFoo = 0.0;
  double dBaz = 0.0;
  // 
  while (dEps > dTol && iter < iMaxiter) {
    
    mGamma_enlarged = Gamma_enlarged(mGamma, vOmega);
    vDelta          = getDelta(mGamma_enlarged, iQ);
    
    //calculate densities
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        
        mLLK(q, t) = dNeyman(vY(t), vNu((int)mIndices(q, 1)), vMu((int)mIndices(q, 0)) * vBeta(vSeason(t)), true);
        
        cLLK((int)mIndices(q, 1), t, (int)mIndices(q, 0)) =  mLLK(q, t) ;
        
      }
    }
    
    //filtering
    fb = FFBS(exp(mLLK.t()), vDelta, mGamma_enlarged, iQ, iT);
    mlalpha = AccessListMat(fb, "mlalpha");
    mlbeta  = AccessListMat(fb, "mlbeta");
    
    dLLK = LogSumExp(mlalpha.col(iT - 1));
    
    mz.zeros();
    mu.zeros();
    
    for (t = 0; t < iT; t++) {
      for (q = 0; q < iQ; q++) {
        
        ca((int)mIndices(q, 1), t, (int)mIndices(q, 0)) = mlalpha(q, t) + mlbeta(q, t) - dLLK;
        
      }
      
      for (l = 0; l < iL; l++) {
        for(j = 0; j < iJ; j++) {
          
          mz(l, t) += exp(ca(l, t, j));
          mu(j, t) += exp(ca(l, t, j));
          
          dFoo = dNeyman(vY(t) + 1.0, vNu(l), vMu(j) * vBeta(vSeason(t)), true);
          
          cK(l, t, j) = (vY(t) + 1)/vNu(l) * exp(dFoo - cLLK(l, t, j));
          
        }
      }
      
    }
    cv.zeros();
    //
    for (t = 1; t < iT; t++) {
      // populate cv
      for (j = 0; j < iQ; j++) {
        for (l = 0; l < iQ; l++) {
          av[(int)mIndices(j, 0)]
          [(int)mIndices(l, 0)]
          [(int)mIndices(j, 1)]
          [(int)mIndices(l, 1)][t] = mGamma_enlarged(l, j) * exp(mlalpha(l,t - 1) + mLLK(j, t) + mlbeta(j, t) - dLLK);
        }
      }
      // populate cv_Omega, vc_B
      for (j1 = 0; j1 < iJ; j1++) {
        for (l1 = 0; l1 < iL; l1++) {
          for (j2 = 0; j2 < iJ; j2++) {
            for (l2 = 0; l2 < iL; l2++) {
              cv(j1, j2, t)     += av[j1][j2][l1][l2][t];
            }
          }
        }
      }
    }
    // Gamma
    for (j = 0; j < iJ; j++) {
      for (l = 0; l < iJ; l++) {
        mGamma_Next(j, l) = accu(cv.tube(l, j));
      }
    }
    // Normalize Gamma
    for (j = 0; j < iJ; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j)/sum(mGamma_Next.row(j));
    }
    // vOmega
    for (l = 0; l < iL; l++) {
      vOmega_Next(l) = accu(mz.row(l).t());
    }
    
    vOmega_Next = vOmega_Next/accu(vOmega_Next);
    
    //Mu
    vMu_Next.zeros();
    for (j = 0; j < iJ; j++) {
      dFoo = 0.0;
      dBaz = 0.0;
      for (t = 0; t < iT; t++) {
        for (l = 0; l < iL; l++) {
          dFoo += (cK(l, t, j) * exp(ca(l, t, j)));
          dBaz += exp(ca(l, t, j));
        }
      }
      vMu_Next(j) = dFoo/dBaz;
    }
    
    //Nu
    vNu_Next.zeros();
    for (l = 0; l < iL; l++) {
      dFoo = 0.0;
      dBaz = 0.0;
      for (t = 0; t < iT; t++) {
        for (j = 0; j < iJ; j++) {
          dFoo += (vY(t) * exp(ca(l, t, j)));
          dBaz += (cK(l, t, j) * exp(ca(l, t, j)));
        }
      }
      vNu_Next(l) = dFoo/dBaz;
    }
    
    //beta
    for (d = 0; d < iD; d++) {
      dFoo = 0.0;
      dBaz = 0.0;
      for (t = 0; t < iT; t++) {
        for (q = 0; q < iQ; q++) {
          dFoo += (mD(t, d) * exp(ca((int)mIndices(q, 1), t, (int)mIndices(q, 0))) * cK((int)mIndices(q, 1), t, (int)mIndices(q, 0)));
          dBaz += (mD(t, d) * exp(ca((int)mIndices(q, 1), t, (int)mIndices(q, 0))) * vMu_Next((int)mIndices(q, 0)));
        }
      }
      vBeta_Next(d) = dFoo/dBaz;
    }
    
    vOmega = vOmega_Next;
    mGamma = mGamma_Next;
    vMu    = vMu_Next;
    vNu    = vNu_Next;
    vBeta  = vBeta_Next; 
    
    LLKSeries(iter) = dLLK;
    
    if(iter > 10) {
      dEps = abs3((LLKSeries(iter) - LLKSeries(iter - 1)))/abs3(LLKSeries(iter - 1));
    }
    
    iter += 1;
    
  }
  
  LLKSeries = LLKSeries.subvec(0, iter - 1);
  
  arma::mat PredictedProb(iT + 1, iQ);
  arma::mat FilteredProb(iT, iQ);
  arma::mat SmoothedProb(iT, iQ);
  
  double dC;
  double dLLK_foo;
  
  PredictedProb.row(0) = vDelta.t();
  
  for (t = 0; t < iT; t++) {
    dC            = max(mlalpha.col(t));
    dLLK_foo      = dC + log(sum(exp(mlalpha.col(t) - dC)));
    FilteredProb.row(t) = exp(mlalpha.col(t).t() - dLLK_foo);
    PredictedProb.row(t + 1) = FilteredProb.row(t) * mGamma_enlarged;
    SmoothedProb.row(t) = exp(mlalpha.col(t).t() + mlbeta.col(t).t() - dLLK);
  }
  
  delete[] av;
  
  List lOut;
  
  lOut["vNu"]    = vNu;
  lOut["vOmega"] = vOmega;
  lOut["vMu"]    = vMu;
  lOut["mGamma"] = mGamma;
  lOut["vBeta"]  = vBeta;
  
  lOut["LLKSeries"] = LLKSeries;
  lOut["dLLK"]      = dLLK;
  lOut["dEps"]      = dEps;
  
  lOut["mLLK"] = mLLK;
  lOut["mz"]   = mz;
  lOut["mu"]   = mu;
  
  lOut["mGamma_enlarged"] = mGamma_enlarged;
  lOut["PredictedProb"]   = PredictedProb;
  lOut["FilteredProb"]    = FilteredProb;
  lOut["SmoothedProb"]    = SmoothedProb;
  
  return lOut;
  
}
