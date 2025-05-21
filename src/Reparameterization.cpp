#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;


//mOmega  is J x L
//mGamma  is J x J
arma::mat Gamma_OmegaZ(arma::mat mGamma_Omega, arma::mat mOmega) {

  int iJ = mGamma_Omega.n_cols;
  int iK = mOmega.n_rows;

  arma::vec vOmega =  vectorise(mOmega);

  arma::mat mOnes = ones(iK, iK);
  arma::vec vOnes = ones(iK * iJ);

  arma::mat mGamma_Star = (vOnes * vOmega.t()) % kron(mGamma_Omega, mOnes);

  return mGamma_Star;

}

//vOmega  is 1 x L
//mGamma  is J x J

//[[Rcpp::export]]
arma::mat Gamma_enlarged(arma::mat mGamma_Omega, arma::vec vOmega) {

  int iJ = mGamma_Omega.n_cols;
  int iK = vOmega.size();

  arma::mat mOmega(iK, iJ);
  for (int j = 0; j < iJ; j++) {
    mOmega.col(j) = vOmega;
  }

  arma::vec vOmega_Vec =  vectorise(mOmega);

  arma::mat mOnes = ones(iK, iK);
  arma::vec vOnes = ones(iK * iJ);

  arma::mat mGamma_Star = (vOnes * vOmega_Vec.t()) % kron(mGamma_Omega, mOnes);

  return mGamma_Star;

}

//[[Rcpp::export]]
arma::vec vDelta_enlarged(arma::vec vDelta, arma::vec vOmega) {

  int iJ = vDelta.size();
  int iM = vOmega.size();

  int iQ = iJ * iM;

  arma::vec vDelta_Star(iQ);

  int iC = 0;
  int j;
  int m;

  for (j = 0; j < iJ; j++) {
    for (m = 0; m < iM; m++) {
      vDelta_Star(iC) = vDelta(j) * vOmega(m);
      iC += 1;
    }
  }

  return vDelta_Star;

}

//[[Rcpp::export]]
arma::mat Indices(int iJ, int iL) {

  arma::mat mIndices = zeros(iJ * iL, 2);

  int j,k;
  int iC = 0;

    for(j = 0; j < iJ; j++) {
      for(k = 0; k < iL; k++) {

        mIndices(iC, 0) = j;
        mIndices(iC, 1) = k;

        iC += 1;

      }
    }

  return mIndices;

}

//[[Rcpp::export]]
arma::mat Gamma_Enlarged_2(arma::mat mGamma_Eta, arma::mat mOmega, arma::mat mGamma_Alpha,
                           arma::mat mIndices) {

  int iQ = mIndices.n_rows;
  int q1;
  int q2;

 // //indices in the enlarged reparameterization
 // // First index for S_t^alpha
 // // Second index for Z_t
 // // Third index for S_t^eta
 //  mIndices = OmegaZB_Indices(iJ, iK, iL);

  arma::mat mGamma_Enlarged(iQ, iQ);

  for (q1 = 0; q1 < iQ; q1++) {
    for (q2 = 0; q2 < iQ; q2++) {
      mGamma_Enlarged(q1, q2) = mOmega((int)mIndices(q2, 1), (int)mIndices(q2, 2)) *
        mGamma_Eta((int)mIndices(q1, 2), (int)mIndices(q2, 2)) *
        mGamma_Alpha((int)mIndices(q1, 0), (int)mIndices(q2, 0));
    }
  }

  return mGamma_Enlarged;

}


//[[Rcpp::export]]
arma::mat Gamma_Enlarged_2_OmegaChain(arma::mat mGamma_Eta, arma::mat mOmega, arma::mat mIndices) {

  int iQ = mIndices.n_rows;
  int q1;
  int q2;

  arma::mat mGamma_Enlarged(iQ, iQ);

  for (q1 = 0; q1 < iQ; q1++) {
    for (q2 = 0; q2 < iQ; q2++) {
      mGamma_Enlarged(q1, q2) = mOmega((int)mIndices(q2, 1), (int)mIndices(q2, 0)) * mGamma_Eta((int)mIndices(q1, 0), (int)mIndices(q2, 0));
    }
  }

  return mGamma_Enlarged;

}
