#ifndef REPARAMETERIZATION_H
#define REPARAMETERIZATION_H

arma::mat Gamma_enlarged(arma::mat mGamma_Omega, arma::vec vOmega);
arma::mat Indices(int iJ, int iL);
arma::vec vDelta_enlarged(arma::vec vDelta, arma::vec vOmega);
arma::mat Gamma_OmegaZB(arma::mat mGamma_Omega, arma::mat mOmega, arma::mat mGamma_B);
arma::mat Gamma_Enlarged_2(arma::mat mGamma_Eta, arma::mat mOmega, arma::mat mGamma_Alpha, arma::mat mIndices);
arma::mat Gamma_Enlarged_2_OmegaChain(arma::mat mGamma_Eta, arma::mat mOmega, arma::mat mIndices) ;
#endif
