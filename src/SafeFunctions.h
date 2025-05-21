#ifndef SAFEFUNCTIONS_H
#define SAFEFUNCTIONS_H

arma::mat GammaCheck(arma::mat mGamma, int iJ);
arma::mat OmegaCheck(arma::mat mOmega, int iK, int iJ);
double LogDensityCheck(double dLLK);
double DensityCheck(double dDensity);
arma::mat SigmaCheck(arma::mat mSigma, int iN, int iK);
arma::mat mKappaCheck(arma::mat mKappa, int iN, int iL) ;
arma::vec DeltaCheck(arma::vec vDelta, int iK);
arma::mat BetaCheck(arma::mat mBeta, int iN, int iD) ;

#endif
