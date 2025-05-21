#ifndef POISSONBINOMIAL_H
#define POISSONBINOMIAL_H

double dPB(double dY, double dX, double dAlpha, double dLambda, arma::vec vLogFactorial_K, bool bLog = true);
double dPB2(double dY, double dX, double dAlpha, double dLambda, arma::vec vLogFactorial_K, bool bLog = true);
double pPB(double dY, double dX, double dAlpha, double dLambda, arma::vec vLogFactorial_K, bool bLog = true);
double pMixPB(double dY, double dX, double dAlpha, arma::vec vLambda, arma::vec vOmega, arma::vec vLogFactorial_K, bool bLog = true);
double rMixPB(double dX, double dAlpha, arma::vec vLambda, arma::vec vOmega);
#endif
