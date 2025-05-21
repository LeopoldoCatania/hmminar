#ifndef SUMMIXPOIS_H
#define SUMMIXPOIS_H
double dSumMixPois(double dY, arma::vec vOmega, arma::vec vNu, double dMu, bool bLog = true, int iTrunc = 100);
double dMixPois(double dY, arma::vec vOmega, arma::vec vNu, bool bLog = true) ;
#endif
