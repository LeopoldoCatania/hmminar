#ifndef UTILS_H
#define UTILS_H
double MixtDensityScale(arma::vec vOmega_log, arma::vec vD_log);
double LogSumExp(arma::vec vLog);
double factorial2(double dX);
double lfactorial2(double dX);
double abs3(double x);
int rando_index(arma::vec p);
arma::vec getDelta(arma::mat gamma, int m);
arma::mat AccessListMat(Rcpp::List list, std::string element_name);
int WhichMax(arma::vec vX);
double LW(double x, double tolerance = 1e-10, int maxit = 50) ;
arma::mat OmegaZB_Indices(int iJ, int iK, int iL);
double LogOfSum(arma::vec vLog);
arma::mat BoundGamma(arma::mat mGamma, double probbound);
arma::vec BoundProbs(arma::vec vP, double probbound);
#endif
