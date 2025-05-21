#ifndef EMMIXPOIS_H
#define EMMIXPOIS_H

Rcpp::List EM_Mixture_Pois(arma::vec vY, int iJ, int maxIter = 1e3, double tol = 1e-8);

#endif
