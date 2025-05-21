#ifndef FFBS_H
#define FFBS_H

Rcpp::List FFBS(arma::mat allprobs, arma::vec delta, arma::mat mGamma, int iJ, int iT);
Rcpp::List FFBS_robust(arma::mat lallprobs, arma::vec delta, arma::mat mGamma, int iJ, int iT);
#endif
