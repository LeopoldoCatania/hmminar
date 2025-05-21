#include <RcppArmadillo.h>
#include "Utils.h"
#include "SafeFunctions.h"

using namespace arma;
using namespace Rcpp;

const double dZeroBound = 1e-20;

//[[Rcpp::export]]
List FFBS(arma::mat allprobs, arma::vec delta, arma::mat mGamma, int iJ, int iT) {

  int i;
  int j;

  for (i = 0; i < iT; i++) {
    for (j = 0; j < iJ; j++) {
      if (allprobs(i, j) < dZeroBound) {

        allprobs(i, j) = dZeroBound;

      }
    }
  }

  arma::mat lalpha = zeros(iJ, iT);
  arma::mat lbeta  = zeros(iJ, iT);

  arma::vec foo(iJ);
  double sumfoo,lscale;

  foo    = delta % allprobs.row(0).t();
  sumfoo = sum(foo);
  lscale = log(sumfoo);
  foo    = foo/sumfoo ;

  lalpha.col(0) = log(foo)+lscale;
  for (i = 1; i < iT; i++) {
    foo           = (foo.t() * mGamma).t() % allprobs.row(i).t();
    sumfoo        = sum(foo);
    lscale        = lscale + log(sumfoo);
    foo           = foo/sumfoo;
    lalpha.col(i) = log(foo) + lscale;
  }

  for (i = 0; i < iJ; i++) {
    foo(i) = 1.0/iJ;
  }
  lscale = log(iJ*1.0);
  for (i = iT - 2; i >= 0; i--) {
    foo          = mGamma * (allprobs.row(i + 1).t() % foo);
    lbeta.col(i) = log(foo) + lscale;
    sumfoo       = sum(foo);
    foo          = foo / sumfoo;
    lscale       = lscale + log(sumfoo);
  }

  List FS;
  FS["mlalpha"] = lalpha;
  FS["mlbeta"]  = lbeta;

  return FS;
}

//[[Rcpp::export]]
List FFBS_robust(arma::mat lallprobs, arma::vec delta, arma::mat mGamma, int iJ, int iT) {

  int i;
  int j;

  mGamma = GammaCheck(mGamma, iJ);
  delta = DeltaCheck(delta, iJ);

  arma::vec ldelta = log(delta);
  arma::mat lGamma = log(mGamma);

  arma::mat lalpha = zeros(iJ, iT);
  arma::mat lbeta  = zeros(iJ, iT);

  arma::vec foo(iJ);
  double lscale;

  arma::vec lfoo(iJ);
  double lsumfoo;

  lfoo = ldelta + lallprobs.row(0).t();
  lsumfoo = LogOfSum(lfoo);
  lscale = lsumfoo;
  lfoo = lfoo - lsumfoo;

  lalpha.col(0) = lfoo + lscale;

  arma::vec vbaz(iJ);

  for (i = 1; i < iT; i++) {

    for (j = 0; j < iJ; j++) {
      vbaz(j) = LogOfSum(lfoo + lGamma.col(j));
    }
    lfoo = vbaz + lallprobs.row(i).t();
    lsumfoo       = LogOfSum(lfoo);
    lscale        = lscale + lsumfoo;
    lfoo           = lfoo - lsumfoo;
    lalpha.col(i) = lfoo + lscale;
  }

  lfoo.fill(-log(1.0*iJ));
  lscale = log(iJ*1.0);

  for (i = iT - 2; i >= 0; i--) {

    vbaz = lallprobs.row(i + 1).t() + lfoo;
    for (j = 0; j < iJ; j++) {
      lfoo(j) = LogOfSum(lGamma.row(j).t() + vbaz);
    }

    lbeta.col(i) = lfoo + lscale;
    lsumfoo       = LogOfSum(lfoo);
    lfoo          = lfoo - lsumfoo;
    lscale       = lscale + lsumfoo;
  }

  List FS;
  FS["mlalpha"] = lalpha;
  FS["mlbeta"]  = lbeta;

  return FS;
}
