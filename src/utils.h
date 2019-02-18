#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
using Eigen::VectorXd;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {
  extern double calc_norm(const VectorXd &x);
  extern double sqr(double x);
}

#endif
