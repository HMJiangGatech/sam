#ifndef UTILS_HPP
#define UTILS_HPP

// #include <Rcpp.h>
// #include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]
#include "eigen3/Eigen/Dense"
using Eigen::VectorXd;

namespace SAM {
  extern double calc_norm(const VectorXd &x);
  extern double sqr(double x);
}

#endif
