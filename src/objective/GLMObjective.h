#ifndef GLMOBJECTIVE_HPP
#define GLMOBJECTIVE_HPP

#include "objective.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
using Eigen::VectorXd;
using Eigen::MatrixXd;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {
  class GLMObjective : public ObjFunction {
  protected:
    Eigen::VectorXd P, W, R;
    double sum_r, sum_w;
    double step_size0;

  public:
    GLMObjective(const double *xmat, const double *y, int n, int d, int p,
                 double step_size0, bool include_intercept);

    double calc_loss(int idx, const VectorXd &deltaBeta);
    VectorXd coordinate_descent(RegFunction *regfunc, int idx);

    void intercept_update();
    void update_gradient(int);

    double get_local_change(const VectorXd &old, int idx);
    double get_local_change_intercept(double old);
    double get_r2();
  };

  class LogisticObjective : public GLMObjective {
  public:
    LogisticObjective(const double *xmat, const double *y, int n, int d, int p,
                      double step_size0, bool include_intercept = false);

    void update_auxiliary();

    double eval();
  };

  class PoissonObjective : public GLMObjective {
  public:
    PoissonObjective(const double *xmat, const double *y, int n, int d, int p,
                     double step_size0, bool include_intercept = false);

    void update_auxiliary();

    double eval();
  };

}


#endif
