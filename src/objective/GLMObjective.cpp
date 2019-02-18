#include "GLMObjective.h"
#include "../utils.h"
#include <stdio.h>
#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <complex>
#include <algorithm>

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {
  const double eps = 1e-4;
  using Eigen::ArrayXd;
  using Eigen::VectorXcd;

  GLMObjective::GLMObjective(const double *xmat, const double *y, int n, int d, int p,
                             double step_size0, bool include_intercept)
    : ObjFunction(xmat, y, n, d, p), P(n), W(n), R(n), sum_r(0), sum_w(0), step_size0(step_size0) {

    if (include_intercept) {
      double avr_y = Y.sum() / n;

      model_param.intercept = log(avr_y / (1 - avr_y));
    }
  }


  VectorXd GLMObjective::coordinate_descent(RegFunction *regfunc, int idx) {
    VectorXd gr = X[idx].transpose() * -R / n;
    VectorXd tmp;

    MatrixXd wXX(p, p);
    wXX.setZero();
    for (int i = 0; i < n; i++)
      wXX += W[i] * X[idx].row(i).transpose() * X[idx].row(i);
    wXX /= n;

    VectorXcd eigenvalues = wXX.eigenvalues();
    double step_size = 0;
    for (int i = 0; i < eigenvalues.size(); i++) {
      step_size = std::max(step_size, eigenvalues[i].real());
    }
    assert(step_size >= 0);

    tmp = regfunc->threshold_p((model_param.beta[idx] - gr / step_size), step_size);
    VectorXd delta_beta = tmp - model_param.beta[idx];
    tmp = (wXX*model_param.beta[idx] - gr);
    VectorXd old_beta = model_param.beta[idx];
    model_param.beta[idx] = regfunc->threshold(tmp)/step_size;

    delta_beta = model_param.beta[idx] - old_beta;

    if (calc_norm(delta_beta) > 1e-8) {
      Xb += X[idx] * delta_beta;

      R -= W.cwiseProduct(X[idx] * delta_beta);
    }

    return model_param.beta[idx];
  }

  void GLMObjective::intercept_update() {
    sum_r = R.sum();
    model_param.intercept += sum_r/sum_w;
    R -= sum_r/sum_w * W;
    sum_r = 0;
  }

  void GLMObjective::update_gradient(int idx) {
    gr[idx] = X[idx].transpose() * (P - Y) / n;
  }

  double GLMObjective::get_local_change(const VectorXd &old, int idx) {
    VectorXd delta_beta = old - model_param.beta[idx];
    VectorXd delta_Xb = X[idx] * delta_beta;
    MatrixXd wXX(p, p);
    for (int i = 0; i < p; i++)
      for (int j = 0; j < p; j++)
        wXX(i, j) = 0;
    for (int i = 0; i < n; i++)
      wXX += W[i] * X[idx].row(i).transpose() * X[idx].row(i);
    VectorXcd eigenvalues = wXX.eigenvalues();
    double max_eigen = 0;
    for (int i = 0; i < eigenvalues.size(); i++) {
      max_eigen = std::max(max_eigen, eigenvalues[i].real());
    }
    return max_eigen * (delta_beta.dot(delta_beta)) / (2 * n);
  }
  double GLMObjective::get_local_change_intercept(double old) {
    double tmp = old - model_param.intercept;
    return (sum_w * tmp * tmp / (2 * n));
  }
  double GLMObjective::get_r2() {
    // NOTE: not needed
    return 0;
  }

  LogisticObjective::LogisticObjective(const double *xmat, const double *y, int n,
                                       int d, int p, double step_size0, bool include_intercept)
    : GLMObjective(xmat, y, n, d, p, step_size0, include_intercept) {
    update_auxiliary();
    for (int i = 0; i < d; i++) update_gradient(i);

    model_param.intercept = 0.0;
    update_auxiliary();

    deviance = fabs(eval());
  }

  void LogisticObjective::update_auxiliary() {
    P = -(Xb.array() + model_param.intercept);
    P = P.array().exp();
    P = (P.array() + 1.0).inverse();
    R = Y - P;

    W = P.array() * -(P.array() - 1);
    sum_w = W.sum();
  }

  double LogisticObjective::eval() {
    double v = 0.0;

    v -= Y.dot((Xb.array() + model_param.intercept).matrix());

    for (int i = 0; i < n; i++)
      if (P[i] > 1e-8) v -= (log(P[i]) - model_param.intercept - Xb[i]);

    return (v / n);
  }

  PoissonObjective::PoissonObjective(const double *xmat, const double *y, int n,
                                     int d, int p, double step_size0, bool include_intercept)
    : GLMObjective(xmat, y, n, d, p, step_size0, include_intercept) {

    model_param.intercept = 0.0;
    update_auxiliary();
    for (int i = 0; i < d; i++) update_gradient(i);


    deviance = fabs(eval());
  }

  void PoissonObjective::update_auxiliary() {
    P = Xb.array() + model_param.intercept;
    P = P.array().exp();
    R = Y - P;
    W = P;
    sum_w = W.sum();
  }

  double PoissonObjective::eval() {
    double v = 0.0;
    for (int i = 0; i < n; i++)
      v = v + P[i] - Y[i] * (model_param.intercept + Xb[i]);
    return (v / n);
  }

}
