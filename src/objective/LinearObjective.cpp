#include <cassert>
#include "LinearObjective.h"
#include <iostream>

namespace SAM {

  //TODO: |beta| > lambda / a

  LinearObjective::LinearObjective(const double *xmat, const double *y, int n, int d, int p, bool include_intercept)
    : ObjFunction(xmat, y, n, d, p), XX(d) {
    r.resize(n);

    if (include_intercept) {
      double avr_y = Y.sum()/n;
      model_param.intercept = avr_y;
    }

    for (int i = 0; i < d; i++) {
      XX[i].resize(p, p);
      XX[i] = X[i].transpose() * X[i] / n;
    }

    r = Y;
    update_auxiliary();

    deviance = fabs(eval());
  }

  VectorXd LinearObjective::coordinate_descent(RegFunction *regfunc,
                                               int idx) {
    Eigen::MatrixXd beta_old = model_param.beta[idx];
    Eigen::MatrixXd tmp = X[idx].transpose() * (r + X[idx] * model_param.beta[idx]) / n;

    model_param.beta[idx] = regfunc->threshold(tmp) * n;

    r = r - X[idx] * (model_param.beta[idx] - beta_old);
    return model_param.beta[idx];
  }

  void LinearObjective::intercept_update() {
    double sum_r = r.sum();
    model_param.intercept = sum_r / n;
  }
  void LinearObjective::update_auxiliary() {
    for (int idx = 0; idx < d; idx++)
      update_gradient(idx);
  }

  void LinearObjective::update_gradient(int idx) {
    gr[idx] = X[idx].transpose() * r / n;
  }

  double LinearObjective::get_local_change(const VectorXd& old, int idx) {
    VectorXd tmp = old - model_param.beta[idx];
    return tmp.transpose() * XX[idx] * tmp;
  }
  double LinearObjective::get_local_change_intercept(double old) {
    double tmp = old - model_param.intercept;
    return fabs(tmp);
  }

  double LinearObjective::eval() {
    double v = 0.0;
    VectorXd pred;
    pred.resize(n);
    for (int i = 0; i < n; i++)
      pred(i) = model_param.intercept;
    for (int i = 0; i < d; i++) {
      pred += X[i] * model_param.beta[i];
    }
    for (int i = 0; i < n; i++) {
      v += sqr(Y(i) - pred(i));
    }
    v = v / n;
    return v;
  }

  double LinearObjective::get_r2() {
    return r.dot(r);
  }


}  // namespace picasso
