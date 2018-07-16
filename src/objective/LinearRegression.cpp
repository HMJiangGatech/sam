#include <cassert>
#include "objective.hpp"

namespace SAM {

  // NOTE: needing some configuration for using STL
  inline double sqr(double x) {
    return x * x;
  }

  LinearRegressionObjective::LinearRegressionObjective(const double *xmat, const double *y, int n, int d, int p, bool include_intercept)
    : ObjFunction(d, p), XX(d) {
    this->d = d;
    this->n = n;
    Y.resize(n);
    for (int i = 0; i < d; i++)
      gr[i].resize(p);

    Xb.resize(n);
    Xb.setZero();

    for (int i = 0; i < n; i++) Y[i] = y[i];

    for (int k = 0; k < d; k++) {
      X[k].resize(n, p);
      for (int j = 0; j < p; j++) {
        for (int i = 0; i < n; i++) X[k](i, j) = xmat[k*p*n+j*n+i];
      }
    }

    r.resize(n);

    if (include_intercept) {
      double avr_y = Y.sum()/n;
      model_param.intercept = avr_y;
    }

    r = Y;
    update_auxiliary();

    for (int i = 0; i < d; i++) {
      XX[i].resize(p, p);
      XX[i] = X[i].transpose() * X[i] / n;
    }

    // saturated fvalue = 0
    deviance = fabs(eval());
  }

  VectorXd LinearRegressionObjective::coordinate_descent(RegFunction *regfunc,
                                                       int idx) {
    Eigen::MatrixXd beta_old = model_param.beta[idx];
    Eigen::MatrixXd tmp = X[idx].transpose() * (r + X[idx] * model_param.beta[idx]);

    model_param.beta[idx] = regfunc->threshold(tmp);

    r = r - X[idx] * (model_param.beta[idx] - beta_old);
    return model_param.beta[idx];
  }

  void LinearRegressionObjective::intercept_update() {
    double sum_r = r.sum();
    model_param.intercept = sum_r / n;
  }
  void LinearRegressionObjective::update_auxiliary() {
    for (int idx = 0; idx < d; idx++)
      update_gradient(idx);
  }

  void LinearRegressionObjective::update_gradient(int idx) {
    gr[idx] = X[idx].transpose() * r / n;
  }

  double LinearRegressionObjective::get_local_change(VectorXd old, int idx) {
    VectorXd tmp = old - model_param.beta[idx];
    return tmp.transpose() * XX[idx] * tmp;
  }
  double LinearRegressionObjective::get_local_change_intercept(double old) {
    double tmp = old - model_param.intercept;
    return fabs(tmp);
  }

  double LinearRegressionObjective::eval() {
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

  double LinearRegressionObjective::get_r2() {
    double res = 0;
    for (int i = 0; i < n; i++) {
      res += r(i) * r(i);
    }
    return res;
  }


}  // namespace picasso
