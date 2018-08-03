#include "GLMObjective.hpp"
#include "../utils.hpp"
#include <stdio.h>
#include <iostream>

namespace SAM {
  const double eps = 1e-5;
  using Eigen::ArrayXd;

  GLMObjective::GLMObjective(const double *xmat, const double *y, int n, int d, int p,
                             double step_size0, bool include_intercept)
    : ObjFunction(xmat, y, n, d, p), P(n), W(n), R(n), sum_r(0), sum_w(0), step_size0(step_size0) {

    if (include_intercept) {
      double avr_y = Y.sum() / n;

      model_param.intercept = log(avr_y / (1 - avr_y));
    }
  }

  double GLMObjective::calc_loss(int idx, const VectorXd &deltaBeta) {
    ArrayXd Xb_p = Xb + X[idx] * deltaBeta;
    ArrayXd P_p = -(Xb_p + model_param.intercept);
    P_p = P_p.exp();
    P_p = 1 / (P_p + 1.0);

    double v = -Y.dot((Xb_p + model_param.intercept).matrix());

    for (int i = 0; i < n; i++)
      if (P_p[i] > 1e-8) v -= (log(P_p[i]) - model_param.intercept - Xb_p[i]);

    return v;
  }

  VectorXd GLMObjective::coordinate_descent(RegFunction *regfunc, int idx) {

    double step_size = step_size0;

    VectorXd tmp;
    double loss0 = calc_loss(idx, VectorXd(p));
    while (step_size > eps) {
      tmp = model_param.beta[idx] - gr[idx] / step_size;
      VectorXd delta_beta = tmp - model_param.beta[idx];
      if (gr[idx].dot(delta_beta) + 0.5*step_size*(delta_beta).dot(delta_beta)+loss0 <= calc_loss(idx, delta_beta)) {
        step_size *= 0.5;
      } else {
        // if (dbg_counter++ < 10) {
        //   printf("%f %f\n", gr[idx].dot(delta_beta) + 0.5*step_size*(delta_beta).dot(delta_beta)+loss0, calc_loss(idx, delta_beta));
        //   printf("tmp:\n");
        //   std::cout << tmp << std::endl;
        // }
        step_size /= 0.5;
        break;
      }
    }

    VectorXd delta_beta = tmp - model_param.beta[idx];
    printf("step_size:%f %f %f\n", step_size, gr[idx].dot(delta_beta) + 0.5*step_size*(delta_beta).dot(delta_beta)+loss0, calc_loss(idx, tmp-model_param.beta[idx]));

    tmp = (model_param.beta[idx] - gr[idx] / step_size) / n;
    VectorXd old_beta = model_param.beta[idx];
    model_param.beta[idx] = regfunc->threshold_p(tmp, step_size) * n;
    // printf("%f %f %f %f\n", calc_norm(model_param.beta[idx]), calc_norm(gr[idx]), calc_norm(tmp), calc_norm(old_beta));

    delta_beta = model_param.beta[idx] - old_beta;
    if (calc_norm(tmp) > 1e-8) {
      // Xb += delta*X[idx*n]
      Xb += X[idx] * delta_beta;

      // r -= delta*w*X
      R = R - W.cwiseProduct(X[idx] * delta_beta);
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
    gr[idx] = X[idx].transpose() * (Y - P) / n;
  }

  double GLMObjective::get_local_change(const VectorXd &old, int idx) {
    VectorXd delta_beta = old - model_param.beta[idx];
    VectorXd delta_Xb = X[idx] * delta_beta;
    return delta_beta.cwiseProduct(W).dot(delta_beta);
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
    // printf("intercept:%f\n", model_param.intercept);
    // printf("P:\n");
    // std::cout << P << std::endl;
    R = Y - P;
    // printf("R:\n");
    // std::cout << R << std::endl;

    W = P.array() * -(P.array() - 1);
    // printf("W:\n");
    // std::cout << W << std::endl;
    sum_w = W.sum();
  }

  double LogisticObjective::eval() {
    double v = 0.0;

    // for (int i = 0; i < n; i++) v -= Y[i] * (model_param.intercept + Xb[i]);
    v -= Y.dot((Xb.array() + model_param.intercept).matrix());

    for (int i = 0; i < n; i++)
      if (P[i] > 1e-8) v -= (log(P[i]) - model_param.intercept - Xb[i]);

    return (v / n);
  }

  PoissonObjective::PoissonObjective(const double *xmat, const double *y, int n,
                                     int d, int p, double step_size0, bool include_intercept)
    : GLMObjective(xmat, y, n, d, p, step_size0, include_intercept) {
    update_auxiliary();
    for (int i = 0; i < d; i++) update_gradient(i);

    model_param.intercept = 0.0;
    update_auxiliary();

    deviance = fabs(eval());
  };

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
