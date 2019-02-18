#ifndef SAM_OBJECTIVE_H
#define SAM_OBJECTIVE_H

#include <cmath>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include "../utils.h"

#include <ctime>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {

  class ModelParam {
  public:
    int d;
    int p;
    vector<VectorXd> beta;
    double intercept;

    ModelParam(int _d, int _p):d(_d), p(_p), beta(_d) {
      for (int i = 0; i < d; i++) {
        beta[i].resize(p);
        beta[i].setZero();
      }
      intercept = 0.0;
    }
  };

  class RegFunction {
  public:
    virtual double threshold(double x) = 0;
    virtual VectorXd threshold(const VectorXd& x) = 0;
    virtual VectorXd threshold_p(const VectorXd &x, double a) = 0;
    virtual void set_param(double lambda, double gamma) = 0;
    virtual double get_lambda() = 0;

    virtual ~RegFunction(){};

    double threshold_l1(double x, double thr) {
      if (x > thr)
        return x - thr;
      else if (x < -thr)
        return x + thr;
      else
        return 0;
    }
    VectorXd threshold_l1(VectorXd x, double thr) {
      double norm = calc_norm(x);
      if (norm <= thr) {
        for (int i = 0; i < (int)x.size(); i++)
          x[i] = 0;
      } else {
        double ratio = (norm - thr) / norm;
        assert(ratio > 0);
        for (int i = 0; i < (int)x.size(); i++)
          x[i] *= ratio;
      }
      return x;
    }
  };

  class RegL1 : public RegFunction {
  private:
    double m_gamma;
    double m_lambda;

  public:
    void set_param(double lambda, double gamma) { m_lambda = lambda; m_gamma = gamma;}
    double get_lambda() { return m_lambda; };
    double threshold(double x) { return threshold_l1(x, m_lambda); }
    VectorXd threshold(const VectorXd &x) {
      return threshold_l1(x, m_lambda);
    }
    VectorXd threshold_p(const VectorXd &x, double a) {
      return threshold_l1(x, m_lambda/a);
    }
  };

  class RegSCAD : public RegFunction {
  private:
    double m_lambda;
    double m_gamma;

  public:
    void set_param(double lambda, double gamma) {
      m_lambda = lambda;
      m_gamma = gamma;
    };
    double get_lambda() { return m_lambda; };

    double threshold(double x) {
      if (fabs(x) > fabs(m_gamma * m_lambda)) {
        return x;
      } else {
        if (fabs(x) > fabs(2 * m_lambda)) {
          return threshold_l1(x, m_gamma * m_lambda / (m_gamma - 1)) /
            (1 - 1 / (m_gamma - 1));
        } else {
          return threshold_l1(x, m_lambda);
        }
      }
    };
  };

  class RegMCP : public RegFunction {
  private:
    double m_lambda;
    double m_gamma;

  public:
    void set_param(double lambda, double gamma) {
      m_lambda = lambda;
      m_gamma = gamma;
    }
    double get_lambda() { return m_lambda; };

    double threshold(double x) {
      if (fabs(x) > fabs(m_gamma * m_lambda)) {
        return x;
      } else {
        if (fabs(x) > fabs(2 * m_lambda)) {
          return threshold_l1(x, m_gamma * m_lambda / (m_gamma - 1)) /
            (1 - 1 / (m_gamma - 1));
        } else {
          return threshold_l1(x, m_lambda);
        }
      }
    }
  };

  class ObjFunction {
  protected:
    int n;  // sample number
    int d;  // sample dimension
    int p;
    int m;

    vector<MatrixXd> X;
    VectorXd Y;

    vector<VectorXd> gr;
    VectorXd Xb;

    ModelParam model_param;

    double deviance;

  public:
    ObjFunction(const double *xmat, const double *y, int n, int d, int p) : n(n), d(d), p(p), m(d*p), X(d), Y(n), gr(d), Xb(n), model_param(d, p) {

      for (int i = 0; i < d; i++)
        gr[i].resize(p);
      Xb.setZero();
      for (int i = 0; i < n; i++) Y(i) = y[i];
      for (int k = 0; k < d; k++) {
        X[k].resize(n, p);
        for (int j = 0; j < p; j++) {
          for (int i = 0; i < n; i++)
            X[k](i, j) = xmat[k*p*n+j*n+i];
        }
      }
    }
    int get_dim() {
      return d;
    }
    int get_p() {
      return p;
    }
    int get_sample_num() {
      return n;
    }

    VectorXd get_grad(int idx) {
      return gr[idx];
    }

    // fabs(null fvalue - saturated fvalue)
    double get_deviance() {
      return (deviance);
    }

    double get_intercept() {
      return model_param.intercept;
    }
    VectorXd get_model_coef(int idx) {
      return model_param.beta[idx];
    }
    void set_intercept(double value) {
      model_param.intercept = value;
    }
    void set_model_coef(VectorXd value, int idx) {
      model_param.beta[idx] = value;
    }

    ModelParam get_model_param() {
      return model_param;
    }
    VectorXd get_model_Xb() const {
      return Xb;
    }

    const ModelParam &get_model_param_ref() {
      return model_param;
    }
    const VectorXd &get_model_Xb_ref() const {
      return Xb;
    }

    // reset model param and also update related aux vars
    void set_model_param(ModelParam &other_param) {
      model_param.d = other_param.d;
      model_param.beta = other_param.beta;
      model_param.intercept = other_param.intercept;
    }

    void set_model_Xb(VectorXd &other_Xb) {
      Xb = other_Xb;
    }

    // coordinate descent
    virtual VectorXd coordinate_descent(RegFunction *regfun, int idx) = 0;

    // update intercept term
    virtual void intercept_update() = 0;

    // update gradient and other aux vars
    virtual void update_auxiliary() = 0;
    virtual void update_gradient(int idx) = 0;

    // compute quadratic change of fvalue on the idx dimension
    virtual double get_local_change(const VectorXd& old, int idx) = 0;
    virtual double get_local_change_intercept(double old) = 0;

    // unpenalized function value
    virtual double eval() = 0;

    virtual ~ObjFunction(){};
    virtual double get_r2() = 0;
  };


}  // namespace picasso

#endif  // SAM_OBJECTIVE_H
