#ifndef SAM_OBJECTIVE_H
#define SAM_OBJECTIVE_H

#include <cmath>
#include "../eigen3/Eigen/Dense"
#include <vector>

#include <ctime>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
namespace SAM {

  class ModelParam {
  public:
    int p;
    int d;
    vector<VectorXd> beta;
    double intercept;

    ModelParam(int _d, int _p):beta(d) {
      d = _d;
      p = _p;
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
    virtual VectorXd threshold(VectorXd x) = 0;
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
      // TODO
      return x;
    }
  };

  class RegL1 : public RegFunction {
  private:
    double m_lambda;
    double m_gamma;

  public:
    void set_param(double lambda, double gamma) { m_lambda = lambda; m_gamma = gamma;}
    double get_lambda() { return m_lambda; };
    double threshold(double x) { return threshold_l1(x, m_lambda); }
    VectorXd threshold(VectorXd x) {
      return threshold_l1(x, m_lambda);
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
    ObjFunction(int d, int p) : d(d), p(p), m(d*p), X(d), gr(d), model_param(d, p) {}
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
    virtual double get_local_change(VectorXd old, int idx) = 0;
    virtual double get_local_change_intercept(double old) = 0;

    // unpenalized function value
    virtual double eval() = 0;

    virtual ~ObjFunction(){};
    virtual double get_r2() = 0;
  };

  class LinearRegressionObjective : public ObjFunction {
  private:
    VectorXd r;
    vector<MatrixXd> XX;

  public:
    LinearRegressionObjective(const double *xmat, const double *y, int n, int d, int p, bool include_intercept);
    VectorXd coordinate_descent(RegFunction *regfunc, int idx);

    void intercept_update();
    void update_auxiliary();
    void update_gradient(int idx);

    double get_local_change(VectorXd old, int idx);
    double get_local_change_intercept(double old);

    double eval();
    double get_r2();


  };

}  // namespace picasso

#endif  // SAM_OBJECTIVE_H
