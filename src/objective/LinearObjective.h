#ifndef LINEAROBJECTIVE_HPP
#define LINEAROBJECTIVE_HPP

#include "objective.h"

namespace SAM {

  class LinearObjective : public ObjFunction {
  private:
    VectorXd r;
    vector<MatrixXd> XX;

  public:
    LinearObjective(const double *xmat, const double *y, int n, int d, int p, bool include_intercept);
    VectorXd coordinate_descent(RegFunction *regfunc, int idx);

    void intercept_update();
    void update_auxiliary();
    void update_gradient(int idx);

    double get_local_change(const VectorXd &old, int idx);
    double get_local_change_intercept(double old);

    double eval();
    double get_r2();

  };
}

#endif
