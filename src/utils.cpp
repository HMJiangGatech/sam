#include "utils.h"

namespace SAM {
  double calc_norm(const VectorXd &x) {
    double res = 0;
    for (int i = 0; i < x.size(); i++)
      res += x(i) * x(i);
    return sqrt(res);
  }
  double sqr(double x) {
    return x * x;
  }
}
