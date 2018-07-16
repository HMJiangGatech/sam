#include "utils.hpp"

namespace SAM {
  double calc_norm(VectorXd x) {
    int res = 0;
    for (int i = 0; i < x.size(); i++)
      res += x(i) * x(i);
    return sqrt(res);
  }
}
