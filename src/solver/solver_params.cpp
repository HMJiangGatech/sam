#include "solver_params.h"

namespace SAM {
  // training parameters
  SolverParams::SolverParams() {
    num_lambda = 100;
    target_lambda = 1e-6;
    reg_type = L1;
    gamma = 3.0;
    num_relaxation_round = 3;
    prec = 1e-4;
    max_iter = 1000;
    include_intercept = true;
    lambdas.clear();
  }

  void SolverParams::set_lambdas(const double *lambda_path, int n) {
    lambdas.resize(n);
    for (int i = 0; i < n; i++) lambdas[i] = lambda_path[i];
    num_lambda = lambdas.size();
    target_lambda = lambdas[num_lambda - 1];
  }

  std::vector<double> SolverParams::get_lambda_path() const {
    return lambdas;
  }

}  // namespace SAM
