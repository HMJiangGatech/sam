#ifndef SAM_ACTNEWTON_HPP
#define SAM_ACTNEWTON_HPP

#include <cmath>
#include <string>

#include "../objective/objective.hpp"
#include "../solver/solver_params.hpp"

namespace SAM {
  class ActNewtonSolver {
  private:
    SolverParams m_param;
    ObjFunction *m_obj;

    std::vector<int> itercnt_path;
    std::vector<ModelParam> solution_path;

  public:
    ActNewtonSolver(ObjFunction *obj, SolverParams param);

    void solve(double *sse, double *func_norm, vector<vector<VectorXd> > &beta_history, int *df);

    const std::vector<int> &get_itercnt_path() const { return itercnt_path; };
    const ModelParam &get_model_param(int i) const { return solution_path[i]; };

    ~ActNewtonSolver() {
      delete m_obj;
      m_obj = nullptr;
    }
  };

}  // namespace SAM

#endif
