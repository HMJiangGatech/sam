#ifndef SAM_ACTNEWTON_HPP
#define SAM_ACTNEWTON_HPP

#include <cmath>
#include <string>

#include "../objective/objective.h"
#include "solver_params.h"

namespace SAM {
  class ActNewtonSolver {
  private:
    SolverParams m_param;
    ObjFunction *m_obj;

    std::vector<int> itercnt_path;

  public:
    std::vector<ModelParam> solution_path;
    ActNewtonSolver(ObjFunction *obj, SolverParams param);

    void solve(double *sse, int *df);

    const std::vector<int> &get_itercnt_path() const { return itercnt_path; };
    const ModelParam &get_model_param(int i) const { return solution_path[i]; };

    ~ActNewtonSolver() {
      delete m_obj;
      m_obj = nullptr;
    }
  };

}  // namespace SAM

#endif
