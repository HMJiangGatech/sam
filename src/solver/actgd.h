#ifndef SAM_ACTGD_H
#define SAM_ACTGD_H

#include <cmath>
#include <string>
#include "../objective/objective.h"
#include "solver_params.h"

namespace SAM {
  class ActGDSolver {
  private:
    SolverParams m_param;
    ObjFunction *m_obj;

    std::vector<int> itercnt_path;
    std::vector<ModelParam> solution_path;

  public:
    ActGDSolver(ObjFunction *obj, SolverParams param);

    void solve();

    const std::vector<int> &get_itercnt_path() const { return itercnt_path; };
    const ModelParam &get_model_param(int i) const { return solution_path[i]; };

    ~ActGDSolver() {
      delete m_obj;
      m_obj = nullptr;
    }
  };

}  // namespace SAM
#endif
