#ifndef PICASSO_SOLVER_PARAMS_H
#define PICASSO_SOLVER_PARAMS_H

#include <vector>

namespace SAM {
  enum RegType { L1, SCAD, MCP };

  // training parameters
  class SolverParams {
  public:
    /*! number of regularization parameters */
    unsigned num_lambda;

    /*! the last paramter on the regularization path */
    double target_lambda;

    /*! type of regularization terms */
    RegType reg_type;

    /*! gamma param for SCAD and MCP regularization */
    double gamma;

    /*！ rounds of relaxation when solving SCAD and MCP penalty */
    unsigned num_relaxation_round;

    /*! precision of optimization */
    double prec;

    /*! max number of iteration for innner loop */
    int max_iter;

    /*! whether or not to add intercept term */
    bool include_intercept;

    std::vector<double> lambdas;

    SolverParams();

    void set_lambdas(const double *lambda_path, int n);

    std::vector<double> get_lambda_path() const;
  };

}  // namespace SAM

#endif
