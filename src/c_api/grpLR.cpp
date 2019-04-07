#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "../utils.h"
#include "../solver/actnewton.h"
#include "../objective/GLMObjective.h"
#include <iostream>

using namespace SAM;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

extern "C" void grpLR(double *A, double *y, double *lambda, int *nnlambda, double *LL0, int *nn, int *dd, int *pp, double *xx, double *aa0, int *mmax_ite, double *tthol, char** regfunc, double *aalpha, double *z, int *df, double *func_norm) {

  double thol = *tthol, L0 = *LL0;
  int nlambda = *nnlambda, n = *nn, d = *dd, p = *pp, max_ite = *mmax_ite;

  vector<MatrixXd> V(d);

  for (int i = 0; i < nlambda; i++) {
    lambda[i] /= n;
  }


  SolverParams *param = new SolverParams();
  param->set_lambdas(lambda, nlambda);
  param->gamma = 3;
  if (strcmp(*regfunc, "MCP") == 0) {
    param->reg_type = MCP;
  } else if (strcmp(*regfunc, "SCAD") == 0) {
    param->reg_type = SCAD;
  } else {
    param->reg_type = L1;
  }
  param->include_intercept = true;
  param->prec = thol;
  param->max_iter = max_ite;
  param->num_relaxation_round = 10;

  ObjFunction *obj = new LogisticObjective(A, y, n, d, p, L0, param->include_intercept);

  ActNewtonSolver solver(obj, *param);

  vector<vector<VectorXd> > beta_history;
  double *sse = new double[nlambda];
  solver.solve(sse, df);

  assert(solver.solution_path.size() == (unsigned int)nlambda);
  for (int i = 0; i < nlambda; i++) {
    ModelParam &model = solver.solution_path[i];
    for (int j = 0; j < d; j++) {
      func_norm[i*d + j] = calc_norm(model.beta[j]);
      for (int k = 0; k < p; k++) {
        xx[i*(d*p+1) + j*p + k] = model.beta[j](k);
      }
    }
    xx[i*(d*p+1) + d*p] = model.intercept;
  }

  delete param;
}
