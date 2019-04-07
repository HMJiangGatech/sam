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
#include "../objective/LinearObjective.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]


using namespace SAM;

extern "C" void grplasso(double *yy, double *XX, double *lambda, int *nnlambda, int *nn, int *dd, int *pp, double *ww, int *mmax_ite, double *tthol, char** regfunc, int *iinput, int *df, double *sse, double *func_norm)
{
  int n,d,p,max_ite,nlambda;
  int input;

  double thol;
  double lambda_max;

  nlambda = *nnlambda;
  n = *nn;
  d = *dd;
  p = *pp;
  max_ite = *mmax_ite;
  thol = *tthol;
  input = *iinput;

  lambda_max = 0;

  vector<MatrixXd> V(d);
  VectorXd y(n);

  for (int i = 0; i < n; i++)
    y(i) = yy[i];


  for (int i = 0; i < d; i++){

    MatrixXd X(n, p);
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < p; k++) {
        X(j,k) = XX[i*n*p + k*n + j];
      }
    }

    //std::cout << "X:" << std::endl << X << std::endl;

    Eigen::JacobiSVD<MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd S = svd.singularValues();

    X = svd.matrixU();

    V[i] = svd.matrixV();
    for (int j = 0; j < p; j++)
      for (int k = 0; k < p; k++)
        V[i](j,k) /= S[k];

    for (int j = 0; j < n; j++) {
      for (int k = 0; k < p; k++) {
        XX[i*n*p + k*n + j] = X(j, k);
      }
    }

    if (input == 0) {
      lambda_max = std::max(lambda_max, calc_norm(X.transpose()*y)/n);
    }
  }

  if (input == 0) {
    for (int i = 0; i < nlambda; i++)
      lambda[i] = lambda[i] * lambda_max;
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

  ObjFunction *obj = new LinearObjective(XX, yy, n, d, p, param->include_intercept);

  ActNewtonSolver solver(obj, *param);

  vector<vector<VectorXd> > beta_history;
  solver.solve(sse, df);

  assert(beta_history.size() == (unsigned int)nlambda);

  //update funcnorm, sse, ww, df
  vector<VectorXd> cur_beta(d, VectorXd::Zero(p));

  assert(solver.solution_path.size() == (unsigned int)nlambda);

  for (int i = 0; i < nlambda; i++) {

    ModelParam &model = solver.solution_path[i];

    for (int j = 0; j < d; j++) {
      func_norm[i*d+j] = calc_norm(model.beta[j]);
      VectorXd res = V[j] * model.beta[j];
      for (int k = 0; k < p; k++) {
        ww[i*d*p + j*p + k] = res(k);
      }
    }
  }

  delete param;
}
