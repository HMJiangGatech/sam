#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "R.h"
#include "math.h"
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#include "solver/actgd.hpp"
#include <vector>
#include "eigen3/Eigen/SVD"
#include "eigen3/Eigen/Dense"
#include "utils.hpp"
#include "solver/actnewton.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace SAM;

extern "C" void grplasso(double *yy, double *XX, double *lambda, int *nnlambda, int *nn, int *dd, int *pp, double *ww, int *mmax_ite, double *tthol, int *iinput, int *df, double *sse, double *func_norm)
{

  int counter,n,d,p,m,max_ite,nlambda;
  int ite_ext,ite_int;
  int s;
  int input;

  int gap_ext,change_ext,back;
  double ilambda,thol,gap_int,gap_tmp0,gap_tmp1;
  double lambda_max,ols_norm;
  double w_norm;
  double *aw_norm;
  double *ols;

  nlambda = *nnlambda;
  n = *nn;
  d = *dd;
  p = *pp;
  m = d*p;
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

    Eigen::JacobiSVD<MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd S = svd.singularValues();

    Eigen::MatrixXd U = svd.matrixU();

    V[i] = svd.matrixV();


    for (int j = 0; j < n; j++) {
      for (int k = 0; k < p; k++) {
        X(j, k) *= S(k);
        XX[i*n*p + k*n + j] = X(j, k);
      }
    }

    if (input == 0) {
      lambda_max = std::max(lambda_max, calc_norm(X.transpose()*y));
    }
  }

  if (input == 0) {
    for (int i = 0; i < nlambda; i++)
      lambda[i] = lambda[i] * lambda_max;
  }

  SolverParams param;
  param.set_lambdas(lambda, nlambda);
  param.reg_type = L1; // TODO: reg_type
  param.include_intercept = true;
  param.prec = thol;
  param.max_iter = max_ite;
  param.num_relaxation_round = 3;

  ObjFunction *obj = new LinearRegressionObjective(XX, yy, n, d, p, param.include_intercept);

  ActNewtonSolver solver(obj, param);

  vector<vector<VectorXd> > beta_history;
  solver.solve(sse, func_norm, beta_history, df);

  for (int i = 0; i < nlambda; i++) {
    for (int j = 0; j < d; j++) {
      beta_history[i][j] = beta_history[i][j] * V[j].transpose();
      for (int k = 0; k < p; k++) {
        ww[i*d*p + j*p + k] = beta_history[i][j](k);
      }
    }
  }

}
