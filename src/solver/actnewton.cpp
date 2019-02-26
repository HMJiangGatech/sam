#include <Rcpp.h>
#include <RcppEigen.h>
#include "actnewton.h"
#include "../objective/objective.h"
#include "solver_params.h"
#include "../utils.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {
  ActNewtonSolver::ActNewtonSolver(ObjFunction *obj, SolverParams param)
    : m_param(param), m_obj(obj) {
    itercnt_path.clear();
    solution_path.clear();
  }


  void ActNewtonSolver::solve(double *sse, int *df) {
    int d = m_obj->get_dim();
    int p = m_obj->get_p();
    int n = m_obj->get_sample_num();

    const std::vector<double> &lambdas = m_param.get_lambda_path();
    itercnt_path.resize(lambdas.size(), 0);

    double dev_thr = m_obj->get_deviance() * m_param.prec;

    // actset_indcat[i] == 1 if i is in the active set
    std::vector<int> actset_indcat(d, 0);
    std::vector<int> actset_indcat_master(d, 0);
    // actset_idx <- which(actset_indcat==1)
    std::vector<int> actset_idx;

    vector<VectorXd> old_coef(d);
    vector<VectorXd> grad(d);
    vector<VectorXd> grad_master(d);

    VectorXd Xb_master(n);

    double max_grad = 0;
    for (int i = 0; i < d; i++) {
      grad[i].resize(p);
      grad_master[i].resize(p);
      grad[i] = m_obj->get_grad(i);
      grad[i] = grad[i].cwiseAbs();
      max_grad = std::max(max_grad, calc_norm(grad[i]));
    }

    // model parameters on the master path
    // each master parameter is relaxed into SCAD/MCP parameter
    ModelParam model_master = m_obj->get_model_param();
    Xb_master = m_obj->get_model_Xb();

    for (int i = 0; i < d; i++) grad_master[i] = grad[i];

    std::vector<double> stage_lambdas(d, 0);
    RegFunction *regfunc = new RegL1();

    for (size_t i = 0; i < lambdas.size(); i++) {
      // start with the previous solution on the master path
      m_obj->set_model_param(model_master);
      m_obj->set_model_Xb(Xb_master);

      for (int j = 0; j < d; j++) {
        grad[j] = grad_master[j];
        actset_indcat[j] = actset_indcat_master[j];
      }

      // init the active set
      double threshold;
      if (i > 0)
        threshold = 2 * lambdas[i] - lambdas[i - 1];
      else
        threshold = 2 * lambdas[i];

      for (int j = 0; j < d; ++j) {
        stage_lambdas[j] = lambdas[i];


        if (calc_norm(grad[j]) > threshold) {
          actset_indcat[j] = 1;
        }
      }


      m_obj->update_auxiliary();
      // loop level 0: multistage convex relaxation
      int loopcnt_level_0 = 0;
      int idx;
      double old_intcpt;
      VectorXd old_beta, updated_coord, beta;
      while (loopcnt_level_0 < (int)m_param.num_relaxation_round) {
        loopcnt_level_0++;

        // loop level 1: active set update
        int loopcnt_level_1 = 0;
        bool terminate_loop_level_1 = true;
        while (loopcnt_level_1 < m_param.max_iter) {
          loopcnt_level_1++;
          terminate_loop_level_1 = true;

          old_intcpt = m_obj->get_intercept();
          for (int j = 0; j < d; j++) old_coef[j] = m_obj->get_model_coef(j);

          // initialize actset_idx
          actset_idx.clear();
          for (int j = 0; j < d; j++)
            if (actset_indcat[j]) {
              regfunc->set_param(stage_lambdas[j], 0.0);
              updated_coord = m_obj->coordinate_descent(regfunc, j);
              if (calc_norm(updated_coord) > 0) actset_idx.push_back(j);
            }


          // loop level 2: proximal newton on active set
          int loopcnt_level_2 = 0;
          bool terminate_loop_level_2 = true;
          // std::cout << "Entered loop 2" << dev_thr << std::endl;
          while (loopcnt_level_2 < m_param.max_iter) {
            loopcnt_level_2++;
            terminate_loop_level_2 = true;

            for (int k = 0; k < (int)actset_idx.size(); k++) {
              idx = actset_idx[k];

              old_beta = m_obj->get_model_coef(idx);
              regfunc->set_param(stage_lambdas[idx], 0.0);

              m_obj->coordinate_descent(regfunc, idx);

              if (m_obj->get_local_change(old_beta, idx) > dev_thr)
                terminate_loop_level_2 = false;
            }

            if (m_param.include_intercept) {
              old_intcpt = m_obj->get_intercept();
              m_obj->intercept_update();
              if (m_obj->get_local_change_intercept(old_intcpt) > dev_thr)
                terminate_loop_level_2 = false;
            }

            if (terminate_loop_level_2) break;
          }

          itercnt_path[i] += loopcnt_level_2;

          terminate_loop_level_1 = true;
          // check stopping criterion 1: fvalue change
          for (int k = 0; k < (int)actset_idx.size(); ++k) {
            idx = actset_idx[k];
            if (m_obj->get_local_change(old_coef[idx], idx) > dev_thr)
              terminate_loop_level_1 = false;
          }
          if ((m_param.include_intercept) &&
              (m_obj->get_local_change_intercept(old_intcpt) > dev_thr))
            terminate_loop_level_1 = false;

          // update p and w
          m_obj->update_auxiliary();

          if (terminate_loop_level_1) break;

          // check stopping criterion 2: active set change
          bool new_active_idx = false;
          for (int k = 0; k < d; k++)
            if (actset_indcat[k] == 0) {
              m_obj->update_gradient(k);
              grad[k] = m_obj->get_grad(k).cwiseAbs();
              if (calc_norm(grad[k]) > stage_lambdas[k]) {
                actset_indcat[k] = 1;
                new_active_idx = true;
              }
            }

          if (!new_active_idx) break;
        }

        if (loopcnt_level_0 == 1) {
          const ModelParam &model_master_ref = m_obj->get_model_param_ref();
          const Eigen::VectorXd &Xb_master_ref = m_obj->get_model_Xb_ref();

          model_master.intercept = model_master_ref.intercept;

          for (int j = 0; j < d; j++) {
            model_master.beta[j] = model_master_ref.beta[j];

            grad_master[j] = grad[j];
            actset_indcat_master[j] = actset_indcat[j];
          }

          for (int j = 0; j < n; j++) Xb_master[j] = Xb_master_ref[j];

        }

        if (m_param.reg_type == L1 || i + 1 != lambdas.size()) break;

        m_obj->update_auxiliary();

        // update stage lambda
        for (int j = 0; j < d; j++) {
          beta = m_obj->get_model_coef(j);

          if (m_param.reg_type == MCP) {
            stage_lambdas[j] = (calc_norm(beta) > lambdas[i] * m_param.gamma)
              ? 0.0
              : lambdas[i] - calc_norm(beta) / m_param.gamma;

          } else if (m_param.reg_type == SCAD)
            stage_lambdas[j] =
              (calc_norm(beta) > lambdas[i] * m_param.gamma)
              ? 0.0
              : ((calc_norm(beta) > lambdas[i])
                 ? ((lambdas[i] * m_param.gamma - calc_norm(beta)) /
                    (m_param.gamma - 1))
                 : lambdas[i]);
          else
            stage_lambdas[j] = lambdas[i];
        }
      }

      solution_path.push_back(m_obj->get_model_param());

      //update sse, df
      df[i] = (int)actset_idx.size();
      sse[i] = m_obj->get_r2();


    }

    delete regfunc;
  }

}  // namespace SAM
