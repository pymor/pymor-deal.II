// This file is part of the pyMOR project (http://www.pymor.org).
// Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "elasticity.hh"

int main() {

  int max_refine = 5;
  /* double mu = 1.; */
  /* double lambda = 1.; */
  ElasticityExample::Parameter param = {{"lambda",
                                         {
                                             1.,
                                         }},
                                        {"mu",
                                         {
                                             1.,
                                         }}};

  std::vector<double> norms;
  std::vector<double> relative_norms;

  const auto reference_refine = max_refine + 3;
  ElasticityExample reference(reference_refine);
  const auto norm = std::bind(std::mem_fn(&ElasticityExample::h1_0_semi_norm), &reference, std::placeholders::_1);
  const auto reference_solution = reference.solve(param);
  const auto reference_norm = norm(reference_solution);

  for (int i = 1; i < max_refine + 1; ++i) {
    ElasticityExample cur_disc(i);
    const auto current_solution = cur_disc.solve(param);
    auto current_solution_refined = cur_disc.transfer_to(reference_refine - i, current_solution);
    current_solution_refined -= reference_solution;
    const auto norm_value = norm(current_solution_refined);
    norms.push_back(norm_value);
    relative_norms.push_back(norm_value / reference_norm);
  }

  std::vector<std::pair<double, double>> eoc;
  for (int i = 1; i < max_refine; ++i) {
    std::cout << std::log(norms[i] / norms[i - 1]) / std::log(1 / 2.) << ", "
              << std::log(relative_norms[i] / relative_norms[i - 1]) / std::log(1 / 2.) << std::endl;
  }
}
