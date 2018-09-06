// This file is part of the pyMOR project (http://www.pymor.org).
// Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef EXAMPLE_HH
#define EXAMPLE_HH

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <functional>

#include "rhs.hh"

using namespace dealii;

class ElasticityExample {

public:
  static constexpr size_t dim{2};

  ElasticityExample(int refine_steps);
  ~ElasticityExample();

  typedef double Number;
  typedef std::map<std::string, std::vector<Number>> Parameter;
  typedef Vector<Number> VectorType;

  void visualize(const VectorType& solution, std::string filename) const;
  VectorType solve(const Parameter& param);

  const SparseMatrix<Number>& lambda_mat() const;
  const SparseMatrix<Number>& mu_mat() const;
  const SparseMatrix<Number>& h1_mat() const;
  const Vector<Number>& rhs() const;

  Number h1_0_semi_norm(const Vector<Number>& v) const;
  Number energy_norm(const Vector<Number>& v) const;

  VectorType transfer_to(int refine_steps, const VectorType& v);

private:
  void setup_system();
  void assemble_h1();
  void assemble_system();
  void _solve(Parameter param, VectorType& solution);

protected:
  void refine_global(int refine_steps = 1);

  Triangulation<dim> triangulation_;
  DoFHandler<dim> dof_handler_;

  FESystem<dim> fe_;

  SparsityPattern sparsity_pattern_;
  SparseMatrix<Number> lambda_system_matrix_, mu_system_matrix_, h1_matrix_;
  Vector<Number> system_rhs_, tmp_data_;
};

#endif // EXAMPLE_HH
