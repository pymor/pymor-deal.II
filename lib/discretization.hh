#ifndef DISCRETIZATION_HH
#define DISCRETIZATION_HH

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <pybind11/pybind11.h>

#include "rhs.hh"

namespace dealii {

class Discretization {
  static constexpr size_t dim{2};

public:
  Discretization(int refine_steps);
  ~Discretization();

  typedef double Number;
  typedef std::map<std::string, std::vector<Number>> Parameter;
  typedef Vector<Number> VectorType;

  void visualize(const VectorType& solution, std::string filename) const;
  VectorType solve(const Parameter& param);

  static pybind11::class_<Discretization> make_py_class(pybind11::module& module);

  const SparseMatrix<Number>& lambda_mat() const;
  const SparseMatrix<Number>& mu_mat() const;
  const Vector<Number>& rhs() const;

private:
  void setup_system();
  void assemble_system(Parameter param);
  void _solve();

protected:
  Triangulation<dim> triangulation_;
  DoFHandler<dim> dof_handler_;

  FESystem<dim> fe_;

  SparsityPattern sparsity_pattern_;
  SparseMatrix<Number> lambda_system_matrix_, mu_system_matrix_, full_system_matrix_;

  Vector<Number> solution_;
  Vector<Number> system_rhs_;
};

#endif // DISCRETIZATION_HH

} // namespace dealii
