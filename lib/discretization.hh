#ifndef DISCRETIZATION_HH
#define DISCRETIZATION_HH

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
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

#include <functional>

#include "rhs.hh"

namespace dealii {

//! Wrapper intended for use w/ CG algorithms
template <class Number>
class MatrixSum : public virtual Subscriptor {

public:
  typedef Number value_type;
  typedef std::vector<const SparseMatrix<Number>*> Matrices;
  MatrixSum(Matrices&& m, std::vector<Number> weights)
    : matrices_(m) {
    assert(m.size() > 0);
    assert(m.size() == weights.size());
    sum_.reinit(matrices_[0]->get_sparsity_pattern());
    sum_.copy_from(*matrices_[0]);
    sum_ *= weights[0];
    for (size_t i = 1; i < m.size(); ++i) {
      sum_.add(weights[i], *matrices_[1]);
    }
  }

  template <class OutVector, class InVector>
  void vmult(OutVector& dst, const InVector& src) const {
    for (auto&& matrix : matrices_) {
      assert(matrix);
      matrix->vmult_add(dst, src);
    }
  }

  SparseMatrix<Number>& sum() {
    return sum_;
  }

  const SparseMatrix<Number>& sum() const{
    return sum_;
  }

private:
  const Matrices matrices_;
  SparseMatrix<Number> sum_;
};

class ElasticityExample {
  static constexpr size_t dim{2};

public:
  ElasticityExample(int refine_steps);
  ~ElasticityExample();

  typedef double Number;
  typedef std::map<std::string, std::vector<Number>> Parameter;
  typedef Vector<Number> VectorType;

  void visualize(const VectorType& solution, std::string filename) const;
  VectorType solve(const Parameter& param);

  static pybind11::class_<ElasticityExample> make_py_class(pybind11::module& module);

  const SparseMatrix<Number>& lambda_mat() const;
  const SparseMatrix<Number>& mu_mat() const;
  const SparseMatrix<Number>& h1_mat() const;
  const Vector<Number>& rhs() const;

  Number h1_0_semi_norm(const Vector<Number>& v) const;

private:
  void setup_system();
  void assemble_h1();
  void assemble_system();
  void _solve(Parameter param, VectorType &solution);

protected:
  Triangulation<dim> triangulation_;
  DoFHandler<dim> dof_handler_;

  FESystem<dim> fe_;

  SparsityPattern sparsity_pattern_;
  SparseMatrix<Number> lambda_system_matrix_, mu_system_matrix_, h1_matrix_;
  Vector<Number> system_rhs_, tmp_data_;
};

#endif // DISCRETIZATION_HH

} // namespace dealii
