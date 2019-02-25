// This file is part of the pyMOR project (http://www.pymor.org).
// Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)


#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename Number>
void bind_vector(pybind11::module& module) {
  typedef dealii::Vector<Number> Vector;
  typedef typename Vector::size_type size_type;
  py::class_<Vector>(module, "Vector")
      .def(py::init<>())
      .def(py::init<const Vector&>())
      .def(py::init<const size_type>())
      .def("swap", &Vector::swap)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self * py::self)
      .def(py::self += py::self)
      .def(py::self -= py::self)
      .def(py::self *= Number())
      .def(py::self /= Number())
      .def("axpy", [](Vector& self, Number a, const Vector& x) { self.add(a, x); })
      .def("norm_sqr", &Vector::norm_sqr)
      .def("mean_value", &Vector::mean_value)
      .def("norm_sqr", &Vector::norm_sqr)
      .def("l1_norm", &Vector::l1_norm)
      .def("l2_norm", &Vector::l2_norm)
      .def("lp_norm", &Vector::lp_norm, py::arg("p"))
      .def("linfty_norm", &Vector::linfty_norm)
      .def("all_zero", &Vector::all_zero)
      .def("norm_sqr", &Vector::norm_sqr)
      .def("size", &Vector::size)
      .def("__getitem__",
           [](const Vector& s, size_type i) {
             if (i >= s.size())
               throw py::index_error();
             return s[i];
           })
      .def("__setitem__",
           [](Vector& s, size_type i, Number v) {
             if (i >= s.size())
               throw py::index_error();
             s[i] = v;
           })
      /// Slicing protocol (optional)
      .def("__getitem__",
           [](const Vector& s, py::slice slice) -> Vector* {
             std::size_t start, stop, step, slicelength;
             if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
               throw py::error_already_set();
             Vector* seq = new Vector(slicelength);
             for (int i = 0; i < slicelength; ++i) {
               (*seq)[i] = s[start];
               start += step;
             }
             return seq;
           })
      .def("__setitem__",
           [](Vector& s, py::slice slice, const Vector& value) {
             std::size_t start, stop, step, slicelength;
             if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
               throw py::error_already_set();
             if ((size_t)slicelength != value.size())
               throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
             for (int i = 0; i < slicelength; ++i) {
               s[start] = value[i];
               start += step;
             }
           })
      .def("__setitem__",
           [](Vector& s, py::slice slice, py::array_t<Number> value) {
             std::size_t start, stop, step, slicelength;
             py::buffer_info info = value.request();
             if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
               throw py::error_already_set();
             if (slicelength != info.size) {
               std::stringstream ss;
               ss << "Left and right hand size of slice assignment have different sizes!";
               ss << slicelength << " vs. " << info.ndim;
               throw std::runtime_error(ss.str());
             }
             for (int i = 0; i < slicelength; ++i) {
               s[start] = *(static_cast<Number*>(info.ptr) + i);
               start += step;
             }
           })
      /// Provide buffer access
      .def_buffer([](Vector& m) -> py::buffer_info {
        return py::buffer_info(&m[0],                                /* Pointer to buffer */
                               sizeof(Number),                       /* Size of one scalar */
                               py::format_descriptor<Number>::value, /* Python struct-style format descriptor */
                               1,                                    /* Number of dimensions */
                               {
                                   m.size(),
                               }, /* Buffer dimensions */
                               {sizeof(Number)});
      })
      .def("__len__", &Vector::size)
      .def("__repr__", [](const Vector& a) {
        std::stringstream ss;
        ss << "<dealii.Vector<Number> with size '" << a.size() << "'>";
        return ss.str();
      });
}

template <typename Number>
void bind_sparse_matrix(pybind11::module& module) {
  typedef dealii::SparseMatrix<Number> Matrix;
  typedef dealii::Vector<Number> Vector;

  auto cg_solve = [](Matrix& self, Vector& solution, const Vector& rhs) {
    dealii::SolverControl solver_control(20000, 1e-12);
    dealii::SolverCG<> solver(solver_control);
    dealii::PreconditionSSOR<> preconditioner;
    preconditioner.initialize(self, 1.2);
    solver.solve(self, solution, rhs, preconditioner);

    // We have made one addition, though: since we suppress output from the
    // linear solvers, we have to print the number of iterations by hand.
    std::cout << "   " << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;
  };

  py::class_<Matrix>(module, "SparseMatrix")
      .def(py::init<>())
      .def(py::init<const dealii::SparsityPattern&>())
      .def(py::self *= Number())
      .def("n", [](const Matrix& mat) { return 0 ? mat.empty() : mat.n(); })
      .def("m", [](const Matrix& mat) { return 0 ? mat.empty() : mat.m(); })
      .def("clear", &Matrix::clear)
      .def("l1_norm", &Matrix::l1_norm)
      .def("linfty_norm", &Matrix::linfty_norm)
      .def("vmult", &Matrix::template vmult<Vector, Vector>)
      .def("Tvmult", &Matrix::template Tvmult<Vector, Vector>)
      .def("get_sparsity_pattern", &Matrix::get_sparsity_pattern, py::return_value_policy::reference)
      .def("add", (void (Matrix::*)(Number, const Matrix&)) & Matrix::template add<Number>)
      .def("copy_from", (Matrix & (Matrix::*)(const Matrix&)) & Matrix::template copy_from<Number>)
      .def("cg_solve", cg_solve);
}

PYBIND11_PLUGIN(pymor_dealii_bindings) {
  py::module m("pymor_dealii_bindings", "Python bindings for deal.II");
  py::class_<dealii::SparsityPattern>(m, "SparsityPattern");
  bind_vector<double>(m);
  bind_sparse_matrix<double>(m);
  return m.ptr();
}
