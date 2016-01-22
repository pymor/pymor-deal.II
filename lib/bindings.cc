#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/base/logstream.h>

#include "discretization.hh"

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace pybind11 {
#define PYBIND11_OVERLOAD_INT_NO_RETURN(class_name, name, ...)                                                         \
  {                                                                                                                    \
    pybind11::gil_scoped_acquire gil;                                                                                  \
    pybind11::function overload = pybind11::get_overload(this, #name);                                                 \
    if (overload)                                                                                                      \
      overload.call(__VA_ARGS__);                                                                                      \
  }

#define PYBIND11_OVERLOAD_NO_RETURN(class_name, name, ...)                                                             \
  PYBIND11_OVERLOAD_INT_NO_RETURN(class_name, name, __VA_ARGS__)                                                       \
  class_name::name(__VA_ARGS__)

} // namespace pybind11

namespace dealii {
template <typename Number>
class PyVector : public Vector<Number> {
  typedef Vector<Number> BaseType;

public:
  using BaseType::Vector;

  virtual void reinit(const typename BaseType::size_type N, const bool fast = true) {
    PYBIND11_OVERLOAD_NO_RETURN(Vector<Number>, reinit, N, fast);
  }

  void axpy(Number a, const Vector<Number>& x) { BaseType::add(a, x); }

  Number* data() { return BaseType::val; }

  // Special iterator data structure for python
  struct PyVectorIterator {
    PyVectorIterator(const PyVector<Number>& vec, py::object ref)
      : vec_(vec)
      , ref_(ref) {}

    Number next() {
      if (index_ == vec_.size())
        throw py::stop_iteration();
      return vec_[index_++];
    }

    const PyVector<Number>& vec_;
    py::object ref_; // keep a reference
    typename BaseType::size_type index_ = 0;
  };

  static pybind11::class_<PyVector<Number>> make_py_class(pybind11::module& module) {
    typedef dealii::Vector<Number> NumberVector;
    typedef dealii::PyVector<Number> PyNumberVector;
    typedef typename PyNumberVector::size_type size_type;
    py::class_<PyNumberVector> cls(module, "Vector");

    cls.def(py::init<>())
        .def(py::init<const PyNumberVector&>())
        .def(py::init<const size_type>())
        .def("reinit", &PyNumberVector::reinit, py::arg("N"), py::arg("fast") = true)
        .def("swap", &NumberVector::swap)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self * py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= Number())
        .def(py::self /= Number())
        .def("axpy", &PyNumberVector::axpy)
        .def("norm_sqr", &NumberVector::norm_sqr)
        .def("mean_value", &NumberVector::mean_value)
        .def("norm_sqr", &NumberVector::norm_sqr)
        .def("l1_norm", &NumberVector::l1_norm)
        .def("l2_norm", &NumberVector::l2_norm)
        .def("lp_norm", &NumberVector::lp_norm, py::arg("p"))
        .def("linfty_norm", &NumberVector::linfty_norm)
        .def("all_zero", &NumberVector::all_zero)
        .def("norm_sqr", &NumberVector::norm_sqr)
        .def("size", &NumberVector::size)
        .def("__getitem__",
             [](const NumberVector& s, size_type i) {
               if (i >= s.size())
                 throw py::index_error();
               return s[i];
             })
        .def("__setitem__",
             [](NumberVector& s, size_type i, Number v) {
               if (i >= s.size())
                 throw py::index_error();
               s[i] = v;
             })
        /// Slicing protocol (optional)
        .def("__getitem__",
             [](const NumberVector& s, py::slice slice) -> NumberVector* {
               py::ssize_t start, stop, step, slicelength;
               if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
                 throw py::error_already_set();
               NumberVector* seq = new NumberVector(slicelength);
               for (int i = 0; i < slicelength; ++i) {
                 (*seq)[i] = s[start];
                 start += step;
               }
               return seq;
             })
        .def("__setitem__",
             [](NumberVector& s, py::slice slice, const NumberVector& value) {
               py::ssize_t start, stop, step, slicelength;
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
             [](NumberVector& s, py::slice slice, py::array_t<Number> value) {
               py::ssize_t start, stop, step, slicelength;
               py::buffer_info info = value.request();
               if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
                 throw py::error_already_set();
               if (slicelength != info.count) {
                 std::stringstream ss;
                 ss << "Left and right hand size of slice assignment have different sizes!";
                 ss << slicelength << " vs. " << info.ndim;
                 throw std::runtime_error(ss.str());
               }
               for (int i = 0; i < slicelength; ++i) {
                 s[start] = *static_cast<double*>(info.ptr + i * sizeof(Number));
                 start += step;
               }
             })
        /// Provide buffer access
        .def_buffer([](PyNumberVector& m) -> py::buffer_info {
          return py::buffer_info(m.data(),                               /* Pointer to buffer */
                                 sizeof(Number),                         /* Size of one scalar */
                                 py::format_descriptor<Number>::value(), /* Python struct-style format descriptor */
                                 1,                                      /* Number of dimensions */
                                 {
                                     m.size(),
                                 }, /* Buffer dimensions */
                                 {sizeof(Number)});
        })
        .def("__len__", &NumberVector::size)
        //      .def("__iter__", [](py::object s) { return typename PyNumberVector::PyVectorIterator(s.cast<const
        //      PyNumberVector &>(), s); })
        .def("__repr__", [](const NumberVector& a) {
          std::stringstream ss;
          ss << "<dealii.Vector<Number> with size '" << a.size() << "'>";
          return ss.str();
        });
    return cls;
  }
};

template <typename Number>
class PySparseMatrix : public SparseMatrix<Number> {
  typedef SparseMatrix<Number> BaseType;

public:
  using BaseType::SparseMatrix;

  virtual void reinit(const SparsityPattern& pattern) {
    PYBIND11_OVERLOAD_NO_RETURN(SparseMatrix<Number>, reinit, pattern);
  }

  Number* data() { return BaseType::val; }

  void cg_solve(Vector<Number>& solution, const Vector<Number>& rhs) {
    SolverControl solver_control(1000, 1e-12);
    SolverCG<> solver(solver_control);
    solver.solve(*this, solution, rhs, PreconditionIdentity());

    // We have made one addition, though: since we suppress output from the
    // linear solvers, we have to print the number of iterations by hand.
    std::cout << "   " << solver_control.last_step() << " CG iterations needed to obtain convergence." << std::endl;
  }

  static pybind11::class_<PySparseMatrix<Number>> make_py_class(pybind11::module& module) {
    typedef dealii::PySparseMatrix<double> PyDoubleMatrix;
    typedef dealii::SparseMatrix<double> DoubleMatrix;
    typedef dealii::Vector<Number> NumberVector;
    return py::class_<PyDoubleMatrix>(module, "SparseMatrix")
        .alias<DoubleMatrix>()
        .def(py::init<>())
        .def(py::init<const SparsityPattern&>())
        .def(py::self *= Number())
        .def("n", [](const DoubleMatrix& mat) { return 0 ? mat.empty() : mat.n(); })
        .def("m", [](const DoubleMatrix& mat) { return 0 ? mat.empty() : mat.m(); })
        .def("clear", &DoubleMatrix::clear)
        .def("l1_norm", &DoubleMatrix::l1_norm)
        .def("linfty_norm", &DoubleMatrix::linfty_norm)
        .def("vmult", &DoubleMatrix::vmult<NumberVector, NumberVector>)
        .def("Tvmult", &DoubleMatrix::Tvmult<NumberVector, NumberVector>)
        .def("get_sparsity_pattern", &DoubleMatrix::get_sparsity_pattern, py::return_value_policy::reference)
        .def("add", (void (DoubleMatrix::*)(Number, const DoubleMatrix&)) & DoubleMatrix::add<Number>)
        .def("copy_from", (DoubleMatrix & (DoubleMatrix::*)(const DoubleMatrix&)) & DoubleMatrix::copy_from<Number>)
        .def("cg_solve", &PyDoubleMatrix::cg_solve);
  }
};

} // namespace dealii

PYBIND11_PLUGIN(pydealii_bindings) {
  typedef double Number;
  typedef dealii::PyVector<Number> PyNumberVector;
  typedef dealii::PySparseMatrix<double> PyDoubleMatrix;

  py::module m("pydealii_bindings", "pybind11 elliptic plugin");
  py::class_<dealii::SparsityPattern> sp(m, "SparsityPattern");
  auto vec = PyNumberVector::make_py_class(m);
  // moving this to fac methods does not compile
  vec.alias<dealii::Vector<Number>>();
  auto mat = PyDoubleMatrix::make_py_class(m);
  auto disc = dealii::ElasticityExample::make_py_class(m);

  return m.ptr();
}
