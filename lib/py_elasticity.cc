#include "elasticity.hh"

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>

// -------- PYTHON BINDINGS -----------------------------------------------------------------------

namespace py = pybind11;

PYBIND11_PLUGIN(dealii_example) {

  py::module m("dealii_example", "deal.II elasticity example");

  py::class_<ElasticityExample>(m, "ElasticityExample")
      .def(py::init<int>(), py::arg("refine_steps") = 4u)
      .def("solve", &ElasticityExample::solve)
      .def("h1_mat", &ElasticityExample::h1_mat, py::return_value_policy::reference_internal)
      .def("mu_mat", &ElasticityExample::mu_mat, py::return_value_policy::reference_internal)
      .def("lambda_mat", &ElasticityExample::lambda_mat, py::return_value_policy::reference_internal)
      .def("visualize", &ElasticityExample::visualize, py::arg("solution"), py::arg("filename"))
      .def("h1_0_semi_norm", &ElasticityExample::h1_0_semi_norm)
      .def("energy_norm", &ElasticityExample::energy_norm)
      .def("transfer_to", &ElasticityExample::transfer_to)
      .def("rhs", &ElasticityExample::rhs, py::return_value_policy::reference_internal);

  return m.ptr();
}
