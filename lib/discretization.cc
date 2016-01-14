#include "discretization.hh"

#include <fstream>

dealii::Discretization::Discretization()
  : dof_handler_(triangulation_)
  , fe_(FE_Q<dim>(1), dim) {}

dealii::Discretization::~Discretization() { dof_handler_.clear(); }

void dealii::Discretization::setup_system() {
  dof_handler_.distribute_dofs(fe_);
  hanging_node_constraints_.clear();
  DoFTools::make_hanging_node_constraints(dof_handler_, hanging_node_constraints_);
  hanging_node_constraints_.close();
  sparsity_pattern_.reinit(dof_handler_.n_dofs(), dof_handler_.n_dofs(), dof_handler_.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_, sparsity_pattern_);

  hanging_node_constraints_.condense(sparsity_pattern_);

  sparsity_pattern_.compress();

  system_matrix_.reinit(sparsity_pattern_);

  solution_.reinit(dof_handler_.n_dofs());
  system_rhs_.reinit(dof_handler_.n_dofs());
}

void dealii::Discretization::assemble_system() {
  QGauss<dim> quadrature_formula(2);

  FEValues<dim> fe_values(fe_, quadrature_formula,
                          update_values | update_gradients | update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe_.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // As was shown in previous examples as well, we need a place where to
  // store the values of the coefficients at all the quadrature points on a
  // cell. In the present situation, we have two coefficients, lambda and
  // mu.
  std::vector<double> lambda_values(n_q_points);
  std::vector<double> mu_values(n_q_points);

  // Well, we could as well have omitted the above two arrays since we will
  // use constant coefficients for both lambda and mu, which can be declared
  // like this. They both represent functions always returning the constant
  // value 1.0. Although we could omit the respective factors in the
  // assemblage of the matrix, we use them here for purpose of
  // demonstration.
  ConstantFunction<dim> lambda(1.), mu(1.);

  // Then again, we need to have the same for the right hand side. This is
  // exactly as before in previous examples. However, we now have a
  // vector-valued right hand side, which is why the data type of the
  // <code>rhs_values</code> array is changed. We initialize it by
  // <code>n_q_points</code> elements, each of which is a
  // <code>Vector@<double@></code> with <code>dim</code> elements.
  RightHandSide<dim> right_hand_side;
  std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim));

  // Now we can begin with the loop over all cells:
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_.begin_active(), endc = dof_handler_.end();
  for (; cell != endc; ++cell) {
    cell_matrix = 0;
    cell_rhs = 0;

    fe_values.reinit(cell);

    // Next we get the values of the coefficients at the quadrature
    // points. Likewise for the right hand side:
    lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
    mu.value_list(fe_values.get_quadrature_points(), mu_values);

    right_hand_side.vector_value_list(fe_values.get_quadrature_points(), rhs_values);

    // Then assemble the entries of the local stiffness matrix and right
    // hand side vector. This follows almost one-to-one the pattern
    // described in the introduction of this example.  One of the few
    // comments in place is that we can compute the number
    // <code>comp(i)</code>, i.e. the index of the only nonzero vector
    // component of shape function <code>i</code> using the
    // <code>fe.system_to_component_index(i).first</code> function call
    // below.
    //
    // (By accessing the <code>first</code> variable of the return value
    // of the <code>system_to_component_index</code> function, you might
    // already have guessed that there is more in it. In fact, the
    // function returns a <code>std::pair@<unsigned int, unsigned
    // int@></code>, of which the first element is <code>comp(i)</code>
    // and the second is the value <code>base(i)</code> also noted in the
    // introduction, i.e.  the index of this shape function within all the
    // shape functions that are nonzero in this component,
    // i.e. <code>base(i)</code> in the diction of the introduction. This
    // is not a number that we are usually interested in, however.)
    //
    // With this knowledge, we can assemble the local matrix
    // contributions:
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int component_i = fe_.system_to_component_index(i).first;

      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
        const unsigned int component_j = fe_.system_to_component_index(j).first;

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) {
          cell_matrix(i, j) +=
              // The first term is (lambda d_i u_i, d_j v_j) + (mu d_i
              // u_j, d_j v_i).  Note that
              // <code>shape_grad(i,q_point)</code> returns the
              // gradient of the only nonzero component of the i-th
              // shape function at quadrature point q_point. The
              // component <code>comp(i)</code> of the gradient, which
              // is the derivative of this only nonzero vector
              // component of the i-th shape function with respect to
              // the comp(i)th coordinate is accessed by the appended
              // brackets.
              ((fe_values.shape_grad(i, q_point)[component_i] * fe_values.shape_grad(j, q_point)[component_j] *
                lambda_values[q_point]) +
               (fe_values.shape_grad(i, q_point)[component_j] * fe_values.shape_grad(j, q_point)[component_i] *
                mu_values[q_point]) +
               // The second term is (mu nabla u_i, nabla v_j).  We
               // need not access a specific component of the
               // gradient, since we only have to compute the scalar
               // product of the two gradients, of which an
               // overloaded version of the operator* takes care, as
               // in previous examples.
               //
               // Note that by using the ?: operator, we only do this
               // if comp(i) equals comp(j), otherwise a zero is
               // added (which will be optimized away by the
               // compiler).
               ((component_i == component_j)
                    ? (fe_values.shape_grad(i, q_point) * fe_values.shape_grad(j, q_point) * mu_values[q_point])
                    : 0)) *
              fe_values.JxW(q_point);
        }
      }
    }

    // Assembling the right hand side is also just as discussed in the
    // introduction:
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      const unsigned int component_i = fe_.system_to_component_index(i).first;

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        cell_rhs(i) += fe_values.shape_value(i, q_point) * rhs_values[q_point](component_i)*fe_values.JxW(q_point);
    }

    // The transfer from local degrees of freedom into the global matrix
    // and right hand side vector does not depend on the equation under
    // consideration, and is thus the same as in all previous
    // examples. The same holds for the elimination of hanging nodes from
    // the matrix and right hand side, once we are done with assembling
    // the entire linear system:
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        system_matrix_.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));

      system_rhs_(local_dof_indices[i]) += cell_rhs(i);
    }
  }

  hanging_node_constraints_.condense(system_matrix_);
  hanging_node_constraints_.condense(system_rhs_);

  // The interpolation of the boundary values needs a small modification:
  // since the solution function is vector-valued, so need to be the
  // boundary values. The <code>ZeroFunction</code> constructor accepts a
  // parameter that tells it that it shall represent a vector valued,
  // constant zero function with that many components. By default, this
  // parameter is equal to one, in which case the <code>ZeroFunction</code>
  // object would represent a scalar function. Since the solution vector has
  // <code>dim</code> components, we need to pass <code>dim</code> as number
  // of components to the zero function as well.
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler_, 0, ZeroFunction<dim>(dim), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values, system_matrix_, solution_, system_rhs_);
}

void dealii::Discretization::_solve() {
  SolverControl solver_control(1000, 1e-12);
  SolverCG<> cg(solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix_, 1.2);

  cg.solve(system_matrix_, solution_, system_rhs_, preconditioner);

  hanging_node_constraints_.distribute(solution_);
}

void dealii::Discretization::visualize(const VectorType& ou, std::string filename) const {
  std::ofstream output(filename);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler_);

  // As said above, we need a different name for each component of the
  // solution function. To pass one name for each component, a vector of
  // strings is used. Since the number of components is the same as the
  // number of dimensions we are working in, the following
  // <code>switch</code> statement is used.
  //
  // We note that some graphics programs have restriction as to what
  // characters are allowed in the names of variables. The library therefore
  // supports only the minimal subset of these characters that is supported
  // by all programs. Basically, these are letters, numbers, underscores,
  // and some other characters, but in particular no whitespace and
  // minus/hyphen. The library will throw an exception otherwise, at least
  // if in debug mode.
  //
  // After listing the 1d, 2d, and 3d case, it is good style to let the
  // program die if we run upon a case which we did not consider. Remember
  // that the <code>Assert</code> macro generates an exception if the
  // condition in the first parameter is not satisfied. Of course, the
  // condition <code>false</code> can never be satisfied, so the program
  // will always abort whenever it gets to the default statement:
  std::vector<std::string> solution_names;
  switch (dim) {
    case 1:
      solution_names.push_back("displacement");
      break;
    case 2:
      solution_names.push_back("x_displacement");
      solution_names.push_back("y_displacement");
      break;
    case 3:
      solution_names.push_back("x_displacement");
      solution_names.push_back("y_displacement");
      solution_names.push_back("z_displacement");
      break;
    default:
      throw ExcNotImplemented();
  }

  // After setting up the names for the different components of the solution
  // vector, we can add the solution vector to the list of data vectors
  // scheduled for output. Note that the following function takes a vector
  // of strings as second argument, whereas the one which we have used in
  // all previous examples accepted a string there. In fact, the latter
  // function is only a shortcut for the function which we call here: it
  // puts the single string that is passed to it into a vector of strings
  // with only one element and forwards that to the other function.
  data_out.add_data_vector(solution_, solution_names);
  data_out.build_patches();
  data_out.write_vtk(output);
}

dealii::Discretization::VectorType dealii::Discretization::solve(const dealii::Parameter& mu) { return VectorType(10); }

namespace py = pybind11;

py::class_<dealii::Discretization> dealii::Discretization::make_py_class(py::module& module) {
  py::class_<dealii::Discretization> disc(module, "Discretization");
  disc.def(py::init<>())
      .def("solve", &dealii::Discretization::solve)
      .def("_solve", &dealii::Discretization::_solve)
      .def("assemble_system", &dealii::Discretization::assemble_system);
  return disc;
}
