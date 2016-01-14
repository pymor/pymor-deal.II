
from pymor.parameters.spaces import CubicParameterSpace
from pymor.vectorarrays.list import ListVectorArray
from pymor.operators.constructions import LincombOperator, VectorFunctional
from pymor.algorithms.basisextension import gram_schmidt_basis_extension
from pymor.parameters.functionals import ProjectionParameterFunctional, ExpressionParameterFunctional
from pymor.reductors.linear import reduce_stationary_affine_linear

from pydealii_bindings import Discretization as CppDiscretization
from pydealii.pymor.operator import DealIIMatrixOperator
from pydealii.pymor.vectorarray import DealIIVector

d = CppDiscretization(refine_steps=6)

param = {"lambda": [1.], "mu": [1.]}
u = d.solve(param)

d.visualize(u, "highdim_solution.vtk")

    # def parameter_functional_factory(x, y):
    #     return ProjectionParameterFunctional(component_name='stiffness',
    #                                          component_shape=(2, 1),
    #                                          coordinates=(args['YBLOCKS'] - y - 1, x),
    #                                          name='diffusion_{}_{}'.format(x, y))
    # parameter_functionals = tuple(parameter_functional_factory(x, y)
    #                               for x, y in product(xrange(args['XBLOCKS']), xrange(args['YBLOCKS'])))
ops = (DealIIMatrixOperator(getattr(d, name)()) for name in ['lambda_mat', 'mu_mat'])
# op = LincombOperator(ops)
rhs = VectorFunctional(ListVectorArray([DealIIVector(d.rhs())]))

# parameter_space = CubicParameterSpace(op.parameter_type, 0.1, 1.)
#
# from pymor.discretizations.basic import StationaryDiscretization
# d = StationaryDiscretization(op, rhs, products=None,
#                              parameter_space=parameter_space)
# reductor = reduce_stationary_affine_linear(error_product=None)
# ext = gram_schmidt_basis_extension

# print('creating analytical problem... ', end='')
# analytical_problem = example_module.Example.AnalyticalProblem(4)
# print('done')
#
# print('discretizing... ', end='')
# discretization = example_module.Example.SimpleDiscretization(analytical_problem)
# print('done')
#
# operator = discretization.get_operator()
# print('parameter type of operator is: {}'.format(operator.parameter_type().report()))
#
# functional = discretization.get_rhs()
# print('parameter type of rhs is: {}'.format(functional.parameter_type().report()))
#
# mu = example_module.Dune.Pymor.Parameter(
#         ['diffusion', 'force', 'neumann', 'dirichlet'],
#         [[1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0]])
#
# print('solving for mu = {}... '.format(mu.report()), end='')
# solution = discretization.create_vector()
# discretization.solve(solution, mu)
# print('done')
#
# name = 'solution'
# filename = 'solution.txt'
# print('writing {} to \'{}\'... '.format(name, filename), end='')
# discretization.visualize(solution, filename, name)
# print('done')