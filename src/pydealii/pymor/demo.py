import timeit

from pymor.algorithms.greedy import greedy
from pymor.parameters.base import Parameter
from pymor.parameters.spaces import CubicParameterSpace
from pymor.vectorarrays.list import ListVectorArray
from pymor.operators.constructions import LincombOperator, VectorFunctional
from pymor.algorithms.basisextension import gram_schmidt_basis_extension
from pymor.parameters.functionals import ProjectionParameterFunctional, ExpressionParameterFunctional, \
    GenericParameterFunctional
from pymor.reductors.linear import reduce_stationary_affine_linear
from pymor.discretizations.basic import StationaryDiscretization
from pymor.tools.timing import Timer

from pydealii_bindings import Discretization as CppDiscretization
from pydealii.pymor.operator import DealIIMatrixOperator
from pydealii.pymor.vectorarray import DealIIVector


class PyVis(object):
    def __init__(self, cpp_disc):
        self._cpp_disc = cpp_disc

    def visualize(self, U, _, filename):
        self._cpp_disc.visualize(U._list[0].impl, filename)

cpp_disc = CppDiscretization(refine_steps=6)

param = {"lambda": [1.], "mu": [1.]}
u = cpp_disc.solve(param)
parameter_type = Parameter(param).parameter_type
cpp_disc.visualize(u, "highdim_solution_cpp.vtk")

lambda_fn, mu_fn = [GenericParameterFunctional(lambda mu: mu[n], Parameter({n: [1.]}).parameter_type) for n in ['lambda', 'mu']]

ops = [DealIIMatrixOperator(getattr(cpp_disc, name)()) for name in ['lambda_mat', 'mu_mat']]
op = LincombOperator(ops, (lambda_fn, mu_fn))
rhs = VectorFunctional(ListVectorArray([DealIIVector(cpp_disc.rhs())]))

py_disc = StationaryDiscretization(op, rhs, products=None,
                                   visualizer=PyVis(cpp_disc),
                                   parameter_space=CubicParameterSpace(parameter_type, 0.5, 1.5))

greedy_data = greedy(py_disc, reduce_stationary_affine_linear, py_disc.parameter_space.sample_uniformly(5),
                     use_estimator=False,
                     extension_algorithm=gram_schmidt_basis_extension, max_extensions=5)
rb_disc, reconstructor = greedy_data['reduced_discretization'], greedy_data['reconstructor']

new_param = {"lambda": [0.8], "mu": [1.2]}
for disc, n in [(cpp_disc, 'cpp'),(rb_disc, 'rb')]:
    with Timer(n):
        solution = disc.solve(new_param)
        disc.visualize(solution, 'param_{}.vtk'.format(n))
