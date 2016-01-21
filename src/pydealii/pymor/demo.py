import timeit
from functools import partial

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

cpp_disc = CppDiscretization(refine_steps=7)

param = {"lambda": [1.], "mu": [1.]}
u = cpp_disc.solve(param)
parameter_type = Parameter(param).parameter_type
cpp_disc.visualize(u, "highdim_solution_cpp.vtk")

lambda_fn, mu_fn = [GenericParameterFunctional(lambda mu: mu[n], Parameter({n: [1.]}).parameter_type) for n in ['lambda', 'mu']]

ops = [DealIIMatrixOperator(getattr(cpp_disc, name)()) for name in ['lambda_mat', 'mu_mat']]
op = LincombOperator(ops, (lambda_fn, mu_fn))
rhs = VectorFunctional(ListVectorArray([DealIIVector(cpp_disc.rhs())]))
viz = PyVis(cpp_disc)
h1_op = DealIIMatrixOperator(cpp_disc.h1_mat(), "h1_0_semi")
py_disc = StationaryDiscretization(op, rhs, products={"h1_0_semi": h1_op},
                                   visualizer=viz,
                                   parameter_space=CubicParameterSpace(parameter_type, 1e-8, 1e7))

greedy_data = greedy(py_disc, reduce_stationary_affine_linear, py_disc.parameter_space.sample_uniformly(3),
                     use_estimator=False,
                     extension_algorithm=gram_schmidt_basis_extension, max_extensions=3)
rb_disc, reconstructor = greedy_data['reduced_discretization'], greedy_data['reconstructor']

for new_param in ({"lambda": [1.], "mu": [1e7]},  {"lambda": [1.], "mu": [1.]},  {"lambda": [1.], "mu": [1e-8]},
                  {"lambda": [1e-8], "mu": [1.]},  {"lambda": [1e7], "mu": [1]}, {"lambda": [1e-8], "mu": [1e-8]}):
    for disc, s in [(cpp_disc, 'cpp'), (py_disc, 'py'), (rb_disc, 'rb')]:
        with Timer(s):
            solution = disc.solve(new_param)
            try:
                disc.visualize(solution, filename='param_{}-{}.vtk'.format(new_param,s))
            except NotImplementedError:
                py_disc.visualize(reconstructor.reconstruct(solution), filename='param_{}-{}.vtk'.format(new_param,s))

            fr = disc.h1_0_semi_norm(solution)
            print('Type: {}\nParam: {}\nNorm: {}'.format(s, new_param, fr))
