import itertools
from pymor.reductors.coercive import CoerciveRBReductor

from pymor.algorithms.greedy import greedy
from pymor.parameters.base import Parameter
from pymor.parameters.spaces import CubicParameterSpace
from pymor.operators.constructions import LincombOperator, VectorFunctional
from pymor.algorithms.gram_schmidt import gram_schmidt
from pymor.parameters.functionals import ProjectionParameterFunctional, ExpressionParameterFunctional
from pymor.discretizations.basic import StationaryDiscretization
from pymor.tools.timing import Timer

from dealii_elasticity import ElasticityExample
from pydealii.pymor.operator import DealIIMatrixOperator


class PyVis(object):
    def __init__(self, cpp_disc):
        self._cpp_disc = cpp_disc

    def visualize(self, U, _, filename):
        self._cpp_disc.visualize(U._list[0].impl, filename)


cpp_disc = ElasticityExample(refine_steps=7)

lambda_fn, mu_fn = [ProjectionParameterFunctional(n, ()) for n in ['lambda', 'mu']]
LOW, HIGH = 1, 10
ops = [DealIIMatrixOperator(getattr(cpp_disc, name)()) for name in ['lambda_mat', 'mu_mat']]
op = LincombOperator(ops, (lambda_fn, mu_fn))
rhs = VectorFunctional(op.range.make_array([cpp_disc.rhs()]))
viz = PyVis(cpp_disc)
h1_op = DealIIMatrixOperator(cpp_disc.h1_mat(), "h1_0_semi")
energy_op = DealIIMatrixOperator(cpp_disc.mu_mat(), "energy")
py_disc = StationaryDiscretization(op, rhs, products={"energy": energy_op},
                                   visualizer=viz,
                                   parameter_space=CubicParameterSpace(op.parameter_type, LOW, HIGH))

coercivity_estimator = ExpressionParameterFunctional("max(mu)", py_disc.parameter_type)
reductor = CoerciveRBReductor(py_disc, product=energy_op, coercivity_estimator=coercivity_estimator)

greedy_data = greedy(py_disc, reductor, py_disc.parameter_space.sample_uniformly(3),
                     use_estimator=True,
                     extension_params={'method': 'gram_schmidt'}, max_extensions=3)
rb_disc = greedy_data['rd']

half = (HIGH - LOW) / 2.
values = itertools.product((LOW, HIGH, half), (LOW, HIGH, half))
for new_param in ({"lambda": [a], "mu": [b]} for a, b in values):
    for disc, s in [(cpp_disc, 'cpp'), (py_disc, 'py'), (rb_disc, 'rb')]:
        with Timer(s):
            solution = disc.solve(new_param)
            try:
                disc.visualize(solution, filename='param_{}-{}.vtk'.format(new_param, s))
            except NotImplementedError:
                py_disc.visualize(reductor.reconstruct(solution), filename='param_{}-{}.vtk'.format(new_param, s))

            fr = disc.energy_norm(solution)
            print('Type: {}\nParam: {}\nNorm: {}'.format(s, new_param, fr))
