# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from pymor.basic import *


# instantiate deal.II model
from dealii_elasticity import ElasticityExample
cpp_disc = ElasticityExample(refine_steps=7)


# wrap as pyMOR discretization
from pymor_dealii.pymor.operator import DealIIMatrixOperator
from pymor_dealii.pymor.vectorarray import DealIIVectorSpace
from pymor_dealii.pymor.gui import DealIIVisualizer


def run():
    d = StationaryDiscretization(
        operator=LincombOperator([DealIIMatrixOperator(cpp_disc.lambda_mat()), DealIIMatrixOperator(cpp_disc.mu_mat())],
                                [ProjectionParameterFunctional('lambda', ()), ProjectionParameterFunctional('mu', ())]),

        rhs=VectorOperator(DealIIVectorSpace.make_array([cpp_disc.rhs()])),

        products={'energy': DealIIMatrixOperator(cpp_disc.mu_mat())},

        visualizer=DealIIVisualizer(cpp_disc)
    )
    d = d.with_(parameter_space=CubicParameterSpace(d.parameter_type, 1, 10))


    # choose reduction method
    reductor = CoerciveRBReductor(
        d,
        product=d.energy_product,
        coercivity_estimator=ExpressionParameterFunctional("max(mu)", d.parameter_type)
    )


    # greedy basis generation
    greedy_data = greedy(d, reductor, d.parameter_space.sample_uniformly(3),
                        use_estimator=True,
                        extension_params={'method': 'gram_schmidt'}, max_extensions=5)


    # get reduced order model
    rd = greedy_data['rd']


    # validate reduced order model
    result = reduction_error_analysis(rd, d, reductor,
                                    test_mus=10, basis_sizes=reductor.RB.dim + 1,
                                    estimator=True, condition=True, error_norms=[d.energy_norm],
                                    plot=True)


    # visualize solution for parameter with maximum reduction error
    mu_max = result['max_error_mus'][0, -1]
    U = d.solve(mu_max)
    U_rb = reductor.reconstruct(rd.solve(mu_max))
    ERR = U - U_rb
    d.visualize([U, U_rb, ERR], legend=['fom', 'rom', 'error'])

    return result


if __name__ == '__main__':
    # print/plot results of validation
    from matplotlib import pyplot as plt
    result = run()
    print(result['summary'])
    plt.show()
