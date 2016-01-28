import logging
from cPickle import dump, dumps
import numpy as np
import matplotlib.pyplot as plt

from pydealii_bindings import ElasticityEoc,ElasticityExample
from plotting import elasticity_error_curves

logging.basicConfig()
MARKERS = ['s', 'x', 'o', 'D', '+', '|', '*', '1', '2', '3', '4', '6', '7']

def elasticity_eoc():
    for lmbda in ( 1., 10., 75., ):
        param = {"lambda": [lmbda], "mu": [1.]}
        cpp_disc = ElasticityEoc(2, 6, param)
        eoc = cpp_disc.run()
        print('absolute eoc for param {}: {}'.format(param, [a[0] for a in eoc]))


def elasticity_calculate_errors(ref_level=9,  levels=(1, 3, 5, 7), steps=10):
    ref_disc = ElasticityExample(ref_level)

    log = logging.getLogger()
    errors = {}
    lambdas = np.linspace(1., 80., steps)
    for lmbda in lambdas:
        param = {"lambda": [lmbda], "mu": [1.]}
        try:
            ref_sol = ref_disc.solve(param)
        except RuntimeError as e:
            log.error('fine solve failed lvl {} -- lambda {}'.format(ref_level, lmbda))
            for coarse_level in levels:
                if not errors.has_key(coarse_level):
                    errors[coarse_level] = []
                errors[coarse_level].append(np.nan)
            continue
        for coarse_level in levels:
            if not errors.has_key(coarse_level):
                errors[coarse_level] = []
            coarse_disc = ElasticityExample(coarse_level)
            try:
                coarse_sol = coarse_disc.solve(param)
            except RuntimeError as e:
                log.error('coarse solve failed lvl {} -- lambda {}'.format(coarse_level, lmbda))
                errors[coarse_level].append(np.nan)
                continue
            coarse_disc.visualize([coarse_sol], ['coarse_lvl_{}_lambda_{}.vtk'.format(coarse_level, lmbda)])
            prolonged = coarse_disc.transfer_to(ref_level - coarse_level, coarse_sol)
            prolonged -= ref_sol
            errors[coarse_level].append(coarse_disc.h1_0_semi_norm(prolonged)/coarse_disc.h1_0_semi_norm(ref_sol))
    for coarse_level in errors.keys():
        plt.plot(lambdas, errors[coarse_level], MARKERS[coarse_level])

    open('elas__{}__ref-{}.dump'.format(levels, ref_level), 'wb').write(dumps((errors, lambdas)))

    elasticity_error_curves(errors, lambdas)

if __name__ == '__main__':
    elasticity_calculate_errors(levels=range(2, 8), steps=20, ref_level=9)
    # elasticity_eoc()