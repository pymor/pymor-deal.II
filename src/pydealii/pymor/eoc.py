import itertools
import numpy as np
import matplotlib.pyplot as plt

from pydealii_bindings import ElasticityEoc,ElasticityExample

LOW, HIGH = 18e-3, 3e-0
from pprint import pprint

for lmbda in ( 1., 10., 75., ):
    param = {"lambda": [lmbda], "mu": [1.]}
    # cpp_disc = ElasticityEoc(2, 6, param)
    # eoc = cpp_disc.run()
    # print('absolute eoc for param {}: {}'.format(param, [a[0] for a in eoc]))

ref_level = 9
coarse_level = 5
ref_disc = ElasticityExample(ref_level)

errors = []
points = np.linspace(10., 30., 7)
for lmbda in points:
    param = {"lambda": [lmbda], "mu": [1.]}
    ref_sol = ref_disc.solve(param)
    coarse_disc = ElasticityExample(coarse_level)
    coarse_sol = coarse_disc.solve(param)
    coarse_disc.visualize([coarse_sol], ["coarse_lambda_{}.vtk".format(lmbda)])
    prolonged = coarse_disc.transfer_to(ref_level - coarse_level, coarse_sol)
    prolonged -= ref_sol
    errors.append(coarse_disc.h1_0_semi_norm(prolonged)/coarse_disc.h1_0_semi_norm(ref_sol))

plt.plot(points, errors)
plt.title('coarse_lvl: {} reference_lvl: {}'.format(coarse_level, ref_level))
plt.xlabel('Lambda')
plt.ylabel('|||u_H-u_h|||')
plt.show()
