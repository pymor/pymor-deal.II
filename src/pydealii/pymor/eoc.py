import itertools

from pydealii_bindings import ElasticityEoc

LOW, HIGH = 18e-3, 3e-0
from pprint import pprint

param = {"lambda": [1.], "mu": [1.]}
cpp_disc = ElasticityEoc(6, 3, param)
eoc = cpp_disc.run()
half = (HIGH - LOW) / 2.
values = itertools.product((LOW, HIGH, half), (LOW, HIGH, half))

print('absolute eoc: {}'.format([a[0] for a in eoc]))
print('relative eoc: {}'.format([a[1] for a in eoc]))