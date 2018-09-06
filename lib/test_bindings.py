# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import pytest
import numpy as np
import pydealii_bindings as dealii

def test_vector():
    v = dealii.Vector(10)
    u = dealii.Vector(10)
    ones = dealii.Vector(10)
    for i in range(len(ones)):
        ones[i] = 1
    w = dealii.Vector(u)
    assert u.size() == w.size() == v.size()
    v[1] = 3
    u[9], u[1] = 3, 3
    assert v != u
    u[9] = 0
    assert v == u
    u[1] = 0
    # currently not working
    # g = iter(u)
    # for val in u:
    #     assert val == 0

    u[1:] = dealii.Vector(9)
    u[:] = ones
    for i in range(len(u)):
        assert u[i] == 1
    v[:] = np.ones((10,), np.double)
    assert v == u

    v.axpy(1.1, u)

    ddones = dealii.Vector(100)
    ddones[:] = np.ones((100,), np.double)
    npdd = np.array(ddones, copy=False)
    assert np.allclose(npdd, np.ones((100,), dtype=np.double))
    npdd += 1.
    ddones /= 2.
    assert np.allclose(npdd, ddones)

if __name__ == "__main__":
    test_vector()
