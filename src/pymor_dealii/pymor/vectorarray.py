# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

try:
    import pymor_dealii_bindings as pd2

    HAVE_DEALII = True
except ImportError:
    HAVE_DEALII = False

import numpy as np
from pymor.vectorarrays.list import ComplexifiedListVectorSpace, CopyOnWriteVector


class DealIIVector(CopyOnWriteVector):
    """Wraps a DealII vector to make it usable with ListVectorArray."""

    def __init__(self, impl):
        self.__auto_init(locals())

    @classmethod
    def from_instance(cls, instance):
        return cls(instance.impl)

    def to_numpy(self, ensure_copy=False):
        if ensure_copy:
            return np.array(self.impl, copy=True, dtype=np.double)
        else:
            self._copy_data_if_needed()
            return np.array(self.impl, copy=False, dtype=np.double)

    @property
    def dim(self):
        return self.impl.size()

    def _copy_data(self):
        self.impl = pd2.Vector(self.impl)

    def _scal(self, alpha):
        if self.dim == 0:
            return
        self.impl *= alpha

    def _axpy(self, alpha, x):
        if x.dim == 0:
            return
        if x is self:
            self.impl *= 1.0 + alpha
        elif x == 1:
            self.impl += x.impl
        elif x == -1:
            self.impl -= x.impl
        else:
            self.impl.axpy(alpha, x.impl)

    def inner(self, other):
        return 0 if self.dim == 0 else self.impl * other.impl

    def l1_norm(self):
        # dealII throws an exception on 0 length norms
        return 0 if self.dim == 0 else self.impl.l1_norm()

    def norm(self):
        return 0 if self.dim == 0 else self.impl.l2_norm()

    def norm2(self):
        return 0 if self.dim == 0 else self.impl.l2_norm() ** 2

    def sup_norm(self):
        return 0 if self.dim == 0 else self.impl.linfty_norm()

    def dofs(self, dof_indices):
        if len(dof_indices) == 0:
            return np.array([], dtype=np.intc)
        assert 0 <= np.min(dof_indices)
        assert np.max(dof_indices) < self.dim
        return np.array([self.impl[i] for i in dof_indices])

    def amax(self):
        A = np.abs(self.to_numpy())
        max_ind = np.argmax(A)
        return max_ind, A[max_ind]


class DealIIVectorSpace(ComplexifiedListVectorSpace):

    real_vector_type = DealIIVector

    def __init__(self, dim):
        self.__auto_init(locals())

    def __eq__(self, other):
        return (
            type(other) is DealIIVectorSpace
            and self.dim == other.dim
        )

    @classmethod
    def space_from_vector_obj(cls, vec):
        return cls(vec.size())

    @classmethod
    def space_from_dim(cls, dim):
        return cls(dim)

    def real_zero_vector(self):
        return DealIIVector(pd2.Vector(self.dim))

    def real_make_vector(self, obj):
        return DealIIVector(obj)

    def real_vector_from_numpy(self, data, ensure_copy=False):
        v = self.real_zero_vector()
        v.to_numpy()[:] = data
        return v
