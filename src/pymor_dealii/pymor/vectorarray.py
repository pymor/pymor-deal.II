# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

try:
    import pymor_dealii_bindings as pd2
    HAVE_DEALII = True
except ImportError:
    HAVE_DEALII = False

import numpy as np

from pymor.vectorarrays.interfaces import VectorSpaceInterface
from pymor.vectorarrays.list import ListVectorSpace, CopyOnWriteVector


class DealIIVector(CopyOnWriteVector):
    """Wraps a DealII vector to make it usable with ListVectorArray."""

    def __init__(self, impl):
        self.impl = impl

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
            self.impl *= 1. + alpha
        else:
            self.impl.axpy(alpha, x.impl)

    def dot(self, other):
        return 0 if self.dim == 0 else self.impl * other.impl

    def l1_norm(self):
        # dealII throws an exception on 0 length norms
        return 0 if self.dim == 0 else self.impl.l1_norm()

    def l2_norm(self):
        return 0 if self.dim == 0 else self.impl.l2_norm()

    def l2_norm2(self):
        return 0 if self.dim == 0 else self.impl.l2_norm()**2

    def sup_norm(self):
        return 0 if self.dim == 0 else self.impl.linfty_norm()

    def dofs(self, dof_indices):
        if len(dof_indices) == 0:
            return np.array([], dtype=np.intc)
        assert 0 <= np.min(dof_indices)
        assert np.max(dof_indices) < self.dim
        return np.array([self.impl[i] for i in dof_indices])

    def amax(self):
        max_ind = np.argmax(self.impl)
        return max_ind, self.impl[max_ind]

    def __add__(self, other):
        return DealIIVector(self.impl + other.impl)

    def __iadd__(self, other):
        self._copy_data_if_needed()
        self.impl += other.impl
        return self

    __radd__ = __add__

    def __sub__(self, other):
        return DealIIVector(self.impl - other.impl)

    def __isub__(self, other):
        self._copy_data_if_needed()
        self.impl -= other.impl
        return self

    def __mul__(self, other):
        return DealIIVector(self.impl * other)

    def __neg__(self):
        return DealIIVector(-self.impl)


class DealIIVectorSpace(ListVectorSpace):

    def __init__(self, dim, id_=None):
        self.dim = dim
        self.id = id_

    def __eq__(self, other):
        return type(other) is DealIIVectorSpace and self.dim == other.dim and self.id == other.id

    @classmethod
    def space_from_vector_obj(cls, vec, id_):
        return cls(vec.size(), id_)

    @classmethod
    def space_from_dim(cls, dim, id_):
        return cls(dim, id_)

    def zero_vector(self):
        return DealIIVector(pd2.Vector(self.dim))

    def make_vector(self, obj):
        return DealIIVector(obj)

    def vector_from_numpy(self, data, ensure_copy=False):
        v = self.zero_vector()
        v.to_numpy()[:] = data
        return v
