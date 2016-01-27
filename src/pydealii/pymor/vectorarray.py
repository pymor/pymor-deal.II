# -*- coding: utf-8 -*-
# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright Holders: Rene Milk, Stephan Rave, Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import absolute_import, division, print_function

try:
    import pydealii_bindings as pd2
    HAVE_DEALII = True
except ImportError:
    HAVE_DEALII = False

import numpy as np

from pymor.vectorarrays.interfaces import VectorSpace
from pymor.vectorarrays.list import ListVectorArray, CopyOnWriteVector


class DealIIVector(CopyOnWriteVector):
    """Wraps a DealII vector to make it usable with ListVectorArray."""

    def __init__(self, impl):
        self.impl = impl

    @classmethod
    def from_instance(cls, instance):
        return cls(instance.impl)

    def _copy_data(self):
        return type(self)(pd2.Vector(self.impl))

    @classmethod
    def make_zeros(cls, subtype):
        dim = subtype
        impl = pd2.Vector(dim)
        return cls(impl)

    @property
    def data(self):
        return np.array(self.impl, copy=False, dtype=np.double)

    @property
    def dim(self):
        return self.impl.size()

    @property
    def subtype(self):
        return self.impl.size()

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

    def sup_norm(self):
        return 0 if self.dim == 0 else self.impl.linfty_norm()

    def components(self, component_indices):
        if len(component_indices) == 0:
            return np.array([], dtype=np.intc)
        assert 0 <= np.min(component_indices)
        assert np.max(component_indices) < self.dim
        return np.array([self.impl[i] for i in component_indices])

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


def DealIIVectorSpace(dim):
    return VectorSpace(ListVectorArray, (DealIIVector, dim))
