# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from pymor.operators.interface import Operator
from pymor.operators.constructions import ZeroOperator

from pymor_dealii.pymor.vectorarray import DealIIVectorSpace
import pymor_dealii_bindings as pd2


class DealIIMatrixOperator(Operator):
    """Wraps a dealII matrix as an |Operator|."""

    linear = True

    def __init__(self, matrix, name=None):
        self.source = DealIIVectorSpace(matrix.m())
        self.range = DealIIVectorSpace(matrix.n())
        self.__auto_init(locals())

    def apply(self, U, mu=None):
        assert U in self.source
        R = self.range.zeros(len(U))
        for u, r in zip(U._list, R._list):
            self.matrix.vmult(r.impl, u.impl)
        return R

    def apply_transpose(self, V, mu=None):
        assert V in self.range
        U = self.source.zeros(len(V))
        for u, r in zip(V._list, U._list):
            self.matrix.Tvmult(r.impl, u.impl)
        return U

    def apply_inverse(self, V, mu=None, least_squares=False):
        assert V in self.range
        if least_squares:
            raise NotImplementedError
        R = self.source.zeros(len(V))
        for r, v in zip(R._list, V._list):
            self.matrix.cg_solve(r.impl, v.impl)
        return R

    def _assemble_lincomb(self, operators, coefficients, identity_shift=0., solver_options=None, name=None):
        if not all(isinstance(op, (DealIIMatrixOperator)) for op in operators):
            return None
        if identity_shift != 0.:
            return None
        assert not solver_options  # linear solver is not yet configurable

        matrix = pd2.SparseMatrix(operators[0].matrix.get_sparsity_pattern())
        matrix.copy_from(operators[0].matrix)
        matrix *= coefficients[0]
        for op, c in zip(operators[1:], coefficients[1:]):
            matrix.add(c, op.matrix)
        return DealIIMatrixOperator(matrix, name=name)
