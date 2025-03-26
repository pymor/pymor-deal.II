# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

import pymor_dealii_bindings as pd2
from pymor.operators.list import LinearComplexifiedListVectorArrayOperatorBase

from pymor_dealii.pymor.vectorarray import DealIIVectorSpace


class DealIIMatrixOperator(LinearComplexifiedListVectorArrayOperatorBase):
    """Wraps a dealII matrix as an |Operator|."""

    def __init__(self, matrix, name=None):
        self.source = DealIIVectorSpace(matrix.m())
        self.range = DealIIVectorSpace(matrix.n())
        self.__auto_init(locals())

    def _real_apply_one_vector(self, u, mu=None, prepare_data=None):
        r = self.range.real_zero_vector()
        self.matrix.vmult(r.impl, u.impl)
        return r

    def _real_apply_inverse_one_vector(
        self, v, mu=None, initial_guess=None, least_squares=False, prepare_data=None
    ):
        if least_squares:
            raise NotImplementedError
        r = self.source.real_zero_vector()
        self.matrix.cg_solve(r.impl, v.impl)
        return r

    def _real_apply_adjoint_one_vector(self, v, mu=None, prepare_data=None):
        r = self.source.real_zero_vector()
        self.matrix.Tvmult(r.impl, v.impl)
        return r

    def _assemble_lincomb(
        self,
        operators,
        coefficients,
        identity_shift=0.0,
        solver_options=None,
        name=None,
    ):
        if not all(isinstance(op, (DealIIMatrixOperator)) for op in operators):
            return None
        if identity_shift != 0.0:
            return None
        assert not solver_options  # linear solver is not yet configurable

        matrix = pd2.SparseMatrix(operators[0].matrix.get_sparsity_pattern())
        matrix.copy_from(operators[0].matrix)
        matrix *= coefficients[0]
        for op, c in zip(operators[1:], coefficients[1:]):
            matrix.add(c, op.matrix)
        return DealIIMatrixOperator(matrix, name=name)
