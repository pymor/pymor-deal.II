# -*- coding: utf-8 -*-
# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright Holders: Rene Milk, Stephan Rave, Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from __future__ import absolute_import, division, print_function

from itertools import izip
from numbers import Number

from pymor.operators.basic import OperatorBase
from pymor.operators.constructions import ZeroOperator

from pydealii.pymor.vectorarray import DealIIVectorSpace
import pydealii_bindings as pd2


class DealIIMatrixOperator(OperatorBase):
    """Wraps a dealII matrix as an |Operator|."""

    linear = True

    def __init__(self, matrix, name=None):
        self.source = DealIIVectorSpace(matrix.m())
        self.range = DealIIVectorSpace(matrix.n())
        self.matrix = matrix
        self.name = name

    def apply(self, U, ind=None, mu=None):
        assert U in self.source
        assert U.check_ind(ind)
        vectors = U._list if ind is None else [U._list[ind]] if isinstance(ind, Number) else [U._list[i] for i in ind]
        R = self.range.zeros(len(vectors))
        for u, r in zip(vectors, R._list):
            self.matrix.vmult(r.impl, u.impl)
        return R

    def apply_adjoint(self, U, ind=None, mu=None, source_product=None, range_product=None):
        assert U in self.range
        assert U.check_ind(ind)
        assert source_product is None or source_product.source == source_product.range == self.source
        assert range_product is None or range_product.source == range_product.range == self.range
        if range_product:
            PrU = range_product.apply(U, ind=ind)._list
        else:
            PrU = U._list if ind is None else [U._list[ind]] if isinstance(ind, Number) else [U._list[i] for i in ind]
        ATPrU = self.source.zeros(len(PrU))
        for u, r in zip(PrU, ATPrU._list):
            self.matrix.Tvmult(r.impl, u.impl)
        if source_product:
            return source_product.apply_inverse(ATPrU)
        else:
            return ATPrU

    def apply_inverse(self, V, ind=None, mu=None, least_squares=False):
        assert V in self.range
        if least_squares:
            raise NotImplementedError
        vectors = V._list if ind is None else [V._list[ind]] if isinstance(ind, Number) else [V._list[i] for i in ind]
        R = self.source.zeros(len(vectors))
        for r, v in zip(R._list, vectors):
            self.matrix.cg_solve(r.impl, v.impl)
        return R

    def assemble_lincomb(self, operators, coefficients, solver_options=None, name=None):
        if not all(isinstance(op, (DealIIMatrixOperator, ZeroOperator)) for op in operators):
            return None
        assert not solver_options  # linear solver is not yet configurable

        matrix = pd2.SparseMatrix(operators[0].matrix.get_sparsity_pattern())
        matrix.copy_from(operators[0].matrix)
        matrix *= coefficients[0]
        for op, c in izip(operators[1:], coefficients[1:]):
            if isinstance(op, ZeroOperator):
                continue
            matrix.add(c, op.matrix)
        return DealIIMatrixOperator(matrix, name=name)
