# This file is part of the pyMOR project (http://www.pymor.org).
# Copyright 2013-2018 pyMOR developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

from pymor.core.base import ImmutableObject
from pymor.vectorarrays.interface import VectorArray


class DealIIVisualizer(ImmutableObject):
    def __init__(self, impl):
        self.impl = impl

    def visualize(
        self,
        U,
        title=None,
        legend=None,
        separate_colorbars=True,
        filename=None,
    ):
        if isinstance(U, VectorArray):
            U = [U]
        if any(len(u) != 1 for u in U):
            raise NotImplementedError
        if not separate_colorbars:
            raise NotImplementedError
        assert title is None or filename is None
        assert legend is None or len(legend) == len(U)

        base_name = title or filename or 'out'
        if len(U) == 1 and not legend:
            filenames = [base_name]
        else:
            legend = legend or list(map(str, list(range(len(U)))))
            filenames = ['_'.join((base_name, l)) for l in legend]

        for u, n in zip(U, filenames):
            uu = u.vectors[0]
            if uu.imag_part is not None:
                self.logger.warning('Imaginary part ignored.')
            self.impl.visualize(uu.real_part.impl, n + '.vtk')
