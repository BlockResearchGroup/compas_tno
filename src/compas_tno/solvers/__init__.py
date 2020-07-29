"""
********************************************************************************
compas_tno.solvers
********************************************************************************

.. currentmodule:: compas_tno.solvers

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Analysis

"""
from __future__ import absolute_import

from .mma_numpy import *
from .solver_MMA import *
from .solver_scipy import *
from .solver_IPOPT import *
from .solver_pyOpt import *
from .solver_MATLAB import *

__all__ = [name for name in dir() if not name.startswith('_')]
