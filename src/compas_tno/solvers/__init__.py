"""
********************************************************************************
compas_tno.solvers
********************************************************************************

.. currentmodule:: compas_tno.solvers

Scipy
========

.. autosummary::
    :toctree: generated/

    run_optimisation_scipy

IPOPT
========

.. autosummary::
    :toctree: generated/

    run_optimisation_ipopt
    Wrapper_ipopt
    Wrapper_ipopt_autodiff

MATLAB
========

.. autosummary::
    :toctree: generated/

    run_optimisation_MATLAB
    run_loadpath_from_form_MATLAB
    call_and_output_CVX_MATLAB

CVXPY
========

.. autosummary::
    :toctree: generated/

    run_optimisation_CVXPY
    run_loadpath_from_form_CVXPY
    call_and_output_CVXPY

PyOpt
========

.. autosummary::
    :toctree: generated/

    run_optimisation_pyOpt

MMA
========

.. autosummary::
    :toctree: generated/

    run_optimisation_MMA
    mma_numpy

Post Processing
===============

.. autosummary::
    :toctree: generated/

    post_process_general

"""
from __future__ import absolute_import

from .mma_numpy import mma_numpy
from .solver_scipy import (
    run_optimisation_scipy
)
from .solver_IPOPT import (
    run_optimisation_ipopt,
    Wrapper_ipopt,
    Wrapper_ipopt_autodiff
)
from .solver_MATLAB import (
    run_optimisation_MATLAB,
    run_loadpath_from_form_MATLAB,
    call_and_output_CVX_MATLAB,
)
from .solver_cvxpy import (
    run_optimisation_CVXPY,
    run_loadpath_from_form_CVXPY,
    call_and_output_CVXPY
)
from .post_process import post_process_general
from .solver_MMA import run_optimisation_MMA
from .solver_pyOpt import run_optimisation_pyOpt


__all__ = [
    'run_optimisation_scipy',

    'run_optimisation_ipopt',
    'Wrapper_ipopt',
    'Wrapper_ipopt_autodiff',

    'mma_numpy',
    'run_optimisation_pyOpt',

    'run_optimisation_MATLAB',
    'run_loadpath_from_form_MATLAB',
    'call_and_output_CVX_MATLAB',

    'run_optimisation_CVXPY',
    'run_loadpath_from_form_CVXPY',
    'call_and_output_CVXPY',

    'post_process_general',

    'run_optimisation_MMA',

]
