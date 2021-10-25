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

MMA
========

.. autosummary::
    :toctree: generated/

    run_optimisation_MMA

MATLAB
========

.. autosummary::
    :toctree: generated/

    run_optimisation_MATLAB

PyOpt
========

.. autosummary::
    :toctree: generated/

    run_optimisation_pyOpt

Post Processing
===============

.. autosummary::
    :toctree: generated/

    post_process_general

"""
from __future__ import absolute_import

from .mma_numpy import *  # noqa: F401 F403
from .solver_scipy import *  # noqa: F401 F403
from .solver_IPOPT import *  # noqa: F401 F403
from .solver_MATLAB import *  # noqa: F401 F403
from .post_process import *  # noqa: F401 F403
from .solver_MMA import *  # noqa: F401 F403
from .solver_pyOpt import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
