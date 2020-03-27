"""
********************************************************************************
compas_tno.algorithms
********************************************************************************

.. currentmodule:: compas_tno.algorithms

.. autosummary::
    :toctree: generated/
    :nosignatures:

    initialise_problem
    initialise_form
    z_from_form
    zlq_from_qid
    zlq_from_q
    q_from_qid
    z_update
    reactions
    set_up_nonlinear_optimisation
    set_up_convex_optimisation
    run_optimisation_scipy
    run_optimisation_MATLAB
    run_optimisation_MMA
    find_independents
    independents_exclude
    independents_include
    inds_incl_excl
    check_independents
    check_horizontal


"""
from __future__ import absolute_import

from .constraints import *
from .derivatives import *
from .equilibrium import *
from .independents import *
from .objectives import *
from .problems import *
from .setup import *

from .solver_MMA import *
from .solver_scipy import *
from .solver_IPOPT import *
from .solver_pyOpt import *
from .solver_MATLAB import *

# Maybe add warnings to this section ?
# try:
#     from .solver_IPOPT import *
# except:
#     print('Warning: Solver IPOPT not available, calling it won\'t be possible.')
# try:
#     from .solver_pyOpt import *
# except:
#     print('Warning: Solver pyOpt not available, calling it won\'t be possible.')
# try:
#     from .solver_MATLAB import *
# except:
#     print('Warning: MATLAB engine not available, calling it won\'t be possible.')


__all__ = [name for name in dir() if not name.startswith('_')]
