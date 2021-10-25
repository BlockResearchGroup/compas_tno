"""
********************************************************************************
compas_tno.problems
********************************************************************************

.. currentmodule:: compas_tno.problems

Set up
======

.. autosummary::
    :toctree: generated/

    set_up_general_optimisation
    set_up_convex_optimisation

Initialisation
==============

.. autosummary::
    :toctree: generated/

    Problem
    initialise_problem
    initialise_form
    initialise_problem_general
    adapt_problem_to_fixed_diagram
    adapt_problem_to_sym_diagram
    adapt_problem_to_sym_and_fixed_diagram

Constraints
===========

.. autosummary::
    :toctree: generated/

    constr_wrapper_general

Jacobian
========

.. autosummary::
    :toctree: generated/

    sensitivities_wrapper_general

Objectives
==========

.. autosummary::
    :toctree: generated/

    f_constant
    f_tight_crosssection
    f_reduce_thk
    f_min_thrust_general
    f_max_thrust_general
    f_bestfit_general
    f_horprojection_general
    f_loadpath_general
    f_complementary_energy
    f_complementary_energy_nonlinear
    f_max_section

Gradients
=========

.. autosummary::
    :toctree: generated/

    gradient_fmin
    gradient_fmax
    gradient_feasibility
    gradient_reduce_thk
    gradient_bestfit
    gradient_loadpath
    gradient_tight_crosssection
    gradient_fmin_general
    gradient_fmax_general
    gradient_bestfit_general
    gradient_horprojection_general
    gradient_loadpath_general
    gradient_complementary_energy
    gradient_complementary_energy_nonlinear
    gradient_max_section

"""
from __future__ import absolute_import

from .constraints import *  # noqa: F401 F403
from .derivatives import *  # noqa: F401 F403
from .jacobian import *  # noqa: F401 F403
from .objectives import *  # noqa: F401 F403
from .problems import *  # noqa: F401 F403
from .callbacks import *  # noqa: F401 F403
from .initialize import *  # noqa: F401 F403
from .setup import *  # noqa: F401 F403
from .proxy import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
