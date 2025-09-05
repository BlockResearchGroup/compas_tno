********************************************************************************
compas_tno.problems
********************************************************************************

.. currentmodule:: compas_tno.problems

Classes
=======

.. autosummary::
    :toctree: generated/

    Problem

Initialisation
==============

.. autosummary::
    :toctree: generated/

    initialise_form
    initialise_problem_general
    adapt_problem_to_fixed_diagram
    adapt_problem_to_sym_diagram
    adapt_problem_to_sym_and_fixed_diagram
    apply_sym_to_form


Starting Points
===============

.. autosummary::
    :toctree: generated/

    startingpoint_sag
    startingpoint_loadpath
    startingpoint_tna
    startingpoint_fdm

Set up
======

.. autosummary::
    :toctree: generated/

    set_up_general_optimisation
    set_up_convex_optimisation

Objectives
==========

.. autosummary::
    :toctree: generated/

    objective_selector
    f_min_thrust
    f_max_thrust
    f_bestfit
    f_horprojection
    f_loadpath_general
    f_complementary_energy
    f_complementary_energy_nonlinear
    f_max_section
    f_constant
    f_reduce_thk
    f_tight_crosssection


Derivatives
============

.. autosummary::
    :toctree: generated/

    d_fobj
    compute_dQ
    gradient_feasibility
    gradient_reduce_thk
    gradient_tight_crosssection
    gradient_fmin
    gradient_fmax
    gradient_bestfit
    gradient_horprojection
    gradient_complementary_energy
    gradient_complementary_energy_nonlinear
    gradient_loadpath
    gradient_max_section


Constraints
===========

.. autosummary::
    :toctree: generated/

    constr_wrapper

Gradients
=========

.. autosummary::
    :toctree: generated/

    d_fobj
    compute_dQ
    gradient_feasibility
    gradient_reduce_thk
    gradient_tight_crosssection
    gradient_fmin
    gradient_fmax
    gradient_bestfit
    gradient_horprojection
    gradient_complementary_energy
    gradient_complementary_energy_nonlinear
    gradient_loadpath
    gradient_max_section

Jacobian
========

.. autosummary::
    :toctree: generated/

    d_fconstr
    sensitivities_wrapper


Bounds Update
=============

.. autosummary::
    :toctree: generated/

    ub_lb_update
    dub_dlb_update
    b_update
    db_update

Callbacks
=========

.. autosummary::
    :toctree: generated/

    callback_save_json
    callback_create_json
    save_geometry_at_iterations
