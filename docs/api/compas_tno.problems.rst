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

    initialize_fdm
    initialize_tna
    initialize_loadpath
    initialise_form
    initialise_problem_general
    adapt_problem_to_fixed_diagram
    adapt_problem_to_sym_diagram
    adapt_problem_to_sym_and_fixed_diagram
    apply_sym_to_form

Set up
======

.. autosummary::
    :toctree: generated/

    set_up_general_optimisation
    set_up_convex_optimisation

Starting Points
===============

.. autosummary::
    :toctree: generated/

    initialize_loadpath
    initialize_tna

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

Proxy
=====

.. autosummary::
    :toctree: generated/

    initialize_loadpath_proxy
    run_NLP_proxy
    run_NLP_proxy2

Callbacks
=========

.. autosummary::
    :toctree: generated/

    callback_save_json
    callback_create_json
    save_geometry_at_iterations
