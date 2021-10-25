.. _optimiser:
********************************************************************************
Optimiser
********************************************************************************

.. currentmodule:: compas_tno.optimisers

.. highlight:: python

This tutorial provides a quick tour of the generation of the ``Optimiser`` object in :mod:`compas_tno.optimisers`.

The ``Optimiser`` object stores in the settings the main information about the optimisation, the main elements and its options are listed below:

Variables
=========

The sets of variables available are described in the table below. They are should be passed to ``optimiser.settings['variables']`` as a list of stings.

.. rst-class:: table table-bordered

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Name
      - Variable
    * - ``q``
      - Force Densities
    * - ``zb``
      - Height of supports
    * - ``xyb``
      - Horizontal position of supports
    * - ``t``
      - Thickness (see min_thk optimisation)
    * - ``lambd``
      - Horizontal Muliplier
    * - ``tub``
      - Additional thickness upper-bound
    * - ``tlb``
      - Additional thickness lower-bound

Objecives
=========

The objeective functions available are described in the table below. They are should be passed to ``optimiser.settings['objective']`` as a string.

.. rst-class:: table table-bordered

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Name
      - Objective
    * - ``min_thrust``
      - Minimise Horizontal Thrust
    * - ``max_thrust``
      - Maximise Horizontal Thrust
    * - ``min_thk``
      - Minimise structural Thickness
    * - ``loadpath``
      - Minimise Loadpath
    * - ``target``
      - Minimise vertical distance to target (bestfit)
    * - ``feasibility``
      - Constant objective function to test feasibility
    * - ``lambd``
      - Maximise Horizontal Multiplier
    * - ``n``
      - Maximise offset from intra/extrados (equivalent to min_thk)
    * - ``Ecomp-linear``
      - Minimise complementary energy (linear)
    * - ``Ecomp-nonlinear``
      - Minimise complementary energy (quadratic)
    * - ``max_section``
      - Minimise the increase in the structural section

Constraints
=========

The constraints available are described in the table below. They are should be passed to ``optimiser.settings['constraints']`` as a list of stings.

.. rst-class:: table table-bordered

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Name
      - Constraint
    * - ``funicular``
      - All members subjected to compressive forces
    * - ``envelope``
      - Height of the vertices within vertical envelope
    * - ``envelopexy``
      - Horizontal movement of vertices within horizontal bounds (non-fixed diagram)
    * - ``reac_bounds``
      - Bounds on the direction and magnitude of the reaction forces (see dome example)

Features
=========

The features available are described in the table below. They are should be passed to ``optimiser.settings['features']`` as a list of stings.

.. rst-class:: table table-bordered

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Name
      - Feature
    * - ``fixed``
      - Fix diagram in plan (limit 'q' variables to the independent edges)
    * - ``sym``
      - Activate symmetry in the diagram ('axis_symmetry' can be passed)
    * - ``adapted-envelope``
      - Deals with the update on the vertical envelope as nodes move horizontally

Solver Selection
=========

The solver name and library should be passed as ``optimiser.settings['solver']`` and ``optimiser.settings['library']``.

.. rst-class:: table table-bordered

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Library/Solver
      - Description
    * - ``Scipy``/``SLSQP``
      - Nonlinear solver SLSQP from library SciPy is used
    * - ``IPOPT``/``IPOPT``
      - Nonlinear solver IPOPT is used
    * - ``PyOpt``/ multiple
      - One of multiple nonlinear solvers available in PyOpt is used
    * - ``MMA``/``MMA``
      - Method of moving assymptotes (MMA) is used

Starting Point
=========

The starting point should be passed as ``optimiser.settings['starting_point']`` and will execute a preconditioning optimisation before the NLOpt.

.. rst-class:: table table-bordered

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Starting Point
      - Description
    * - ``loadpath``
      - Use the loadpath (for fixed diagram) as starting point (default)
    * - ``relax``
      - Relax the diagram using the current force densities (FDM)
    * - ``sag``
      - Relax the diagram increasing the force densities in the boundary members (FDM)
    * - ``tna``
      - Compute the horizontal graphical equilibrium and returns the thrust network
    * - ``current``
      - Use the current state (force densities and geometry) to initiate the optimisation
