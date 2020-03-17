********************************************************************************
Tutorial
********************************************************************************

Overview
============

.. figure:: /_images/chart.png
    :figclass: figure
    :class: figure-img img-fluid


The package interprets the assessment problem with four distinct elements.

The :mod:`FormDiagram` element deals with the topology and geometry of the connected graph representing the projection of the path of the forces within the structure.

The :mod:`Shape` element deals with the geometrical data acquired from the existing structure's upper and lower surfaces (extrados and intrados).

The :mod:`Optimiser` element deals with the parameters for the optimisation, such as the definition of the solver, the objective function and the constraints to be applied.

The :mod:`Analysis` element manages the information form all of the above and set up and perform the optimisation.

Now check the Examples to understand how to set up and solve a optimisation problem.
