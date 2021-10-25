.. _tutorial:
********************************************************************************
Tutorial
********************************************************************************

Datastructure
=============

.. figure:: /_images/chart.png
    :figclass: figure
    :class: figure-img img-fluid


The datastructure of the package interprets the assessment problem with four distinct elements.

The :ref:`FormDiagram <form>` element deals with the topology and geometry of the connected graph representing the projection of the path of the forces within the structure.

The :ref:`Shape <shape>` element deals with the geometrical data acquired from the existing structure's upper and lower surfaces (extrados and intrados).

The :ref:`Optimiser <optimiser>` element deals with the parameters for the optimisation, such as the definition of the solver, the objective function and the constraints to be applied.

The :ref:`Analysis <analysis>` element manages the information form all of the above and set up and perform the optimisation.

Each of the elements is explained in detail in this tutorial. A series of :ref:`examples <examples>` are available to showcase possible applications.

Steps of the tutorial
--------

.. toctree::
    :maxdepth: 1
    :glob:

    tutorial/form
    tutorial/shape
    tutorial/optimiser
    tutorial/analysis
