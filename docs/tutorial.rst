.. _tutorial:

********************************************************************************
Tutorial
********************************************************************************

Datastructure
=============

.. figure:: /_images/chart.png
    :figclass: figure
    :class: figure-img img-fluid


The Datastructure of the package interprets the assessment problem with four distinct objects.

The :ref:`FormDiagram <form>` element deals with the topology and geometry of the connected graph representing the projection of the path of the forces within the structure.

The :ref:`Envelope <envelope>` element deals with the geometrical data acquired from the existing structure's upper and lower surfaces (extrados and intrados).

The :ref:`Optimiser <optimiser>` element deals with the parameters for the optimisation, such as the definition of the solver, the objective function and the constraints to be applied.

The :ref:`Analysis <analysis>` element manages the information form all of the above and set up and perform the optimisation.

These elements are created in the workflow defined :ref:`here <workflow>`. To illustrate the workflow a series of :ref:`examples <examples>` are available to showcase possible applications.

Steps of the tutorial
---------------------

.. toctree::
    :maxdepth: 1
    :glob:

    tutorial/0_workflow
    tutorial/1_form
    tutorial/2_envelope
    tutorial/3_optimiser
    tutorial/4_analysis
