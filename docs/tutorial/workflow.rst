.. _workflow:

********************************************************************************
Workflow
********************************************************************************

This page guides you over the main step of ``compas_tno`` workflow described in the image below:

.. figure:: ../_images/workflow.png
    :figclass: figure
    :class: figure-img img-fluid

The steps are numered herein.

1. Create the :ref:`FormDiagram <form>` representing the forces within the structure.
2. Create the masonry :ref:`Shape <shape>` representing the geometric bounds of the structure.
3. Assign boundary conditions to the :ref:`FormDiagram <form>` representing supports, loads, envelope, bounds on axial forces, etc.
4. Compute the starting point for the optimisation (see options :ref:`here <optimiser>`).
5. Define the :ref:`Optimiser <optimiser>` settings.
6. Solve the :ref:`Analysis <analysis>` and display the results.
