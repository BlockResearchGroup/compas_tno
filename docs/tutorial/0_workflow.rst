.. _workflow:

********************************************************************************
Workflow
********************************************************************************

This page guides you over the main step of **COMPAS TNO** workflow which is summarised in the image below:

.. figure:: ../_images/workflow.png
    :figclass: figure
    :class: figure-img img-fluid

The steps are numered herein.

1. The :ref:`FormDiagram <form>` defines the flow of forces in the structure.
2. The :ref:`Shape <shape>` object defines the geometry of the masonry to be analysed.
3. The :ref:`Optimiser <optimiser>` object stores the settings that will be necessary to perform the optimisation.
4. The :ref:`Analysis <analysis>` gathers the form diagram, shape and optimiser objects, performing preconditioning operations and running the optimisation.
5. The Solution is obtained after the optimisation
