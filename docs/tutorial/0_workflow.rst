.. _workflow:

********************************************************************************
Workflow
********************************************************************************

This page guides you over the main step of **COMPAS TNO** workflow which is summarised in the image below:

.. figure:: ../_images/workflow.png
    :figclass: figure
    :class: figure-img img-fluid

The steps are numered herein.

1. The :class:`FormDiagram <compas_tna.diagrams.FormDiagram>` defines the flow of forces in the structure.
2. The :class:`Envelope <compas_tna.envelope.Envelope>` object defines the geometry of the masonry to be analysed.
3. The :class:`Optimiser <compas_tno.optimisers.Optimiser>` object stores the settings that will be necessary to perform the optimisation.
4. The :class:`Analysis <compas_tno.analysis.Analysis>` gathers the form diagram, shape and optimiser objects, performing preconditioning operations and running the optimisation.
5. The Solution is obtained after the optimisation
