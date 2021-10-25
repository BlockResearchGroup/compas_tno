.. _shape:

********************************************************************************
Shape
********************************************************************************

.. currentmodule:: compas_tno.shapes

.. currentmodule:: compas_tno.shapes.Shape

.. highlight:: python

This tutorial provides a quick tour of the generation of ``Shape`` in :mod:`compas_tno.shapes`.

The Shape can be created in a general way using the methods ``from_pointcloud`` or ``from_meshes`` see :mod:`Shape <compas_tno.shapes.Shape>` documentation.
TNO offers a library of shapes based on commom masonry geometries which will be described herein.

Diagrams Library
================

.. figure:: ../_images/shapes.png
    :figclass: figure
    :class: figure-img img-fluid

The library is accessed with a dictionary containing infomation about the type, density, fixity etc.

Rectangular shapes
-------------------

The data dictionary for rectangular diagrams look like this:

.. code-block:: Python

    data = {
        'type': 'crossvault',
        'thk': 0.5,
        'discretisation': [20, 20],
        'xy_span': [[0.0, 10.0], [0.0, 10.0]],
        't': 0.0,
    }

The ``type`` can be changed for one of the layouts depicted in the Figure above. The parameter ``t`` define the lower-bound of nodes that have no projection in the intrados.

Circular shapes
----------------

The data dictionary for circular diagrams look like this:

.. code-block:: Python

    data = {
        'type': 'dome',
        'thk': 0.15,
        'discretisation': [20, 40],
        't' : 0.0,
        'center': [5.0, 5.0],
        'radius': 5.0,
    }

The ``type`` can be changed for one of the layouts depicted in the Figure above. The parameter ``t`` define the lower-bound of nodes that have no projection in the intrados.
