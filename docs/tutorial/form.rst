.. _form:
********************************************************************************
FormDiagram
********************************************************************************

.. currentmodule:: compas_tno.diagrams

.. currentmodule:: compas_tno.diagrams.FormDiagram

.. highlight:: python

This tutorial provides a quick tour of the generation of ``FormDiagram`` in :mod:`compas_tno.diagrams`.

The diagram can be created in a general way using the methods ``from_lines`` or ``from_meshes`` see :mod:`FormDiagram <compas_tno.diagrams.FormDiagram>` documentation.
TNO offers a library of diagrams based on commom layouts which will be described herein.

Diagrams Library
================

.. figure:: ../_images/formdiagrams.png
    :figclass: figure
    :class: figure-img img-fluid

The library is accessed with a dictionary containing infomation about the type, density, fixity etc.

Rectangular diagrams
----------

The data dictionary for rectangular diagrams look like this:

.. code-block:: Python

    data = {
        'type': 'cross_fd',
        'xy_span': [[0, 10], [0, 10]],
        'discretisation': 10,
        'fix': 'corners',
    }

The ``type`` can be changed for one of the layouts depicted in the Figure above. The fixity can be modified to ``all`` or ``corners``.

Circular diagrams
----------

The data dictionary for circular diagrams look like this:

.. code-block:: Python

    data = {
        'type': 'radial_fd',
        'center': [5.0, 5.0],
        'radius': 5.0,
        'discretisation': [8, 12],
        'r_oculus': 1.25,
        'diagonal': True,
        'partial_diagonal': 'right',
    }

The ``type`` can be changed for one of the layouts depicted in the Figure above. The diagonals can have the sense modified ``right`` or ``left``.
