"""
********************************************************************************
compas_tno.viewers
********************************************************************************

.. currentmodule:: compas_tno.viewers

Viewer
========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Viewer


Annimation
==========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    animation_from_optimisation
    animation_from_section

"""

from __future__ import absolute_import

from compas_view2.objects import Object
from compas_view2.objects import MeshObject
from compas_tno.shapes.meshdos import MeshDos
from compas_tno.diagrams.form import FormDiagram
from compas_tno.diagrams.force import ForceDiagram

from .animation import (
    animation_from_optimisation,
    animation_from_section
)
from .viewer import Viewer

Object.register(MeshDos, MeshObject)
Object.register(FormDiagram, MeshObject)
Object.register(ForceDiagram, MeshObject)


__all__ = [
    'Viewer',
    'animation_from_optimisation',
    'animation_from_section'
]
