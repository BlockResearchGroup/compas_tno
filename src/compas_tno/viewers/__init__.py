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

"""

from __future__ import absolute_import

from compas_view2.objects import Object
from compas_view2.objects import MeshObject
from compas_tno.shapes.meshdos import MeshDos
from compas_tno.diagrams.form import FormDiagram

from .animation import *  # noqa: F401 F403
from .viewer import *  # noqa: F401 F403

Object.register(MeshDos, MeshObject)
Object.register(FormDiagram, MeshObject)


__all__ = [name for name in dir() if not name.startswith('_')]
