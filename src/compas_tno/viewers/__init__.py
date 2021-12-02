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

from compas_view2.objects.object import Object
from compas_view2.objects.meshobject import MeshObject
from compas_tno.shapes.meshdos import MeshDos

from .animation import *  # noqa: F401 F403
from .viewer import *  # noqa: F401 F403

Object.register(MeshDos, MeshObject)

__all__ = [name for name in dir() if not name.startswith('_')]
