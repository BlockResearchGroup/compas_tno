"""
********************************************************************************
compas_tno.shapes
********************************************************************************

.. currentmodule:: compas_tno.shapes

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Shape

"""
from __future__ import absolute_import

from .meshdos import MeshDos
from .shape import Shape
from .rectangular_topology import *

import compas

if not compas.IPY:
    from .circular_arch import *
    from .crossvault import *
    from .dome import *
    from .general import *
    from .pavillionvault import *
    from .pointed_arch import *
    from .pointed_crossvault import *
    from .shells import *

__all__ = [name for name in dir() if not name.startswith('_')]
