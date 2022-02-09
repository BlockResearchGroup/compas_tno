"""
********************************************************************************
compas_tno.shapes
********************************************************************************

.. currentmodule:: compas_tno.shapes

Shape
=====

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Shape

MeshDos
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    MeshDos

Shape topologies
================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    rectangular_topology


"""
from __future__ import absolute_import

from .meshdos import MeshDos
from .shape import Shape
from .rectangular_topology import rectangular_topology

import compas

if not compas.IPY:
    from .circular_arch import *  # noqa: F401 F403
    from .crossvault import *  # noqa: F401 F403
    from .dome import *  # noqa: F401 F403
    from .general import *  # noqa: F401 F403
    from .pavillionvault import *  # noqa: F401 F403
    from .pointed_arch import *  # noqa: F401 F403
    from .pointed_crossvault import *  # noqa: F401 F403
    from .shells import *  # noqa: F401 F403

# __all__ = [name for name in dir() if not name.startswith('_')]

__all__ = [
    'MeshDos',
    'Shape',

    'rectangular_topology',

]
