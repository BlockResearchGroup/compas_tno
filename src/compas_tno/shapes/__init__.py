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

from .meshdos import *  # noqa: F401 F403
from .shape import *  # noqa: F401 F403
from .rectangular_topology import *  # noqa: F401 F403

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

__all__ = [name for name in dir() if not name.startswith('_')]
