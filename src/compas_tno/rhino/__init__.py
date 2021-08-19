"""
********************************************************************************
compas_tno.rhino
********************************************************************************

.. currentmodule:: compas_tno.rhino

.. autosummary::
    :toctree: generated/
    :nosignatures:


"""
from __future__ import absolute_import

import compas

from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser

from .diagramartist import DiagramArtist  # noqa: F401
from .formartist import FormArtist  # noqa: F401
from .forceartist import ForceArtist  # noqa: F401
from .shapeartist import ShapeArtist

from .diagramobject import DiagramObject  # noqa: F401
from .formobject import FormObject  # noqa: F401
from .forceobject import ForceObject  # noqa: F401
from .shapeobject import ShapeObject
from .optimiserobject import OptimiserObject  # noqa: F401

from .scene import Scene  # noqa: F401
from .settings import SettingsForm  # noqa: F401

if compas.IPY:

    from compas_rhino.objects import BaseObject
    from compas_rhino.artists import BaseArtist

    DiagramObject.register(FormDiagram, FormObject)
    DiagramArtist.register(FormDiagram, FormArtist)

    DiagramObject.register(ForceDiagram, ForceObject)
    DiagramArtist.register(ForceDiagram, ForceArtist)

    BaseObject.register(Shape, ShapeObject)
    BaseArtist.register(Shape, ShapeArtist)

    BaseObject.register(Optimiser, OptimiserObject)
    BaseArtist.register(Optimiser, DiagramArtist)

__all__ = [name for name in dir() if not name.startswith('_')]
