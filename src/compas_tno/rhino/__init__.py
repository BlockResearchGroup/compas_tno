"""
********************************************************************************
compas_tno.rhino
********************************************************************************

.. currentmodule:: compas_tno.rhino

.. autosummary::
    :toctree: generated/
    :nosignatures:

    constr_wrapper


"""
from __future__ import absolute_import

from compas_tno.diagrams import FormDiagram
from compas_tno.diagrams import ForceDiagram

from .diagramartist import DiagramArtist  # noqa: F401
from .formartist import FormArtist  # noqa: F401
from .forceartist import ForceArtist  # noqa: F401

from .diagramobject import DiagramObject  # noqa: F401
from .formobject import FormObject  # noqa: F401
from .forceobject import ForceObject  # noqa: F401

DiagramArtist.register(FormDiagram, FormArtist)
DiagramArtist.register(ForceDiagram, ForceArtist)

DiagramObject.register(FormDiagram, FormObject)
DiagramObject.register(ForceDiagram, ForceObject)

from .scene import Scene  # noqa: F401
from .settings import SettingsForm  # noqa: F401

from .shapeartist import *  # WIP

__all__ = [name for name in dir() if not name.startswith('_')]
