"""
********************************************************************************
compas_tno.diagrams
********************************************************************************

.. currentmodule:: compas_tno.diagrams


Diagrams
========

.. autosummary::
    :toctree: generated/

    FormDiagram
    ForceDiagram


Rectangular diagrams
====================

.. autosummary::
    :toctree: generated/

    create_cross_form
    create_cross_diagonal
    create_cross_with_diagonal
    create_fan_form
    create_ortho_form


Circular diagrams
=================

.. autosummary::
    :toctree: generated/

    create_circular_radial_form
    create_circular_radial_spaced_form
    create_circular_spiral_form


Linear diagrams
===============

.. autosummary::
    :toctree: generated/

    create_arch_form_diagram
    create_linear_form_diagram


"""
from __future__ import absolute_import

from .force import *  # noqa: F401 F403
from .form import *  # noqa: F401 F403
from .graph import *  # noqa: F401 F403

from .diagram_arch import *  # noqa: F401 F403
from .diagram_circular import *  # noqa: F401 F403
from .diagram_rectangular import *  # noqa: F401 F403

__all__ = [name for name in dir() if not name.startswith('_')]
