"""
********************************************************************************
compas_tno.algorithms
********************************************************************************

.. currentmodule:: compas_tno.algorithms


Equilibrium
===========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    z_from_form
    z_update
    xyz_from_q
    reactions


Independents
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    find_independents
    independents_exclude
    independents_include
    check_independents
    check_horizontal


Smoothing
=========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    constrained_smoothing
    apply_sag


Graphic Statics
===============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    form_update_with_parallelisation
    force_update_from_form
    reciprocal_from_form


"""
from __future__ import absolute_import

import compas

if not compas.IPY:
    from .equilibrium import *  # noqa: F401 F403
    # from .equilibrium_pytorch import *  # noqa: F401 F403
    from .independents import *  # noqa: F401 F403
    from .smoothing import *  # noqa: F401 F403
    from .graphstatics import *  # noqa: F401 F403


__all__ = [name for name in dir() if not name.startswith('_')]
