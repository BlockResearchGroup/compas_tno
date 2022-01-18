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

    equilibrium_fdm
    vertical_equilibrium_fdm
    q_from_qid
    xyz_from_q
    xyz_from_xopt
    compute_reactions


Independents
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    find_independents
    independents_exclude
    independents_include
    inds_incl_excl
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
    from .equilibrium import (
        equilibrium_fdm,
        vertical_equilibrium_fdm,
        q_from_qid,
        xyz_from_q,
        compute_reactions,
        xyz_from_xopt
    )
    from .independents import (
        find_independents,
        independents_exclude,
        independents_include,
        inds_incl_excl,
        check_independents,
        check_horizontal
    )
    from .smoothing import (
        constrained_smoothing,
        apply_sag
    )
    from .graphstatics import (
        form_update_with_parallelisation,
        force_update_from_form,
        reciprocal_from_form
    )
    # from .equilibrium_pytorch import *  # noqa: F401 F403


__all__ = [
    'equilibrium_fdm',
    'vertical_equilibrium_fdm',
    'q_from_qid',
    'xyz_from_q',
    'compute_reactions',
    'xyz_from_xopt',

    'find_independents',
    'independents_exclude',
    'independents_include',
    'inds_incl_excl',
    'check_independents',
    'check_horizontal',

    'constrained_smoothing',
    'apply_sag',

    'form_update_with_parallelisation',
    'force_update_from_form',
    'reciprocal_from_form'
]
