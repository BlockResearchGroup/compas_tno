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
    q_from_variables
    xyz_from_q
    xyz_from_xopt
    weights_from_xyz
    weights_from_xyz_dict
    compute_reactions
    equilibrium_residual


Independents
============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    find_independents_backward
    find_independents_forward
    find_independents_QR
    find_independents
    independents_exclude
    independents_include
    inds_incl_excl
    check_independents
    check_horizontal_loads


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
        q_from_variables,
        xyz_from_q,
        compute_reactions,
        xyz_from_xopt,
        weights_from_xyz,
        weights_from_xyz_dict,
        equilibrium_residual
    )
    from .independents import (
        find_independents_forward,
        find_independents_backward,
        find_independents_QR,
        find_independents,
        independents_exclude,
        independents_include,
        inds_incl_excl,
        check_independents,
        check_horizontal_loads
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
    'q_from_variables',
    'xyz_from_q',
    'compute_reactions',
    'xyz_from_xopt',
    'weights_from_xyz',
    'weights_from_xyz_dict',
    'equilibrium_residual',

    'find_independents_forward',
    'find_independents_backward',
    'find_independents_QR',
    'find_independents',
    'independents_exclude',
    'independents_include',
    'inds_incl_excl',
    'check_independents',
    'check_horizontal_loads',

    'constrained_smoothing',
    'apply_sag',

    'form_update_with_parallelisation',
    'force_update_from_form',
    'reciprocal_from_form'
]
