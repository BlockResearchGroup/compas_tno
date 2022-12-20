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

Circular Arch
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    arch_shape
    arch_shape_polar
    arch_ub_lb_update
    arch_dub_dlb
    arch_b_update
    arch_db

Pointed Arch
=============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    pointed_arch_shape
    pointed_arch_ub_lb_update
    pointed_arch_dub_dlb

Cross Vault
===========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    cross_vault_highfields_proxy
    cross_vault_highfields
    crossvault_ub_lb_update
    crossvault_middle_update
    crossvault_dub_dlb

Pointed Cross Vault
===================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    pointed_vault_heightfields_proxy
    pointed_vault_heightfields
    pointed_vault_middle_update
    pointed_vault_ub_lb_update
    pointed_vault_dub_dlb

Dome
====

.. autosummary::
    :toctree: generated/
    :nosignatures:

    dome_heightfields_proxy
    set_dome_heighfield
    set_dome_with_spr
    set_dome_polar_coord
    geom_dome
    dome_zt_update
    dome_ub_lb_update
    dome_dub_dlb
    dome_b_update
    dome_db
    dome_b_update_with_n
    dome_db_with_n

Pavillion Vault
===============

.. autosummary::
    :toctree: generated/
    :nosignatures:

    pavillion_vault_highfields_proxy
    pavillion_vault_highfields
    pavillionvault_ub_lb_update
    pavillionvault_dub_dlb
    pavillionvault_b_update
    pavillionvault_db

Shells
================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    domical_vault
    parabolic_shell_highfields
    parabolic_shell_middle_update
    parabolic_shell_ub_lb_update

General Meshes
================

.. autosummary::
    :toctree: generated/
    :nosignatures:

    general_ub_lb_update_with_t_middle_constant
    general_db_with_t_middle_constant
    general_ub_lb_update_with_t_middle_variable
    general_db_with_t_middle_variable
    general_ub_lb_update_with_t_intrados
    general_db_with_t_intrados
    general_ub_lb_update_with_s
    general_dub_dlb_with_s
    general_ub_lb_update_with_n
    general_dub_dlb_with_n

"""
from __future__ import absolute_import

from .meshdos import MeshDos
from .shape import Shape
from .rectangular_topology import rectangular_topology

import compas

if not compas.IPY:
    from .circular_arch import (
        arch_shape,
        arch_shape_polar,
        arch_ub_lb_update,
        arch_dub_dlb,
        arch_b_update,
        arch_db
    )
    from .crossvault import (
        cross_vault_highfields_proxy,
        cross_vault_highfields,
        crossvault_ub_lb_update,
        crossvault_middle_update,
        crossvault_dub_dlb,
    )
    from .dome import (
        dome_heightfields_proxy,
        set_dome_heighfield,
        set_dome_with_spr,
        set_dome_polar_coord,
        geom_dome,
        dome_zt_update,
        dome_ub_lb_update,
        dome_dub_dlb,
        dome_b_update,
        dome_db,
        dome_b_update_with_n,
        dome_db_with_n
    )
    from .general import (
        general_ub_lb_update_with_t_middle_constant,
        general_db_with_t_middle_constant,
        general_ub_lb_update_with_t_middle_variable,
        general_db_with_t_middle_variable,
        general_ub_lb_update_with_t_intrados,
        general_db_with_t_intrados,
        general_ub_lb_update_with_s,
        general_dub_dlb_with_s,
        general_ub_lb_update_with_n,
        general_dub_dlb_with_n
    )
    from .pavillionvault import (
        pavillion_vault_highfields_proxy,
        pavillion_vault_highfields,
        pavillionvault_ub_lb_update,
        pavillionvault_dub_dlb,
        pavillionvault_b_update,
        pavillionvault_db
    )
    from .pointed_arch import (
        pointed_arch_shape,
        pointed_arch_ub_lb_update,
        pointed_arch_dub_dlb
    )
    from .pointed_crossvault import (
        pointed_vault_heightfields_proxy,
        pointed_vault_heightfields,
        pointed_vault_middle_update,
        pointed_vault_ub_lb_update,
        pointed_vault_dub_dlb
    )
    from .shells import (
        domical_vault,
        parabolic_shell_highfields,
        parabolic_shell_middle_update,
        parabolic_shell_ub_lb_update
    )


__all__ = [
    'MeshDos',
    'Shape',

    'rectangular_topology',

    'arch_shape',
    'arch_shape_polar',
    'arch_ub_lb_update',
    'arch_dub_dlb',
    'arch_b_update',
    'arch_db',

    'cross_vault_highfields_proxy',
    'cross_vault_highfields',
    'crossvault_ub_lb_update',
    'crossvault_middle_update',
    'crossvault_dub_dlb',

    'dome_heightfields_proxy',
    'set_dome_heighfield',
    'set_dome_with_spr',
    'set_dome_polar_coord',
    'geom_dome',
    'dome_zt_update',
    'dome_ub_lb_update',
    'dome_dub_dlb',
    'dome_b_update',
    'dome_db',
    'dome_b_update_with_n',
    'dome_db_with_n',

    'general_ub_lb_update_with_t_middle_constant',
    'general_db_with_t_middle_constant',
    'general_ub_lb_update_with_t_middle_variable',
    'general_db_with_t_middle_variable',
    'general_ub_lb_update_with_t_intrados',
    'general_db_with_t_intrados',
    'general_ub_lb_update_with_s',
    'general_dub_dlb_with_s',
    'general_ub_lb_update_with_n',
    'general_dub_dlb_with_n',

    'pavillion_vault_highfields_proxy',
    'pavillion_vault_highfields',
    'pavillionvault_ub_lb_update',
    'pavillionvault_dub_dlb',
    'pavillionvault_b_update',
    'pavillionvault_db',

    'pointed_arch_shape',
    'pointed_arch_ub_lb_update',
    'pointed_arch_dub_dlb',

    'pointed_vault_heightfields_proxy',
    'pointed_vault_heightfields',
    'pointed_vault_middle_update',
    'pointed_vault_ub_lb_update',
    'pointed_vault_dub_dlb',

    'domical_vault',
    'parabolic_shell_highfields',
    'parabolic_shell_middle_update',
    'parabolic_shell_ub_lb_update'
]
