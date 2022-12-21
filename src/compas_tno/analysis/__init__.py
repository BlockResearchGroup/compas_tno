"""
********************************************************************************
compas_tno.analysis
********************************************************************************

.. currentmodule:: compas_tno.analysis

Classes
=======

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Analysis

Routines
========

.. autosummary::
    :toctree: generated/
    :nosignatures:

    limit_analysis_GSF
    thk_minmax_GSF
    max_n_minmax_GSF

"""

from .analysis import (
    Analysis
)

from .routines import (
    limit_analysis_GSF,
    thk_minmax_GSF,
    max_n_minmax_GSF
)

__all__ = [
    'Analysis',

    'limit_analysis_GSF',
    'thk_minmax_GSF',
    'max_n_minmax_GSF'
]
