"""
********************************************************************************
compas_tno
********************************************************************************

.. currentmodule:: compas_tno


.. toctree::
    :maxdepth: 1

    compas_tno.algorithms
    compas_tno.analysis
    compas_tno.blender
    compas_tno.diagrams
    compas_tno.optimisers
    compas_tno.plottters
    compas_tno.rhino
    compas_tno.shapes
    compas_tno.utilities
    compas_tno.viewers

"""

from __future__ import print_function

import os

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Block Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'
__version__ = '0.1.0'


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, '../../'))
DATA = os.path.abspath(os.path.join(HOME, 'data'))
DOCS = os.path.abspath(os.path.join(HOME, 'docs'))
TEMP = os.path.abspath(os.path.join(HOME, 'temp'))


__all__ = ['HOME', 'DATA', 'DOCS', 'TEMP']
