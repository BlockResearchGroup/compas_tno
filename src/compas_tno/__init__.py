"""
********************************************************************************
compas_tno
********************************************************************************

.. currentmodule:: compas_tno


.. toctree::
    :maxdepth: 1

    compas_tno.algorithms
    compas_tno.analysis
    compas_tno.diagrams
    compas_tno.optimisers
    compas_tno.shapes
    compas_tno.plotters
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



def get(filename):
    """Get the full path to one of the sample data files.
    Parameters
    ----------
    filename : str
        The name of the data file.
        The following are available.
        *
    Returns
    -------
    str
        The full path to the specified file.
    Notes
    -----
    The file name should be specified relative to the ``compas_tno`` sample data folder.
    This folder is only locally available if you installed ``compas_tno`` from source,
    or if you are working directly with the source.
    In all other cases, the function will get the corresponding files direcly from
    the GitHub repo, at https://raw.githubusercontent.com/BlockResearchGroup/compas_tno/master/data.
    Examples
    --------
    The ``compas_tno.get`` function is meant to be used in combination with the static
    constructors of the data structures.
    .. code-block:: python
        import compas_tno
        from compas_tno.diagrams import FormDiagram
        form = FormDiagram.from_obj(compas.get('faces.obj'))
    """
    filename = filename.strip('/')

    localpath = os.path.abspath(os.path.join(DATA, filename))

    if os.path.exists(localpath):
        return localpath
    else:
        return "https://raw.githubusercontent.com/BlockResearchGroup/compas_tno/master/data/{}".format(filename)