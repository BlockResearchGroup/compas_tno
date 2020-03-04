
from numpy import abs
from numpy import argmin
from numpy import array
from numpy import float64
from numpy import dot
from numpy import hstack
from numpy import isnan
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import sqrt
from numpy import sum
from numpy import vstack
from numpy import zeros
from numpy import ones
from numpy import append
from numpy.linalg import pinv
from numpy.linalg import matrix_rank
from numpy.random import rand
from numpy.random import randint

from scipy.sparse.linalg import spsolve
from scipy.optimize import fmin_slsqp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import factorized
from compas.numerical import normrow
from compas.numerical import normalizerow

from compas.numerical import connectivity_matrix
from compas.numerical import devo_numpy
from compas.numerical import equilibrium_matrix
from compas.numerical import normrow
from compas.numerical import nonpivots
from compas.numerical.linalg import _chofactor
from compas_plotters import MeshPlotter
from compas.utilities import geometric_key
from compas.numerical.linalg import spsolve_with_known

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'update_forcediagram',
    'recalculate_qs',
]

def update_forcediagram(form,force):

    n = form.number_of_vertices()
    print('Vertices on form {0}'.format(n))
    k_i = form.key_index()
    print(len(k_i))
    xyz = zeros((n, 3))
    for key in form.vertices():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
    xy = xyz[:, :2]

    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True})]
    # edges = [[k_i[u], k_i[v]] for u, v in form.edges()]
    C	 = connectivity_matrix(edges, 'csr')
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True})]
    # q = [attr['q'] for u, v, attr in form.edges(True)]
    q = [form.edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})]
    Q = diags(q)
    uv = C.dot(xy)

    _k_i = force.key_index()
    _known = []
    for x in force.fixed():
        _known.append(_k_i[x])

    _n = force.number_of_vertices()
    _xyz = zeros((_n, 3))
    for key in force.vertices():
        i = _k_i[key]
        _xyz[i, :] = force.vertex_coordinates(key)
    _xy = _xyz[:, :2]

    _edges = force.ordered_edges(form)
    _C = connectivity_matrix(_edges, 'csr')
    _Ct = _C.transpose()

    _xy = spsolve_with_known(_Ct.dot(_C), _Ct.dot(Q).dot(uv), _xy, _known)

    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    return force

def recalculate_qs(form,force):

    k_i	 = form.key_index()
    uv_i	= form.uv_index()
    vcount  = form.number_of_vertices()
    anchors = list(form.anchors())
    fixed   = list(form.fixed())
    fixed   = set(anchors + fixed)
    fixed   = [k_i[key] for key in fixed]
    free	= list(set(range(vcount)) - set(fixed))
    edges   = [(k_i[u], k_i[v]) for u, v in form.edges()]
    xyz	 = array(form.get_vertices_attributes('xyz'), dtype='float64')
    for i in range(vcount):
        xyz[i,2] = 0.0
    C	   = connectivity_matrix(edges, 'csr')
    _scale = force.scale
    _xyz   = array(force.get_vertices_attributes('xyz'), dtype='float64')
    for i in range(_xyz.shape[1]):
        _xyz[i,2] = 0.0
    _edges = force.ordered_edges(form)
    _C	 = connectivity_matrix(_edges, 'csr')
    uvw  = C.dot(xyz)
    _uvw = _C.dot(_xyz)
    l	= normrow(uvw)
    _l   = normrow(_uvw)
    f	= _scale * _l
    q	= _l / l

    print('q from force')

    return q
