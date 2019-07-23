
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
from compas_tna.utilities import apply_bounds
from compas_tna.utilities import parallelise_sparse
from compas_tna.utilities import parallelise_nodal

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax
from compas_tna.equilibrium import vertical_from_q
from compas_tna.utilities import rot90
from compas.geometry import angle_vectors_xy

from compas.numerical import connectivity_matrix
from compas.numerical import devo_numpy
from compas.numerical import equilibrium_matrix
from compas.numerical import normrow
from compas.numerical import nonpivots
from compas.numerical.linalg import _chofactor
from compas_plotters import MeshPlotter
from compas.utilities import geometric_key
from copy import deepcopy
from compas.geometry.distance import distance_point_point_xy


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'zlq_from_qid',
    'z_from_form',
    'horizontal_check',
    'update_tna',
    'scale_fdm',
    'scale_form',
    'update_form',
    'paralelise_form',
    'reactions'
]

def z_from_form(mesh):

    from compas.numerical import fd_numpy

    # preprocess

    k_i   = mesh.key_index()
    xyz   = mesh.get_vertices_attributes(('x', 'y', 'z'))
    loads = mesh.get_vertices_attributes(('px', 'py', 'pz'))
    q     = mesh.get_edges_attribute('q')
    fixed = mesh.vertices_where({'is_fixed': True})
    fixed = [k_i[k] for k in fixed]
    edges = [(k_i[u], k_i[v]) for u, v in mesh.edges()]

    # compute equilibrium
    # update the mesh geometry

    xyz, q, f, l, r = fd_numpy(xyz, edges, fixed, q, loads)

    for key, attr in mesh.vertices(True):
        index = k_i[key]
        attr['x'] = xyz[index, 0]
        attr['y'] = xyz[index, 1]
        attr['z'] = xyz[index, 2]

    return mesh

def zlq_from_qid(qid, args):

    q, ind, dep, Edinv, Ei, C, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y = args
    # q, ind, dep, Edinv, Ei, C, Ci, Cit, U, V, p, px, py, pz, tol, z, free, planar, lh, sym, *_ = args
    # q, ind, dep, Edinv, Ei, C, Ci, Cit, p, pz, z, free, planar, lh2, sym = args[:-5]
    q[ind, 0] = qid
    q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
    q_ = 1 * q
    q[sym] *= 0

    if not planar:
        z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz)
    l2 = lh + C.dot(z)**2

    return z, l2, q, q_

def horizontal_check(form, plot = False):

    # Mapping

    k_i  = form.key_index()
    uv_i = form.uv_index()

    # Vertices and edges

    n     = form.number_of_vertices()
    fixed = [k_i[key] for key in form.fixed()]
    rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
    free  = list(set(range(n)) - set(fixed) - set(rol))


    # Co-ordinates and loads

    xyz = zeros((n, 3))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)

    xy = xyz[:, :2]
    px = px[free]
    py = py[free]
    pz = pz[free]

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cf  = C[:, fixed]
    Cit = Ci.transpose()
    uvw = C.dot(xyz)
    U   = uvw[:, 0]
    V   = uvw[:, 1]
    q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

    # Horizontal checks

    Rx = Cit.dot(U * q.ravel()) - px.ravel()
    Ry = Cit.dot(V * q.ravel()) - py.ravel()
    R  = sqrt(Rx**2 + Ry**2)
    Rmax = max(R)

    if plot:
        eq_node = {key: R[k_i[key]] for key in form.vertices_where({'is_fixed': False})}
        plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
        plotter.draw_vertices(text=eq_node)
        plotter.draw_edges()
        plotter.show()

    return Rmax

def update_tna(form, delete_face=True, plots=True, save=False):

    if delete_face:
        form.delete_face(0)

    corners = list(form.vertices_where({'is_fixed': True}))
    form.set_vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
    form.update_boundaries(feet=2)

    # for key in form.edges_where({'is_external': True}):
    #     form.set_edge_attribute(key,'q',value=x_reaction)
    #     form.set_edge_attribute(key,'fmin',value=x_reaction)
    #     form.set_edge_attribute(key,'fmax',value=x_reaction)
    #     form.set_edge_attribute(key,'lmin',value=1.00)
    #     form.set_edge_attribute(key,'lmax',value=1.00)

    for u, v in form.edges_where({'is_external': False}):
        qi = form.get_edge_attribute((u,v),'q')
        a = form.vertex_coordinates(u)
        b = form.vertex_coordinates(v)
        lh = distance_point_point_xy(a,b)
        form.set_edge_attribute((u,v),'fmin',value=qi*lh)
        form.set_edge_attribute((u,v),'fmax',value=qi*lh)
        form.set_edge_attribute((u,v),'lmin',value=lh)
        form.set_edge_attribute((u,v),'lmax',value=lh)

    force = ForceDiagram.from_formdiagram(form)

    horizontal(form,force,display=False)
    # Vertical?

    if plots:
        plot_force(force, radius=0.05).show()
        plot_form(form, radius=0.05).show()

    st = 'discretize/02_complete_1div_complete'

    if save:
        force.to_obj('/Users/mricardo/compas_dev/compas_loadpath/data/'+ st +'_force.obj')
        form.to_json('/Users/mricardo/compas_dev/compas_loadpath/data/'+ st +'_tna.json')


    return form, force

def scale_fdm(form, r):

    uv_i = form.uv_index()
    q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True, 'is_external' : False})])[:, newaxis]
    q = q * r

    for u, v in form.edges_where({'is_external': False}):
        if form.get_edge_attribute((u,v),'is_edge') is True:
            i = uv_i[(u,v)]
            [qi] = q[i]
            form.set_edge_attribute((u,v),'q',value=qi)

    form = z_from_form(form)

    return form

def scale_form(form,r):

    k_i     = form.key_index()
    uv_i    = form.uv_index()
    vcount  = len(form.vertex)
    anchors = list(form.anchors())
    fixed   = list(form.fixed())
    fixed   = set(anchors + fixed)
    fixed   = [k_i[key] for key in fixed]
    free    = list(set(range(vcount)) - set(fixed))
    edges   = [(k_i[u], k_i[v]) for u, v in form.edges_where({'is_edge': True})]
    xyz     = array(form.get_vertices_attributes('xyz'), dtype=float64)
    p       = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)
    q       = [attr.get('q', 1.0) for u, v, attr in form.edges_where({'is_edge': True}, True)]
    q       = array(q, dtype=float64).reshape((-1, 1))
    C       = connectivity_matrix(edges, 'csr')
    Ci      = C[:, free]
    Cf      = C[:, fixed]
    Cit     = Ci.transpose()
    Ct      = C.transpose()

    q = q * r
    Q = diags([q.ravel()], [0])

    A       = Cit.dot(Q).dot(Ci)
    B       = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A,p[free, 2] - B.dot(xyz[fixed, 2]))

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z']  = xyz[index, 2]

    for u, v, attr in form.edges_where({'is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['q'] = q[index, 0]

    return form

def update_form(form,q):

    k_i     = form.key_index()
    uv_i    = form.uv_index()
    vcount  = len(form.vertex)
    anchors = list(form.anchors())
    fixed   = list(form.fixed())
    fixed   = set(anchors + fixed)
    fixed   = [k_i[key] for key in fixed]
    free    = list(set(range(vcount)) - set(fixed))
    edges   = [(k_i[u], k_i[v]) for u, v in form.edges_where({'is_edge': True})]
    xyz     = array(form.get_vertices_attributes('xyz'), dtype=float64)
    p       = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)
    q       = array(q, dtype=float64).reshape((-1, 1))
    C       = connectivity_matrix(edges, 'csr')
    Ci      = C[:, free]
    Cf      = C[:, fixed]
    Cit     = Ci.transpose()
    Ct      = C.transpose()
    Q = diags([q.ravel()], [0])

    A       = Cit.dot(Q).dot(Ci)
    B       = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A,p[free, 2] - B.dot(xyz[fixed, 2]))

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z']  = xyz[index, 2]

    for u, v, attr in form.edges_where({'is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['q'] = q[index, 0]

    return form

def paralelise_form(form, force, q, alpha = 1.0, kmax = 100, plot = None, display = False):

    # Update constraints in edges of Form

    uv_i = form.uv_index()

    for u,v in form.edges_where({'is_edge': True, 'is_external': False}):
        i = uv_i[(u,v)]
        key = (u,v)
        form.set_edge_attribute(key,'q',value=q[i])
        # print(q[i])
        a = form.vertex_coordinates(u)
        b = form.vertex_coordinates(v)
        lh = distance_point_point_xy(a,b)
        f_target = q[i]*lh
        # print(f_target)
        form.set_edge_attribute(key,'fmin',value=f_target)
        form.set_edge_attribute(key,'fmax',value=f_target)
        # form.set_edge_attribute(key,'lmin',value=lh)
        # form.set_edge_attribute(key,'lmax',value=lh)

    if plot:
        plot_form(form).show()
        plot_force(force).show()

    # Initialize

    k_i     = form.key_index()
    uv_i    = form.uv_index()
    vcount  = len(form.vertex)
    anchors = list(form.anchors())
    fixed   = list(form.fixed())
    fixed   = set(anchors + fixed)
    fixed   = [k_i[key] for key in fixed]
    free    = list(set(range(vcount)) - set(fixed))
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True})]
    xy    = array(form.get_vertices_attributes('xy'), dtype=float64)
    lmin  = array([attr.get('lmin', 1e-7) for u, v, attr in form.edges_where({'is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    lmax  = array([attr.get('lmax', 1e+7) for u, v, attr in form.edges_where({'is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    fmin  = array([attr.get('fmin', 1e-7) for u, v, attr in form.edges_where({'is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    fmax  = array([attr.get('fmax', 1e+7) for u, v, attr in form.edges_where({'is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    C     = connectivity_matrix(edges, 'csr')
    Ct    = C.transpose()
    CtC   = Ct.dot(C)
    Ci      = C[:, free]
    Cf      = C[:, fixed]
    Cit     = Ci.transpose()

    # force = ForceDiagram.from_formdiagram(form)

    _k_i   = force.key_index()
    _fixed = list(force.fixed())
    _fixed = [_k_i[key] for key in _fixed]
    _fixed = _fixed or [0]
    _edges = force.ordered_edges(form)
    _xy    = array(force.get_vertices_attributes('xy'), dtype=float64)
    _C     = connectivity_matrix(_edges, 'csr')
    _Ct    = _C.transpose()
    _Ct_C  = _Ct.dot(_C)

    _xy[:] = rot90(_xy, +1.0)

    uv  = C.dot(xy)
    _uv = _C.dot(_xy)
    l   = normrow(uv)
    _l  = normrow(_uv)

    t   = alpha * normalizerow(uv) + (1 - alpha) * normalizerow(_uv)

    # Paralelize

    for k in range(kmax):
        # apply length bounds
        apply_bounds(l, lmin, lmax)
        apply_bounds(_l, fmin, fmax)
        # print, if allowed
        if display:
            print(k)
        if alpha != 1.0:
            # if emphasis is not entirely on the form
            # update the form diagram
            xy = parallelise_sparse(CtC, Ct.dot(l * t), xy, fixed, 'CtC')
            uv = C.dot(xy)
            l  = normrow(uv)
        if alpha != 0.0:
            # if emphasis is not entirely on the force
            # update the force diagram
            _xy = parallelise_sparse(_Ct_C, _Ct.dot(_l * t), _xy, _fixed, '_Ct_C')
            _uv = _C.dot(_xy)
            _l  = normrow(_uv)

    f = _l
    q = (f / l).astype(float64)
    q = q.reshape(-1, 1)
    # print('Final qs')
    # print(q)

    _xy[:] = rot90(_xy, -1.0)

    a = [angle_vectors_xy(uv[i], _uv[i], deg=True) for i in range(len(edges))]

    # print(a)

    for key, attr in form.vertices(True):
        i = k_i[key]
        attr['x'] = xy[i, 0]
        attr['y'] = xy[i, 1]
    for u, v, attr in form.edges_where({'is_edge': True}, True):
        i = uv_i[(u, v)]
        attr['q'] = q[i, 0]
        attr['f'] = f[i, 0]
        attr['l'] = l[i, 0]
        attr['a'] = a[i]


    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    # Update Z

    xyz     = array(form.get_vertices_attributes('xyz'), dtype=float64)
    q =     q.reshape(-1, 1)
    Q =     diags([q.ravel()], [0])
    p       = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)

    A       = Cit.dot(Q).dot(Ci)
    A_solve = factorized(A)
    B       = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = A_solve(p[free, 2] - B.dot(xyz[fixed, 2]))

    l  = normrow(C.dot(xyz))
    f  = q * l
    r  = C.transpose().dot(Q).dot(C).dot(xyz) - p

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z']  = xyz[index, 2].tolist()
        attr['rx'] = r[index, 0]
        attr['ry'] = r[index, 1]
        attr['rz'] = r[index, 2]
    for u, v, attr in form.edges_where({'is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['f'] = f[index, 0]

    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    if plot:
        plot_form(form).show()
        plot_force(force).show()

    return form, force

def reactions(form, show=False):

    from compas.geometry.basic import scale_vector
    from compas.geometry.basic import normalize_vector
    from compas.geometry.basic import sum_vectors
    from compas.geometry.basic import norm_vector
    from compas.geometry.distance import distance_point_point_xy

    eq_node = {}
    sum_res = 0

    for key, nbrs in iter(form.halfedge.items()):
        if form.vertex[key]['is_fixed'] is True:
            vectors = []
            for nbr, face in iter(nbrs.items()):
                qi = form.get_edge_attribute((key, nbr), 'q')
                a = form.vertex_coordinates(key)
                b = form.vertex_coordinates(nbr)
                f = distance_point_point_xy(a,b)*qi
                vectors.append(scale_vector(normalize_vector([a[0]-b[0],a[1]-b[1],0]),f))
            sum_vec = sum_vectors(vectors)
            eq_node[key] = norm_vector(sum_vec)
            sum_res += eq_node[key]
            form.set_vertex_attribute(key, 'rx', value=sum_vec[0])
            form.set_vertex_attribute(key, 'ry', value=sum_vec[1])
            form.set_vertex_attribute(key, 'r', value=norm_vector(sum_vec))
            # if len(vectors) > 1:
            rx = sum_vec[0]
            ry = sum_vec[1]
            r = norm_vector(sum_vec)

    # print('Reaction forces are rx: {0:.3f} and ry: {1:.3f} - Total: {2:.3f}'.format(rx,ry, r))
    if show:
        plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
        plotter.draw_vertices(text=eq_node)
        plotter.draw_edges()
        plotter.show()

    return



