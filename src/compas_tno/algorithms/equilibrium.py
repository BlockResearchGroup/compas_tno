
from numpy import array
from numpy import float64
from numpy import newaxis
from numpy import zeros
from numpy import hstack

from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import splu
from scipy.sparse import diags
from scipy.sparse.linalg import factorized
from compas.numerical import normrow
from compas.numerical import normalizerow
from compas.numerical import spsolve_with_known
from compas_tna.utilities import apply_bounds
from compas_tna.utilities import parallelise_sparse
from compas.numerical import fd_numpy

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas_tna.equilibrium import horizontal
from compas_tna.utilities import rot90
from compas.geometry import angle_vectors_xy

from compas.numerical import connectivity_matrix
from compas_plotters import MeshPlotter
from compas.geometry import distance_point_point_xy

from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_force


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'zlq_from_qid',
    'zlq_from_q',
    'zq_from_qid',
    'q_from_qid',
    'z_update',
    'xyz_from_q',
    'z_from_form',
    'update_tna',
    'force_update_from_form',
    'update_form',
    'paralelise_form',
    'reactions',
    'apply_sag'
]


def z_from_form(mesh):
    """ Relaxation of Form-Diagram. FDM with 'q's stored in the form (All coordinates can change).

    Parameters
    ----------
    form : obj
        The FormDiagram.

    Returns
    -------
    form : obj
        The relaxed form diagram.

    """

    # preprocess

    k_i = mesh.key_index()
    xyz = mesh.vertices_attributes(('x', 'y', 'z'))
    loads = mesh.vertices_attributes(('px', 'py', 'pz'))
    q = mesh.edges_attribute('q')
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
    """ Calculate z's from independent edges.

    Parameters
    ----------
    qid : list
        Force densities of the independent edges.
    args : tuple
        Arrays and matrices relevant to the operation.


    Returns
    -------
    z : array
        Heights of the nodes
    l2 : array
        Lenghts squared
    q : array
        Force densities without symetrical edges (q[sym] = 0)
    q_ : array
        Force densities with symetrical edges

    """

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym = args[:22]
    q[ind] = array(qid).reshape(-1, 1)
    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    q_ = 1 * q
    q[sym] *= 0

    # if not planar:
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))
    l2 = lh + C.dot(z)**2

    return z, l2, q, q_


def zq_from_qid(qid, args):
    """ Calculate z's from independent edges.

    Parameters
    ----------
    qid : list
        Force densities of the independent edges.
    args : tuple
        Arrays and matrices relevant to the operation.


    Returns
    -------
    z : array
        Heights of the nodes
    q : array
        Force densities

    """

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed = args[:20]
    q[ind] = array(qid).reshape(-1, 1)
    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

    return z, q


def q_from_qid(qid, args):
    """ Calculate q's from all qid's.

    Parameters
    ----------
    qid : list
        Force densities of the independent edges.
    args : tuple
        Arrays and matrices relevant to the operation.


    Returns
    -------
    q : array
        Force densities on all edges.

    """

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p = args[:14]
    q[ind] = array(qid).reshape(-1, 1)
    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))

    return q


def zlq_from_q(q, args):
    """ Calculate z's from all q's.

    Parameters
    ----------
    q : array
        Force densities of all edges.
    args : tuple
        Arrays and matrices relevant to the operation.


    Returns
    -------
    z : array
        Heights of the nodes
    l2 : array
        Lenghts squared
    q : array
        Force densities without symetrical edges (q[sym] = 0)
    q_ : array
        Force densities with symetrical edges

    """

    q_old, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym = args[:22]
    q_ = 1 * q
    q[sym] *= 0
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))
    l2 = lh + C.dot(z)**2

    return z, l2, q, q_


def xyz_from_q(q, Pi, Xb, Ci, Cit, Cb):
    """ Calculate coordinates's xyz the q's.

    Parameters
    ----------
    q : array (m)
        Force densities of all edges.
    Pi : array (ni x 3)
        External forces applied to the free vertices.
    Xb : array (nb x 3)
        Position of the fixed vertices.
    Ci : array (m x ni)
        Connectivity matrix on the free vertices
    Cit : array (ni x m)
        Transpose of connectivity matrix on the free vertices
    Cb : array (m x nb)
        Connectivity matrix on the fixed vertices

    Returns
    -------
    Xfree : array (ni)
        x, y, z coordinates of the nodes

    """

    CiQCb = Cit.dot(diags(q.flatten())).dot(Cb)
    CiQCi = Cit.dot(diags(q.flatten())).dot(Ci)
    SPLU_D = splu(CiQCi)
    Xfree = SPLU_D.solve(Pi - CiQCb.dot(Xb))

    return Xfree


def z_update(form):
    """ Built-in update of the heights in the Form-Diagram (only z-coordinates can change).

    Parameters
    ----------
    form : obj
        The FormDiagram.

    Returns
    -------
    form : obj
        The scaled form diagram.

    """

    k_i = form.key_index()
    uv_i = form.uv_index()
    vcount = len(form.vertex)
    anchors = list(form.anchors())
    fixed = list(form.fixed())
    fixed = set(anchors + fixed)
    fixed = [k_i[key] for key in fixed]
    free = list(set(range(vcount)) - set(fixed))
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    xyz = array(form.get_vertices_attributes('xyz'), dtype=float64)
    p = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)
    q = [attr.get('q', 1.0) for u, v, attr in form.edges_where({'_is_edge': True}, True)]
    q = array(q, dtype=float64).reshape((-1, 1))
    C = connectivity_matrix(edges, 'csr')
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cit = Ci.transpose()

    Q = diags([q.ravel()], [0])

    A = Cit.dot(Q).dot(Ci)
    B = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A, p[free, 2] - B.dot(xyz[fixed, 2]))

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z'] = xyz[index, 2]

    return form


def update_tna(form, delete_face=True, plots=False, save=False):

    if delete_face:
        form.delete_face(0)

    corners = list(form.vertices_where({'is_fixed': True}))
    print(form)
    form.set_vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
    form.update_boundaries(feet=2)

    # for key in form.edges_where({'_is_external': True}):
    #     form.edge_attribute(key,'q',value=x_reaction)
    #     form.edge_attribute(key,'fmin',value=x_reaction)
    #     form.edge_attribute(key,'fmax',value=x_reaction)
    #     form.edge_attribute(key,'lmin',value=1.00)
    #     form.edge_attribute(key,'lmax',value=1.00)

    for u, v in form.edges_where({'_is_external': False}):
        qi = form.edge_attribute((u, v), 'q')
        a = form.vertex_coordinates(u)
        b = form.vertex_coordinates(v)
        lh = distance_point_point_xy(a, b)
        form.edge_attribute((u, v), 'fmin', value=qi*lh)
        form.edge_attribute((u, v), 'fmax', value=qi*lh)
        form.edge_attribute((u, v), 'lmin', value=lh)
        form.edge_attribute((u, v), 'lmax', value=lh)

    force = ForceDiagram.from_formdiagram(form)
    horizontal(form, force, display=False)
    # Vertical?

    if plots:
        plot_force(force, form, radius=0.05).show()
        plot_form(form, radius=0.05).show()

    # st = 'discretize/02_complete_1div_complete'

    if save:
        force.to_obj('/Users/mricardo/compas_dev/compas_loadpath/data/' + st + '_force.obj')
        form.to_json('/Users/mricardo/compas_dev/compas_loadpath/data/' + st + '_tna.json')

    return form, force


def force_update_from_form(force, form):
    """Update the force diagram after modifying the (force densities of) the form diagram.

    Parameters
    ----------
    force : :class:`ForceDiagram`
        The force diagram on which the update is based.
    form : :class:`FormDiagram`
        The form diagram to update.

    Reference
    -------

    See ``compas_tna`` package.

    Returns
    -------
    None
        The form and force diagram are updated in-place.
    """
    # --------------------------------------------------------------------------
    # form diagram
    # --------------------------------------------------------------------------
    vertex_index = form.vertex_index()

    xy = array(form.xy(), dtype=float64)
    edges = [[vertex_index[u], vertex_index[v]] for u, v in form.edges_where({'_is_edge': True})]
    C = connectivity_matrix(edges, 'csr')
    Q = diags([form.q()], [0])
    uv = C.dot(xy)
    # --------------------------------------------------------------------------
    # force diagram
    # --------------------------------------------------------------------------
    _vertex_index = force.vertex_index()

    _known = [_vertex_index[force.anchor()]]
    _xy = array(force.xy(), dtype=float64)
    _edges = force.ordered_edges(form)
    _edges[:] = [(_vertex_index[u], _vertex_index[v]) for u, v in _edges]
    _C = connectivity_matrix(_edges, 'csr')
    _Ct = _C.transpose()
    print(C.shape, _C.shape, Q.shape, uv.shape)
    # --------------------------------------------------------------------------
    # compute reciprocal for given q
    # --------------------------------------------------------------------------
    _xy = spsolve_with_known(_Ct.dot(_C), _Ct.dot(Q).dot(uv), _xy, _known)
    # --------------------------------------------------------------------------
    # update force diagram
    # --------------------------------------------------------------------------
    for vertex, attr in force.vertices(True):
        index = _vertex_index[vertex]
        attr['x'] = _xy[index, 0]
        attr['y'] = _xy[index, 1]

def update_form(form, q):

    k_i = form.key_index()
    uv_i = form.uv_index()
    vcount = len(form.vertex)
    anchors = list(form.anchors())
    fixed = list(form.fixed())
    fixed = set(anchors + fixed)
    fixed = [k_i[key] for key in fixed]
    free = list(set(range(vcount)) - set(fixed))
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    xyz = array(form.get_vertices_attributes('xyz'), dtype=float64)
    p = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)
    q = array(q, dtype=float64).reshape((-1, 1))
    C = connectivity_matrix(edges, 'csr')
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cit = Ci.transpose()
    Ct = C.transpose()
    Q = diags([q.ravel()], [0])

    A = Cit.dot(Q).dot(Ci)
    B = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A, p[free, 2] - B.dot(xyz[fixed, 2]))

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z'] = xyz[index, 2]

    for u, v, attr in form.edges_where({'_is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['q'] = q[index, 0]

    return form


def paralelise_form(form, force, q, alpha=1.0, kmax=100, plot=None, display=False):

    # Update constraints in edges of Form

    uv_i = form.uv_index()

    for u, v in form.edges_where({'_is_edge': True, '_is_external': False}):
        i = uv_i[(u, v)]
        key = (u, v)
        form.edge_attribute(key, 'q', value=q[i])
        # print(q[i])
        a = form.vertex_coordinates(u)
        b = form.vertex_coordinates(v)
        lh = distance_point_point_xy(a, b)
        f_target = q[i]*lh
        # print(f_target)
        form.edge_attribute(key, 'fmin', value=f_target)
        form.edge_attribute(key, 'fmax', value=f_target)
        # form.edge_attribute(key,'lmin',value=lh)
        # form.edge_attribute(key,'lmax',value=lh)

    if plot:
        plot_form(form).show()
        plot_force(force, form).show()

    # Initialize

    k_i = form.key_index()
    uv_i = form.uv_index()
    vcount = len(form.vertex)
    anchors = list(form.anchors())
    fixed = list(form.fixed())
    fixed = set(anchors + fixed)
    fixed = [k_i[key] for key in fixed]
    free = list(set(range(vcount)) - set(fixed))
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True})]
    xy = array(form.get_vertices_attributes('xy'), dtype=float64)
    lmin = array([attr.get('lmin', 1e-7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    lmax = array([attr.get('lmax', 1e+7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    fmin = array([attr.get('fmin', 1e-7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    fmax = array([attr.get('fmax', 1e+7) for u, v, attr in form.edges_where({'_is_edge': True}, True)], dtype=float64).reshape((-1, 1))
    C = connectivity_matrix(edges, 'csr')
    Ct = C.transpose()
    CtC = Ct.dot(C)
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cit = Ci.transpose()

    # force = ForceDiagram.from_formdiagram(form)

    _k_i = force.key_index()
    _fixed = list(force.fixed())
    _fixed = [_k_i[key] for key in _fixed]
    _fixed = _fixed or [0]
    _edges = force.ordered_edges(form)
    _xy = array(force.get_vertices_attributes('xy'), dtype=float64)
    _C = connectivity_matrix(_edges, 'csr')
    _Ct = _C.transpose()
    _Ct_C = _Ct.dot(_C)

    _xy[:] = rot90(_xy, +1.0)

    uv = C.dot(xy)
    _uv = _C.dot(_xy)
    l = normrow(uv)
    _l = normrow(_uv)

    t = alpha * normalizerow(uv) + (1 - alpha) * normalizerow(_uv)

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
            l = normrow(uv)
        if alpha != 0.0:
            # if emphasis is not entirely on the force
            # update the force diagram
            _xy = parallelise_sparse(_Ct_C, _Ct.dot(_l * t), _xy, _fixed, '_Ct_C')
            _uv = _C.dot(_xy)
            _l = normrow(_uv)

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
    for u, v, attr in form.edges_where({'_is_edge': True}, True):
        i = uv_i[(u, v)]
        attr['q'] = q[i, 0]
        attr['f'] = f[i, 0]
        attr['_l'] = l[i, 0]
        attr['a'] = a[i]

    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    # Update Z

    xyz = array(form.get_vertices_attributes('xyz'), dtype=float64)
    q = q.reshape(-1, 1)
    Q = diags([q.ravel()], [0])
    p = array(form.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)

    A = Cit.dot(Q).dot(Ci)
    A_solve = factorized(A)
    B = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = A_solve(p[free, 2] - B.dot(xyz[fixed, 2]))

    l = normrow(C.dot(xyz))
    f = q * l
    r = C.transpose().dot(Q).dot(C).dot(xyz) - p

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z'] = xyz[index, 2].tolist()
        attr['_rx'] = r[index, 0]
        attr['_ry'] = r[index, 1]
        attr['_rz'] = r[index, 2]
    for u, v, attr in form.edges_where({'_is_edge': True}, True):
        index = uv_i[(u, v)]
        attr['f'] = f[index, 0]

    for key, attr in force.vertices(True):
        i = _k_i[key]
        attr['x'] = _xy[i, 0]
        attr['y'] = _xy[i, 1]

    if plot:
        plot_form(form).show()
        plot_force(force, form).show()

    return form, force


def reactions(form, plot=False):
    """ Compute and plot the reaction on the supports.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to calculate the reactions.
    plot : bool
        True/False to plot the reactions.


    Returns
    -------

    """

    # Mapping

    k_i = form.key_index()
    i_k = form.index_key()

    # Vertices and edges

    n = form.number_of_vertices()
    fixed = [k_i[key] for key in form.fixed()]
    rol = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    free = list(set(range(n)) - set(fixed) - set(rol))

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    px = zeros((n, 1))
    py = zeros((n, 1))
    pz = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)

    # C and E matrices

    C = connectivity_matrix(edges, 'csr')
    uvw = C.dot(xyz)
    U = uvw[:, 0]
    V = uvw[:, 1]
    W = uvw[:, 2]
    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]

    # Horizontal checks

    Rx = C.transpose().dot(U * q.ravel()) - px.ravel()
    Ry = C.transpose().dot(V * q.ravel()) - py.ravel()
    Rz = C.transpose().dot(W * q.ravel()) - pz.ravel()

    eq_node = {key: [round(Rx[k_i[key]], 1), round(Ry[k_i[key]], 1)] for key in form.vertices_where({'is_fixed': True})}

    for i in fixed:
        key = i_k[i]
        form.vertex_attribute(key, '_rx', value=Rx[i])
        form.vertex_attribute(key, '_ry', value=Ry[i])
        form.vertex_attribute(key, '_rz', value=Rz[i])
        if plot:
            print('Reactions in key: {0} are:'.format(key))
            print(Rx[i], Ry[i], Rz[i])

    for key in form.vertices_where({'rol_x': True}):
        i = k_i[key]
        eq_node[key] = round(Rx[i], 2)
        form.vertex_attribute(key, '_rx', value=Rx[i])
        if plot:
            print('Reactions in Partial X-Key: {0} :'.format(key))
            print(Rx[i])

    for key in form.vertices_where({'rol_y': True}):
        i = k_i[key]
        eq_node[key] = round(Ry[i], 2)
        form.vertex_attribute(key, '_ry', value=Ry[i])
        if plot:
            print('Reactions in Partial Y-Key: {0} :'.format(key))
            print(Ry[i])

    if plot:
        plotter = MeshPlotter(form, figsize=(10, 7), fontsize=8)
        plotter.draw_vertices(text=eq_node)
        plotter.draw_edges()
        plotter.show()

    return


def apply_sag(form, boundary_force=10.0):

    for u, v in form.edges():
        form.edge_attribute((u, v), 'q', 1.0)

    for u, v in form.edges_on_boundary():
        form.edge_attribute((u, v), 'q', boundary_force)

    z_from_form(form)

    for key in form.vertices():
        form.vertex_attribute(key, 'z', 0.0)

    return form
