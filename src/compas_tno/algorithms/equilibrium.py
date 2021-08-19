
from numpy import array
from numpy import float64
from numpy import newaxis
from numpy import zeros

from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import splu
from scipy.sparse import diags

from compas.numerical import fd_numpy
from compas.numerical import connectivity_matrix

from compas_plotters import MeshPlotter


__all__ = [
    'zlq_from_qid',
    'zlq_from_q',
    'zq_from_qid',
    'q_from_qid',
    'z_update',
    'xyz_from_q',
    'z_from_form',
    'reactions',
    'apply_sag'
]


def z_from_form(mesh):
    """ Relaxation of Form-Diagram. FDM with 'q's stored in the form (All coordinates can change).

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram with attributes relevant to the analysis (force densities, loads, fixities, etc).

    Returns
    -------
    form : FormDiagram
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


def scale_fdm(form, r):
    """ scale the FormDiagram of a factor r using FDM (all coordinates can change).

    Parameters
    ----------
    form : obj
        The FormDiagram.
    r : float
        The scaling factor on force densities.

    Returns
    -------
    form : obj
        The scaled form diagram.

    """

    uv_i = form.uv_index()
    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True, '_is_external': False})])[:, newaxis]
    q = q * r

    for u, v in form.edges_where({'_is_external': False}):
        if form.edge_attribute((u, v), '_is_edge') is True:
            i = uv_i[(u, v)]
            [qi] = q[i]
            form.edge_attribute((u, v), 'q', value=qi)

    form = z_from_form(form)

    return form


def scale_form(form, r):
    """ Scale the FormDiagram of a factor r using built-in FDM (only z-coordinates can change).

    Parameters
    ----------
    form : obj
        The FormDiagram.
    r : float
        The scaling factor on force densities.

    Returns
    -------
    form : obj
        The scaled form diagram.

    """

    from numpy import float64
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve

    k_i = form.key_index()
    uv_i = form.uv_index()
    vcount = len(form.vertex)
    anchors = list(form.anchors())
    fixed = list(form.fixed())
    fixed = set(anchors + fixed)
    fixed = [k_i[key] for key in fixed]
    free = list(set(range(vcount)) - set(fixed))
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    xyz = array(form.vertices_attributes('xyz'), dtype=float64)
    p = array(form.vertices_attributes(('px', 'py', 'pz')), dtype=float64)
    q = [form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]
    q = array(q, dtype=float64).reshape((-1, 1))
    C = connectivity_matrix(edges, 'csr')
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cit = Ci.transpose()

    # TODO: Change to use function zq_from_qid

    q = q * r
    Q = diags([q.ravel()], [0])

    A = Cit.dot(Q).dot(Ci)
    B = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A, p[free, 2] - B.dot(xyz[fixed, 2]))

    i = 0
    for key in form.vertices():
        form.vertex_attribute(key, 'z', xyz[i, 2])
        i = i + 1

    i = 0
    for u, v in form.edges_where({'_is_edge': True}):
        form.edge_attribute((u, v), 'q', q[i, 0])
        i = i + 1

    return


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
        pz = form.vertex_attribute(key, 'pz')
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


def apply_sag(form, boundary_force=10.0, signe_compression=-1.0):

    for u, v in form.edges():
        form.edge_attribute((u, v), 'q', signe_compression*1.0)

    for u, v in form.edges_on_boundary():
        form.edge_attribute((u, v), 'q', signe_compression*boundary_force)

    z_from_form(form)

    for key in form.vertices():
        form.vertex_attribute(key, 'z', 0.0)

    return form
