
from numpy import array
from numpy import float64
from numpy import newaxis
from numpy import zeros
from numpy import hstack

from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import splu
from scipy.sparse.linalg import lsqr
from scipy.sparse import diags

from compas.numerical import fd_numpy
from compas.numerical import connectivity_matrix

from compas_plotters import MeshPlotter


__all__ = [
    'z_from_form',
    'z_update',
    'zlq_from_qid',
    'zq_from_qid',
    'q_from_qid',
    'xyz_from_q',
    'reactions',
]


def z_from_form(form):
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

    k_i = form.key_index()
    xyz = form.vertices_attributes(('x', 'y', 'z'))
    loads = form.vertices_attributes(('px', 'py', 'pz'))
    q = form.edges_attribute('q')
    fixed = form.vertices_where({'is_fixed': True})
    fixed = [k_i[k] for k in fixed]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges()]

    # compute equilibrium
    # update the mesh geometry

    xyz, q, f, l, r = fd_numpy(xyz, edges, fixed, q, loads)

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['x'] = xyz[index, 0]
        attr['y'] = xyz[index, 1]
        attr['z'] = xyz[index, 2]

    return form


def z_update(form, zmax=None):
    """ Built-in update of the heights in the Form-Diagram (only z-coordinates can change).

    Parameters
    ----------
    form : obj
        The FormDiagram.
    zmax : float, optional
        The maximum height of the diagram.
        By default no maximum is set.

    Returns
    -------
    form : obj
        The scaled form diagram.

    """

    k_i = form.key_index()
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

    Q = diags([q.ravel()], [0])

    A = Cit.dot(Q).dot(Ci)
    B = Cit.dot(Q).dot(Cf)

    xyz[free, 2] = spsolve(A, p[free, 2] - B.dot(xyz[fixed, 2]))

    if zmax:
        factor = max(xyz[free, 2])/zmax
        q = q * factor
        xyz[free, 2] = xyz[free, 2]/factor

        for index, edge in enumerate(form.edges_where({'_is_edge': True})):
            form.edge_attribute(edge, 'q', q[index].item())

    for key, attr in form.vertices(True):
        index = k_i[key]
        attr['z'] = xyz[index, 2]

    return form


def zlq_from_qid(qid, args):  # deprecated remove on next clean up
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


def zq_from_qid(qid, args):  # deprecated remove on next clean up
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


def q_from_qid(qid, args):  # deprecated remove on next clean up
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
    try:
        SPLU_D = splu(CiQCi)
        Xfree = SPLU_D.solve(Pi - CiQCb.dot(Xb))
    except BaseException:
        A = CiQCi
        b = Pi - CiQCb.dot(Xb)
        resx = lsqr(A, b[:, 0])
        resy = lsqr(A, b[:, 1])
        resz = lsqr(A, b[:, 2])
        print('* Warning: System might be bad conditioned - recured to LSQR')
        xfree = resx[0].reshape(-1, 1)
        yfree = resy[0].reshape(-1, 1)
        zfree = resz[0].reshape(-1, 1)
        Xfree = hstack([xfree, yfree, zfree])

    return Xfree


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
    None
        Form Diagram is modified in place.

    """

    # Mapping

    k_i = form.key_index()
    i_k = form.index_key()

    # Vertices and edges

    n = form.number_of_vertices()
    fixed = [k_i[key] for key in form.fixed()]
    # rol = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]

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
