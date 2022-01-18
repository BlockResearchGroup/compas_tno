
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


def equilibrium_fdm(form):
    """ Compute equilibrium of the form diagram using the force density method (FDM) with 'q's stored in the form (All coordinates can change).

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


def vertical_equilibrium_fdm(form, zmax=None):
    """ Compute equilibrium of the form diagram using the force density method (FDM) and update only the z-coordinates.

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram.
    zmax : float, optional
        The maximum height of the diagram.
        By default no maximum is set.

    Returns
    -------
    form : FormDiagram
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


def q_from_qid(q, ind, Edinv, Ei, ph):
    r""" Calculate q's from all qid's.

    Parameters
    ----------
    q : array [mx1]
        Force density vector.
    ind : list [k]
        List with the indices of the k independent edges.
    Edinv : array [nix(m - k)]
        Inverse of equilibrium matrix sliced to the dependent edges.
    Ei : array [nix(m - k)]
        Equilibrium matrix sliced to the independent edges.
    ph : array [nix1]
        Stack of the horizontal loads applied to the free nodes [pxi, pyi].

    Returns
    -------
    q : array
        Force densities on all edges.

    Notes
    -------
    Solver of thhe following equation:
    $\mathbf{q}_\mathrm{id} = - mathbf{Ed}_\mathrm{id}$(\mathbf{E}_\mathrm{i} * \mathbf{q}_\mathrm{i} - \mathbf{p}_\mathrm{h})

    Reference
    ---------
    Block and Lachauer, 2014...

    """

    m = len(q)
    dep = list(set(range(m)) - set(ind))
    q[dep] = Edinv.dot(- ph + Ei.dot(q[ind]))

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


def compute_reactions(form, plot=False):
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

    return


def xyz_from_xopt(variables, M):
    """Compute the nodal position and relevant parameters and from the variables of the optimisation and the class of matrices (M)

    Parameters
    ----------
    variables : array
        The ``n`` variables of the optimisation process
    M : class
        The relevant matrices to be stored

    Returns
    -------
    M
        The updated relevant matrices
    """

    if isinstance(M, list):
        M = M[0]

    # variables
    k = M.k  # number of force variables
    n = M.n  # number of vertices
    nb = len(M.fixed)  # number of fixed vertices
    t = M.shape.datashape['t']

    qid = variables[:k]
    check = k
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[check:check + 2*nb]
        check = check + 2*nb
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
        nbxy = nb
    if 'zb' in M.variables:
        zb = variables[check: check + nb]
        check = check + nb
        M.X[M.fixed, [2]] = zb.flatten()
    if 't' in M.variables or 'n' in M.variables:
        thk = variables[check: check + 1]
        check = check + 1
    elif 'lambd' in M.variables:
        lambd = variables[check: check + 1]
        M.P[:, [0]] = lambd * M.px0
        M.P[:, [1]] = lambd * M.py0
        check = check + 1
    if 'tub' in M.variables:
        tub = variables[check: check + n]
        M.tub = tub
        check = check + n
    if 'tlb' in M.variables:
        tlb = variables[check: check + n]
        M.tlb = tlb
        check = check + n
    if 'tub_reac' in M.variables:
        tub_reac = variables[check: check + 2*nb].reshape(-1, 1)
        M.tub_reac = tub_reac
        check = check + 2*nb

    q = M.q
    Xfixed = M.X[M.fixed]

    # update geometry
    M.X[M.free] = xyz_from_q(q, M.P[M.free], Xfixed, M.Ci, M.Cit, M.Cb)

    return M
