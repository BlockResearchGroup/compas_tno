from numpy import array
from numpy import zeros
from numpy import vstack
from numpy import hstack
from numpy import newaxis
from numpy import identity
from numpy.linalg import pinv
from numpy.linalg import matrix_rank

from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import vstack as svstack

from compas.numerical import connectivity_matrix
from compas.numerical import normrow

from compas.utilities import geometric_key

from compas_tno.algorithms import find_independents
from compas_tno.algorithms import check_independents
from compas_tno.algorithms import check_horizontal

from compas_tno.utilities import apply_radial_symmetry
from compas_tno.utilities import apply_symmetry_from_axis
from compas_tno.utilities import find_sym_axis_in_rect_patterns
from compas_tno.utilities import build_symmetry_transformation

import time


class Problem():

    """
    The ``Problem`` class stores the matrices used in the optimisation. These are listed as parameters of the class and described below.

    Parameters
    ----------

    ``q`` : array(m x 1)
        The vector of force densities
    m : int
        The number of edges
    n : int
        The number of vertices
    nb : int
        The number of fixed vertices
    E : array(2n x m)
        The horizontal equilibrium matrix
    C : array(n x m)
        The connectivity matrix
    Ct : array(m x n)
        The transposed connectivity matrix
    Ci : array(m x ni)
        The sliced connectivity matrix with regards to the internal vertices
    Cit : array(ni x m)
        The transpose of Ci
    Cb : array(nb x m)
        The sliced connectivity matrix with regards to the constrained (support) vertices
    U : array(m x m)
        The diagonal matrix of coorrdinate differences in the x-direction
    V : array(m x m)
        The diagonal matrix of coorrdinate differences in the y-direction
    P : array(n x 3)
        The applied external loads in all nodes
    free : list
        The list with the index of the free vertices
    fixed : list
        The list with the index of the fixed vertices
    phfree : array(ni x 2)
        The applied horizonal loads on the free vertices
    ph : array(ni x 2)
        A column vector with the the applied horizonal loads on the free vertices
    lb : array(n x 1)
        The lower-bound (intrados) limit for the nodes
    ub : array(n x 1)
        The upper-bound (extrados) limit for the nodes
    lb0 : array(n x 1)
        The lower-bound (intrados) limit for the nodes at the starting point (can be modified if thickness in the objective)
    ub0 : array(n x 1)
        The upper-bound (extrados) limit for the nodes at the starting point (can be modified if thickness in the objective)
    s : array(n x 1)
        The middle surface of the nodes
    X : array(n x 3)
        The nodal position of the vertices of the network
    x0 : array(n x 1)
        The x-position of the vertices in the network
    y0 : array(n x 1)
        The y-position of the vertices in the network
    free_x : list
        index of the vertices free to move in x
    free_y : list
        index of the vertices free to move in y
    rol_x : list
        index of the vertices constrained partially on x.
    rol_y : list
        index of the vertices constrained partially on y.
    Citx : array
        Connectivity matrix sliced transposed on the vertices free on x.
    City : array
        Connectivity matrix sliced transposed on the vertices free on y.
    Cbtx : array
        Connectivity matrix sliced transposed on the vertices partially fixed on x.
    Cbty : array
        Connectivity matrix sliced transposed on the vertices partially fixed on y.
    xlimits : array(n x 1)
        Limits on the x-direction in which the nodes can move
    ylimits : array(n x 1)
        Limits on the y-direction in which the nodes can move
    qmin : array(m x 1)
        Lower-bounds of the force densities in the edges
    qmax : array(m x 1)
        Upper-bounds of the force densities in the edges
    k_i : dict
        key-index dictionary
    uv_i : dict
        uv-index dictionary
    i_uv : dict
        index-uv dictionary
    ind : list
        List with the index of the independent edges
    k : int
        Number of independents in the problem
    dep : list
        List with the index of the dependent edges
    B : array(m x k)
        Matrix transforming the force densities in the independent edges to all force densities

    """

    def __init__(self):
        pass

    @classmethod
    def from_formdiagram(cls, form, indset=None, printout=False, find_inds=False):
        problem = cls()

        return problem


def initialise_problem(form, indset=None, printout=None, find_inds=True, tol=0.001, c=0.5):  # old function -> delete
    """ Initialise the problem for a given Form-Diagram and return the set of matrices and vectors to optimise.

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram.
    printout : bool
        Control output.
    find_inds : bool
        If True will calculate the independents
    indset : list
        Independent set to use. If empty the independents are calculated normally.
    tol: float
        Tolerance for check the independent edges equilibrium

    Returns
    -------
    args
        List of matrices and vectors used to perform the optimisation.

    """

    # Mapping

    k_i = form.key_index()
    uv_i = form.uv_index()

    # Vertices and edges

    n = form.number_of_vertices()
    m = len(list(form.edges_where({'_is_edge': True})))
    fixed = [k_i[key] for key in form.fixed()]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    sym = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
    free = list(set(range(n)) - set(fixed))

    # Constraints

    lb_ind = []
    ub_ind = []
    lb = []
    ub = []
    for key in form.vertices():
        i = k_i[key]
        lbi = form.vertex_attribute(key, 'lb')
        ubi = form.vertex_attribute(key, 'ub')
        pz = form.vertex_attribute(key, 'pz')  # WIP: take out nodes unloaded from the constraints
        if lb is not None and pz != 0.0:
            lb.append(lbi)
            lb_ind.append(i)
        if ub is not None and pz != 0.0:
            ub.append(ubi)
            ub_ind.append(i)
    lb = array(lb).reshape(-1, 1)
    ub = array(ub).reshape(-1, 1)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x = zeros((n, 1))
    y = zeros((n, 1))
    z = zeros((n, 1))
    s = zeros((n, 1))
    px = zeros((n, 1))
    py = zeros((n, 1))
    pz = zeros((n, 1))
    s = zeros((n, 1))
    w = zeros((n, 1))
    xlimits = zeros((n, 2))
    ylimits = zeros((n, 2))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i] = vertex.get('x')
        y[i] = vertex.get('y')
        z[i] = vertex.get('z')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i] = vertex.get('target', 0)
        w[i] = vertex.get('weight', 1.0)  # weight used in case of fiting...
        xlimits[i, 0] = vertex.get('xmin')
        xlimits[i, 1] = vertex.get('xmax')
        ylimits[i, 0] = vertex.get('ymin')
        ylimits[i, 1] = vertex.get('ymax')
    Wfree = diags(w[free].flatten())
    xy = xyz[:, :2]

    # If Rols are set up

    rol_x = []
    rol_y = []

    for key in form.vertices_where({'rol_x': True}):
        rol_x.append(k_i[key])

    for key in form.vertices_where({'rol_y': True}):
        rol_y.append(k_i[key])

    free_x = list(set(free) - set(rol_x))
    free_y = list(set(free) - set(rol_y))

    # C and E matrices

    C = connectivity_matrix(edges, 'csr')
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cftx = C[:, rol_x].transpose()
    Cfty = C[:, rol_y].transpose()
    Ct = C.transpose()
    Cit = Ci.transpose()
    Citx = Ct[free_x, :]  # .toarray()
    City = Ct[free_y, :]  # .toarray()
    uvw = C.dot(xyz)
    U = diags(uvw[:, 0].flatten())
    V = diags(uvw[:, 1].flatten())
    E = svstack((Citx.dot(U), City.dot(V))).toarray()

    start_time = time.time()

    # Independent and dependent branches

    if find_inds:

        if indset:
            ind = []
            for u, v in form.edges_where({'_is_edge': True}):
                if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
                    ind.append(uv_i[(u, v)])
        else:
            ind = find_independents(E)

        k = len(ind)
        dep = list(set(range(m)) - set(ind))

        elapsed_time = time.time() - start_time

        if printout:
            print('Shape Equilibrium Matrix: ', E.shape)
            print('Rank Equilibrium Matrix: ', matrix_rank(E))
            print('Found {0} independents'.format(k))
            print('# Vertices: {0} | #UB: {1} | #LB: {2}'.format(form.number_of_vertices(), len(ub), len(lb)))
            print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

        for u, v in form.edges_where({'_is_edge': True}):
            form.edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

        Edinv = -csr_matrix(pinv(E[:, dep]))
        Ei = E[:, ind]

    else:
        k = m
        ind = list(set(range(m)))
        dep = []
        Edinv = []
        Ei = []

    # Set-up

    lh = normrow(C.dot(xy))**2
    p = vstack([px[free_x], py[free_y]])
    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]

    if any(p) is True:
        check_hor = check_horizontal(E, p)
        if printout:
            if check_hor:
                print('Horizontal Loads can be taken!')
            else:
                print('Horizontal Loads are not suitable for this FD!')

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k,
            lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, xlimits, ylimits)
    args_inds = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Citx, City, Cf, U, V, p, px, py, pz, z, free_x, free_y, fixed, lh, sym, k)
    if find_inds is True:
        checked = check_independents(args_inds, tol=tol)
        if checked:
            if printout:
                print('Independents checked!')
            pass
        else:
            print('Warning: independent edges not equilibrated')

    return args


def initialise_form(form, printout=False):
    """ Initialise the problem for a Form-Diagram and return the FormDiagram with independent edges assigned and the matrices relevant to the equilibrium problem.

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram.
    printout : bool
        Whether or not printouts with the main resuslts should appear in the screen.

    Returns
    -------
    M : obj
        The object with all the matrices essential to the analysis.

    Note
    -------
    The FormDiagram is updated in place. Check ``initialise_problem_general`` for more info.

    """

    i_uv = form.index_uv()

    M = initialise_problem_general(form, printout=printout)
    adapt_problem_to_fixed_diagram(M, form, printout=printout)
    ind = M.ind

    form.update_default_edge_attributes({'is_ind': False})
    gkeys = []
    for i in ind:
        u, v = i_uv[i]
        gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
        form.edge_attribute((u, v), 'is_ind', value=True)
    form.attributes['indset'] = gkeys

    return M


def initialise_problem_torch(form, indset=None, printout=None, find_inds=True, tol=0.001):
    """ Initialise the problem for a given Form-Diagram and return the set of matrices and vectors to optimise.

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram.
    printout : bool
        Control output.
    find_inds : bool
        If True will calculate the independents
    indset : list
        Independent set to use. If empty the independents are calculated normally.
    tol: float
        Tolerance for check the independent edges equilibrium

    Returns
    -------
    args
        List of matrices and vectors used to perform the optimisation.

    """

    import torch as th

    # Mapping

    k_i = form.key_index()
    i_k = form.index_key()
    uv_i = form.uv_index()

    # Vertices and edges

    n = form.number_of_vertices()
    m = len(list(form.edges_where({'_is_edge': True})))
    fixed = [k_i[key] for key in form.fixed()]
    rol = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    sym = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
    free = list(set(range(n)) - set(fixed) - set(rol))

    # Constraints

    lb_ind = []
    ub_ind = []
    lb = []
    ub = []
    for i in range(n):
        key = i_k[i]
        if form.vertex_attribute(key, 'lb', None):
            lb_ind.append(i)
            lb.append(form.vertex_attribute(key, 'lb'))
        if form.vertex_attribute(key, 'ub', None):
            ub_ind.append(i)
            ub.append(form.vertex_attribute(key, 'ub'))

    lb = array(lb)
    ub = array(ub)
    lb.shape = (len(lb), 1)
    ub.shape = (len(ub), 1)

    # Co-ordinates and loads

    xyz = th.tensor(zeros((n, 3)))
    x = th.tensor(zeros((n, 1)))
    y = th.tensor(zeros((n, 1)))
    z = th.tensor(zeros((n, 1)))
    s = th.tensor(zeros((n, 1)))
    px = th.tensor(zeros((n, 1)))
    py = th.tensor(zeros((n, 1)))
    pz = th.tensor(zeros((n, 1)))
    s = th.tensor(zeros((n, 1)))
    w = th.tensor(zeros((n, 1)))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i] = vertex.get('x')
        y[i] = vertex.get('y')
        z[i] = vertex.get('z')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i] = vertex.get('target', 0)
        w[i] = vertex.get('weight', 1.0)  # weight used in case of fiting...

    Wfree = diags(w[free].flatten())
    xy = xyz[:, :2]

    # If Rols are set up

    rol_x = []
    rol_y = []

    for key in form.vertices_where({'rol_x': True}):
        rol_x.append(k_i[key])

    for key in form.vertices_where({'rol_y': True}):
        rol_y.append(k_i[key])

    free_x = list(set(free) - set(rol_x))
    free_y = list(set(free) - set(rol_y))

    # C and E matrices

    C = th.tensor(connectivity_matrix(edges, 'array'))
    Ci = C[:, free]
    Cf = C[:, fixed]
    Cftx = C[:, rol_x].transpose
    Cfty = C[:, rol_y].transpose
    Ct = C.transpose
    Cit = Ci.transpose
    Citx = Ct[free_x, :]  # .toarray()
    City = Ct[free_y, :]  # .toarray()
    uvw = C.dot(xyz)
    U = diags(uvw[:, 0].flatten())
    V = diags(uvw[:, 1].flatten())
    E = svstack((Citx.dot(U), City.dot(V)), dtype='csr').toarray()
    print('Equilibrium Matrix Shape: ', E.shape)

    start_time = time.time()

    # Independent and dependent branches

    if find_inds:

        if indset:
            ind = []
            for u, v in form.edges_where({'_is_edge': True}):
                if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
                    ind.append(uv_i[(u, v)])
        else:
            ind = find_independents(E)

        k = len(ind)
        dep = list(set(range(m)) - set(ind))

        elapsed_time = time.time() - start_time

        if printout:
            print('Shape Equilibrium Matrix: ', E.shape)
            print('Found {0} independents'.format(k))
            print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

        for u, v in form.edges_where({'_is_edge': True}):
            form.edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

        Edinv = -csr_matrix(pinv(E[:, dep]))
        Ei = E[:, ind]

    else:
        k = m
        ind = list(set(range(m)))
        dep = []
        Edinv = []
        Ei = []

    # Set-up

    lh = normrow(C.dot(xy))**2
    p = vstack([px[free_x], py[free_y]])
    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]

    if any(p):
        check_hor = check_horizontal(E, p)
        if check_hor:
            print('Horizontal Loads can be taken!')
        else:
            print('Horizontal Loads are not suitable for this FD!')

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k,
            lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty)
    args_inds = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Citx, City, Cf, U, V, p, px, py, pz, z, free_x, free_y, fixed, lh, sym, k)
    if find_inds is True:
        checked = check_independents(args_inds, tol=tol)
        if checked:
            print('Independents checked!')
            pass
        else:
            print('Warning: independent edges not equilibrated')

    return args


def initialise_problem_general(form, printout=None, tol=0.001):
    """ Initialise the problem for a given Form-Diagram building the main matrices used in the subsequent analysis.

    Parameters
    ----------
    form : FormDiagram
        The FormDiagram.
    printout : bool
        Control output.
    find_inds : bool
        If True will calculate the independents
    indset : list
        Independent set to use. If empty the independents are calculated normally.
    tol: float
        Tolerance for check the independent edges equilibrium

    Returns
    -------
    args
        List of matrices and vectors used to perform the optimisation.

    """

    # Mapping

    k_i = form.key_index()
    uv_i = form.uv_index()
    i_uv = form.index_uv()

    # Vertices and edges

    n = form.number_of_vertices()
    m = len(list(form.edges_where({'_is_edge': True})))
    fixed = [k_i[key] for key in form.fixed()]
    nb = len(fixed)
    edges = [(k_i[u], k_i[v]) for u, v in form.edges_where({'_is_edge': True})]
    free = list(set(range(n)) - set(fixed))

    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})]).reshape(-1, 1)  # review need of 'is_edge': True
    qmax = array([form.edge_attribute((u, v), 'qmax') for u, v in form.edges_where({'_is_edge': True})]).reshape(-1, 1)
    qmin = array([form.edge_attribute((u, v), 'qmin') for u, v in form.edges_where({'_is_edge': True})]).reshape(-1, 1)

    # Co-ordinates, loads and constraints

    xyz = zeros((n, 3))
    x = zeros((n, 1))
    y = zeros((n, 1))
    z = zeros((n, 1))
    s = zeros((n, 1))
    px = zeros((n, 1))
    py = zeros((n, 1))
    pz = zeros((n, 1))
    lb = zeros((n, 1))
    ub = zeros((n, 1))
    xlimits = zeros((n, 2))
    ylimits = zeros((n, 2))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i] = vertex.get('x')
        y[i] = vertex.get('y')
        z[i] = vertex.get('z')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i] = vertex.get('target', 0)
        xlimits[i, 0] = vertex.get('xmin', None)
        xlimits[i, 1] = vertex.get('xmax', None)
        ylimits[i, 0] = vertex.get('ymin', None)
        ylimits[i, 1] = vertex.get('ymax', None)
        lb[i] = vertex.get('lb', None)
        ub[i] = vertex.get('ub', None)

    # Partial supports, or rollers.

    rol_x = []
    rol_y = []

    for key in form.vertices_where({'rol_x': True}):
        rol_x.append(k_i[key])
    for key in form.vertices_where({'rol_y': True}):
        rol_y.append(k_i[key])

    free_x = list(set(free) - set(rol_x))
    free_y = list(set(free) - set(rol_y))

    # C and E matrices

    C = connectivity_matrix(edges, 'csr')
    Ci = C[:, free]
    Cb = C[:, fixed]
    Cbtx = C[:, rol_x].transpose()
    Cbty = C[:, rol_y].transpose()
    Ct = C.transpose()
    Cit = Ci.transpose()
    Citx = Ct[free_x, :]
    City = Ct[free_y, :]
    uvw = C.dot(xyz)
    U = diags(uvw[:, 0].flatten())
    V = diags(uvw[:, 1].flatten())
    E = svstack((Citx.dot(U), City.dot(V))).toarray()

    # Settings of a free to move problem 'all q are variable'

    ind = list(range(m))
    dep = []
    k = m
    B = identity(m)

    # Create Class and build matrices

    problem = Problem()

    problem.q = q
    problem.m = m
    problem.n = n
    problem.nb = nb
    problem.E = E
    problem.C = C
    problem.Ct = Ct
    problem.Ci = Ci
    problem.Cit = Cit
    problem.Cb = Cb
    problem.U = U
    problem.V = V
    problem.P = hstack([px, py, pz])
    problem.free = free
    problem.fixed = fixed
    problem.phfree = hstack([px[free], py[free]])  # deprecated?
    problem.ph = vstack([px[free], py[free]])
    problem.lb = lb
    problem.ub = ub
    problem.lb0 = lb
    problem.ub0 = ub
    problem.s = s
    problem.X = hstack([x, y, z])
    problem.x0 = x
    problem.y0 = y
    problem.free_x = free_x
    problem.free_y = free_y
    problem.rol_x = rol_x
    problem.rol_y = rol_y
    problem.Citx = Citx
    problem.City = City
    problem.Cbtx = Cbtx
    problem.Cbty = Cbty
    problem.xlimits = xlimits
    problem.ylimits = ylimits
    problem.qmin = qmin
    problem.qmax = qmax
    problem.k_i = k_i
    problem.uv_i = uv_i
    problem.i_uv = i_uv
    problem.ind = ind
    problem.k = k
    problem.dep = dep
    problem.B = B

    return problem


def adapt_problem_to_fixed_diagram(problem, form, printout=False):
    """ Adapt the problem assuming that the form diagram is fixed in plan."""

    ind = []

    start_time = time.time()

    # Independent and dependent branches

    if form.attributes['indset']:
        ind = []
        for u, v in form.edges_where({'_is_edge': True}):
            if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in form.attributes['indset']:
                ind.append(problem.uv_i[(u, v)])
    else:
        ind = find_independents(problem.E)  # see if it can be improved with crs matrix

    k = len(ind)
    dep = list(set(range(problem.m)) - set(ind))

    elapsed_time = time.time() - start_time

    if printout:
        print('Reduced problem to {0} force variables'.format(k))
        print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    gkeys = []
    for u, v in form.edges_where({'_is_edge': True}):
        if problem.uv_i[(u, v)] in ind:
            form.edge_attribute((u, v), 'is_ind', True)
            gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
        else:
            form.edge_attribute((u, v), 'is_ind', False)
    form.attributes['indset'] = gkeys

    Edinv = -csr_matrix(pinv(problem.E[:, dep]))
    Ei = csr_matrix(problem.E[:, ind])
    B = zeros((problem.m, k))
    B[dep] = Edinv.dot(Ei).toarray()
    B[ind] = identity(k)

    problem.ind = ind
    problem.k = k
    problem.dep = dep
    problem.B = B
    problem.Edinv = Edinv
    problem.Ei = Ei

    return


def adapt_problem_to_sym_diagram(problem, form, list_axis_symmetry=None, center=None, correct_loads=True, printout=False):
    """ Adapt the problem assuming that the form diagram is symmetric."""

    start_time = time.time()

    apply_sym_to_form(form, list_axis_symmetry, center, correct_loads)

    Esym = build_symmetry_transformation(form, printout=False)
    mapsym = form.build_symmetry_map()
    ind = sorted(list(mapsym.values()))

    k = len(ind)
    dep = list(set(range(problem.m)) - set(ind))

    elapsed_time = time.time() - start_time

    if printout:
        print('Reduced problem to {0} force variables'.format(k))
        print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    for u, v in form.edges_where({'_is_edge': True}):
        form.edge_attribute((u, v), 'is_ind', True if problem.uv_i[(u, v)] in ind else False)

    B = Esym

    problem.ind = ind
    problem.k = k
    problem.dep = dep
    problem.B = B

    return


def adapt_problem_to_sym_and_fixed_diagram(problem, form, list_axis_symmetry=None, center=None, correct_loads=True, printout=False):
    """ Adapt the problem assuming that the form diagram is symmetric and fixed in plane."""

    start_time = time.time()

    adapt_problem_to_fixed_diagram(problem, form, printout)

    apply_sym_to_form(form, list_axis_symmetry, center, correct_loads)

    ind = problem.ind
    k = problem.k
    i_uv = problem.i_uv
    uv_i = problem.uv_i
    unique_sym = []
    ind_reduc = []

    for index in ind:
        u, v = i_uv[index]
        index_sym = form.edge_attribute((u, v), 'sym_key')
        if index_sym not in unique_sym:
            unique_sym.append(index_sym)
            ind_reduc.append(index)

    k_sym = len(ind_reduc)
    Bsym = zeros((k, k_sym))

    j = 0
    for key in unique_sym:
        edges_sym = [uv_i[(u, v)] for u, v in form.edges_where({'sym_key': key})]
        col = zeros((k, 1))
        for i, ind_i in enumerate(ind):
            if ind_i in edges_sym:
                col[i] = 1.0
        Bsym[:, [j]] = col
        j += 1

    B = problem.B.dot(Bsym)
    ind = ind_reduc
    k = k_sym
    dep = list(set(range(problem.m)) - set(ind))

    elapsed_time = time.time() - start_time

    if printout:
        print('Reduced problem to {0} force variables'.format(k))
        print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    for u, v in form.edges_where({'_is_edge': True}):
        form.edge_attribute((u, v), 'is_ind', True if problem.uv_i[(u, v)] in ind else False)

    problem.ind = ind
    problem.k = k
    problem.dep = dep
    problem.B = B

    return


def apply_sym_to_form(form, list_axis_symmetry=None, center=None, correct_loads=True):
    """ Apply symmetry to the form diagram."""

    if not list_axis_symmetry:
        data = form.parameters
        if data['type'] in ['radial_fd', 'radial_spaced_fd', 'spiral_fd']:
            apply_radial_symmetry(form, center=center, correct_loads=correct_loads)
        elif data['type'] in ['cross_fd', 'fan_fd', 'cross_diagonal', 'cross_with_diagonal', 'ortho_fd']:
            list_axis_symmetry = find_sym_axis_in_rect_patterns(data)
            print('Axis of Symmetry identified:', list_axis_symmetry)
        else:
            print('Symmetry applied to a unknown form diagram type')
            raise NotImplementedError

    if list_axis_symmetry:
        apply_symmetry_from_axis(form, list_axis_symmetry=list_axis_symmetry, correct_loads=correct_loads)

    return
