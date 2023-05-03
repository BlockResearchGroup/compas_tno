from numpy import array
from numpy import zeros
from numpy import vstack
from numpy import hstack
from numpy import identity
from numpy import asarray
from numpy.linalg import pinv
from numpy.linalg import svd

from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import vstack as svstack

from compas.numerical import connectivity_matrix

from compas.utilities import reverse_geometric_key

from compas.geometry import Point
from compas.geometry import distance_point_point_xy

from compas_tno.algorithms import check_independents
from compas_tno.algorithms import check_horizontal_loads
from compas_tno.algorithms import find_independents

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

    None

    Attributes
    ----------

    q : array(m x 1)
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

        self.q = None
        self.m = None
        self.n = None
        self.ni = None
        self.nb = None
        self.E = None
        self.C = None
        self.Ct = None
        self.Ci = None
        self.Cit = None
        self.Cb = None
        self.U = None
        self.V = None
        self.P = None
        self.free = None
        self.fixed = None
        self.ph = None
        self.lb = None
        self.ub = None
        self.lb0 = None
        self.ub0 = None
        self.s = None
        self.X = None
        self.x0 = None
        self.y0 = None
        self.free_x = None
        self.free_y = None
        self.rol_x = None
        self.rol_y = None
        self.Citx = None
        self.City = None
        self.Cbtx = None
        self.Cbty = None
        self.xlimits = None
        self.ylimits = None
        self.qmin = None
        self.qmax = None
        self.k_i = None
        self.uv_i = None
        self.i_uv = None
        self.ind = None
        self.k = None
        self.dep = None
        self.B = None

        pass

    @classmethod
    def from_formdiagram(cls, form, indset=None, printout=False, find_inds=False):
        """Initialise the problem object form a form diagram."""
        problem = cls()

        return problem


def initialise_form(form, find_inds=True, method='SVD', printout=False, tol=None):
    """ Initialise the problem for a Form-Diagram and return the FormDiagram with independent edges assigned and the matrices relevant to the equilibrium problem.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The FormDiagram
    find_inds : bool, optional
        Whether or not independents should be found (fixed diagram), by default True
    method : str, optional
        Method to find independent edges, the default is 'SVD'. More options to come.
    printout : bool, optional
        Whether or not prints should appear on the screen, by default False
    tool : float, optional
        Tolerance od independents, by default None

    Returns
    -------
    M : :class:`~compas_tno.problems.Problem`
        The object with all the matrices essential to the analysis.

    Note
    -------
    The FormDiagram is updated in place. Check ``initialise_problem_general`` for more info.

    """

    M = initialise_problem_general(form)

    if find_inds:
        adapt_problem_to_fixed_diagram(M, form, method=method, printout=printout, tol=tol)

    return M


def initialise_problem_general(form):
    """ Initialise the problem for a given Form-Diagram building the main matrices used in the subsequent analysis.

    Parameters
    ----------
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The FormDiagram.

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
    ni = len(free)

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
        if abs(s[i]) < 1e-6:
            s[i] = 0.0
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
    d = zeros((m, 1))

    # Permutation matrix to fixed and free keys
    # with this matrix z = [Pmatrix]z will have the first ni nodes as free
    # and the last nb nodes as fixed

    Pfree = zeros((n, ni))
    Pfixed = zeros((n, nb))
    Pfree[free, :] = identity(ni)
    Pfixed[fixed, :] = identity(nb)
    Pmatrix = hstack([Pfree, Pfixed])

    # Create Class and build matrices

    problem = Problem()

    problem.q = q
    problem.m = m
    problem.n = n
    problem.ni = ni
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
    problem.ph = vstack([px[free_x], py[free_y]])
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
    problem.d = d
    problem.Pmatrix = Pmatrix
    # problem.Bfixed = Bfixed

    return problem


def adapt_problem_to_fixed_diagram(problem, form, method='SVD', printout=False, tol=None):
    """Adapt the problem assuming that the form diagram is fixed in plan.

    Parameters
    ----------
    problem : :class:`~compas_tno.problems.Problem`
        Matrices of the problem
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to be analysed
    method : str, optional
        Method to find independent edges, the default is 'SVD'. More options to come.
    printout : bool, optional
        If prints should show in the screen, by default False
    tol : float, optional
        Tolerance of the singular values, by default None

    """

    ind = []
    tol_old_ind = 1e-3

    start_time = time.time()

    # Independent and dependent branches

    if form.attributes['indset']:
        # check if it is a string and "restaure the points"
        indset = [a if not isinstance(a, str) else reverse_geometric_key(a) for a in form.attributes['indset']]
        ind = []

        for u, v in form.edges_where({'_is_edge': True}):
            index = problem.uv_i[(u, v)]
            edgemid = Point(*(form.edge_midpoint(u, v)[:2] + [0]))
            for pt in indset:
                if distance_point_point_xy(edgemid, pt) < tol_old_ind:
                    ind.append(index)
                    break
            if index in ind:
                indset.remove(pt)

        if printout:
            print('Loaded {} previous independents'.format(len(form.attributes['indset'])))
            print('Found {} independents in the new pattern'.format(len(ind)))
        if len(form.attributes['indset']) != len(ind):
            print('Did not match problem inds')
            ind = find_independents(problem.E, method=method, tol=tol)
    else:
        ind = find_independents(problem.E, method=method, tol=tol)

    k = len(ind)
    dep = list(set(range(problem.m)) - set(ind))

    elapsed_time = time.time() - start_time

    if printout:
        print('Reduced problem to {0} force variables with ind. edges'.format(k))
        print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    points = []
    for u, v in form.edges_where({'_is_edge': True}):
        if problem.uv_i[(u, v)] in ind:
            form.edge_attribute((u, v), 'is_ind', True)
            points.append(Point(*(form.edge_midpoint(u, v)[:2] + [0])))
        else:
            form.edge_attribute((u, v), 'is_ind', False)
    form.attributes['indset'] = points

    rcond = 1e-17
    if tol:
        rcond = tol
    Ed = problem.E[:, dep]
    Edinv = -csr_matrix(pinv(problem.E[:, dep], rcond=rcond))
    Ei = csr_matrix(problem.E[:, ind])
    B = zeros((problem.m, k))
    B[dep] = Edinv.dot(Ei).toarray()
    B[ind] = identity(k)

    d = zeros((problem.m, 1))
    d[dep] = -Edinv.dot(problem.ph)  # q = Bqi + d | d = Ed(-1)*ph

    if any(problem.ph):
        check_hor = check_horizontal_loads(problem.E, problem.ph)
        if check_hor:
            print('Horizontal Loads can be taken!')
        else:
            print('Horizontal Loads are not suitable for this FD!')

    problem.ind = ind
    problem.k = k
    problem.dep = dep
    problem.B = B
    problem.Edinv = Edinv
    problem.Ed = Ed
    problem.Ei = Ei
    problem.d = d
    problem.d0 = d
    problem.time_inds = elapsed_time

    return


def adapt_problem_to_sym_diagram(problem, form, list_axis_symmetry=None, center=None, correct_loads=True, printout=False):
    """Adapt the problem assuming that the form diagram is symmetric.

    Parameters
    ----------
    problem : :class:`~compas_tno.problems.Problem`
        The problem with matrices for calculation
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to analyse
    list_axis_symmetry : [list], optional
        List of the axis of symmetry to consider, by default None
    center : [list], optional
        The center of the pattern if it's a circula pattern, by default None
    correct_loads : bool, optional
        If update should be done in the applied loads regarding the symmetry, by default True
    printout : bool, optional
        If prints should show in the screen, by default False

    """

    start_time = time.time()

    apply_sym_to_form(form, list_axis_symmetry, center, correct_loads)

    Esym = build_symmetry_transformation(form, printout=False)
    mapsym = form.build_symmetry_map()
    ind = sorted(list(mapsym.values()))

    k = len(ind)
    dep = list(set(range(problem.m)) - set(ind))

    elapsed_time = time.time() - start_time

    if printout:
        print('Reduced problem to {0} force variables by SYM'.format(k))
        print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    for u, v in form.edges_where({'_is_edge': True}):
        form.edge_attribute((u, v), 'is_ind', True if problem.uv_i[(u, v)] in ind else False)

    B = Esym

    problem.ind = ind
    problem.k = k
    problem.dep = dep
    problem.B = B

    return


def adapt_problem_to_sym_and_fixed_diagram(problem, form, method='SVD', list_axis_symmetry=None, center=None, correct_loads=True, printout=False, tol=None):
    """ Adapt the problem assuming that the form diagram is symmetric and fixed in plane.

    Parameters
    ----------
    problem : :class:`~compas_tno.problems.Problem`
        The problem with matrices for calculation
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to analyse
    method : str, optional
        Method to find independent edges, the default is 'SVD'. More options to come.
    list_axis_symmetry : [list], optional
        List of the axis of symmetry to consider, by default None
    center : [list], optional
        The center of the pattern if it's a circula pattern, by default None
    correct_loads : bool, optional
        If update should be done in the applied loads regarding the symmetry, by default True
    printout : bool, optional
        If prints should show in the screen, by default False

    """

    start_time = time.time()

    adapt_problem_to_fixed_diagram(problem, form, method=method, printout=printout, tol=tol)

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
    """ Apply symmetry to the form diagram.

    Parameters
    ----------
    problem : :class:`~compas_tno.problems.Problem`
        The problem with matrices for calculation
    form : :class:`~compas_tno.diagrams.FormDiagram`
        The form diagram to analyse
    list_axis_symmetry : [list], optional
        List of the axis of symmetry to consider, by default None
    center : [list], optional
        The center of the pattern if it's a circula pattern, by default None
    correct_loads : bool, optional
        If update should be done in the applied loads regarding the symmetry, by default True

    """

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


def check_bad_independents(problem, form, printout=False, eps=10e+6):

    keep_inds = []

    for i_ind in range(len(problem.ind)):
        column_i = abs(problem.B[:, i_ind]).flatten()
        max_col = max(column_i)
        if max_col < eps:
            keep_inds.append(i_ind)

    if len(keep_inds) == problem.m:
        return
    else:
        print('original inds:', len(problem.ind), problem.ind)
        print('Updating the relevant independent')
        inds = []
        for i_ind in keep_inds:
            inds.append(problem.ind[i_ind])

        q0 = zeros((problem.m, 1))
        form.edges_attribute('is_ind', False)
        for i, edge in enumerate(form.edges_where({'_is_edge': True})):
            if i in inds:
                form.edge_attribute(edge, 'is_ind', True)
            elif i in problem.ind:
                q0[i] = form.edge_attribute(edge, 'q')

        print(q0)

        problem.ind = inds
        problem.dep = list(set(range(problem.m)) - set(inds))
        problem.B = problem.B[:, keep_inds]
        problem.d += q0
        problem.k = len(inds)
        print('final inds', len(problem.ind), problem.ind)

    return


def plot_svds(M, tol=None):
    """ Plot a diagram of SVD.

    Parameters
    ----------
    M : :class:`~compas_tno.problems.Problem`
        The problem with matrices for calculation
    tol : float (optional)
        The tolerance to find the independents
    """

    import matplotlib.pyplot as plt
    from numpy.linalg import matrix_rank

    k = len(M.ind)
    n, m = M.E.shape
    # mn_min = min(n, m)

    _, s, _ = svd(asarray(M.E))
    # _, s, _ = svds(M.E, k=min(M.E.shape), solver='propack')
    print('max/min singular vectors E', max(s), min(s), len(s))
    print('Shape E: {} | #ind: {} | rank : {}:'.format(M.E.shape, k, matrix_rank(M.E, tol=tol)))

    # mn = max(m - n, n - m)
    mn = max(m - n, 0)
    zs = k - mn
    print('# null-SVD:', zs)

    # if zs + 3 > len(s):
    #     print('Last {} SVs: {}'.format(len(s[-(zs + 3):]), s[-(zs + 3):]))
    # print('ALL SVS:', s)

    check = check_independents(M)
    print('Check independents:', check)
    print('Shape/ Rank Ed:', M.Ed.shape, matrix_rank(M.Ed))

    if zs > 0:
        first_zero = s[len(s)-zs]
        last_non_zero = s[len(s)-zs - 1]

        x_zs = list(range(len(s)))[-zs:]
        y_zs = s[-zs:]

        print('Last Non Zero:', last_non_zero)
        print('First Zero SV:', first_zero)

        # len_zero = len(s[s<1.0])
        lin_x = len(s) - zs

        porcentage_key = (last_non_zero - first_zero)/last_non_zero
        print('percentage key is:', porcentage_key)

        _, ax = plt.subplots()
        ax.plot(s)
        ax.scatter(x_zs, y_zs)
        ax.plot([lin_x, lin_x], [0, 1.0], color='black')
        ax.plot()
        plt.show()

    else:
        print('EQ. Matrix is FULL rank: (k = m-n)')
        _, ax = plt.subplots()
        # ax.set_ylim(2.7, 3.0)
        ax.plot(s)
        plt.show()

    return
