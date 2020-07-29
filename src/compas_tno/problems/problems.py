from numpy import array
from numpy import zeros
from numpy import vstack
from numpy import newaxis
from numpy.linalg import pinv
from numpy.linalg import matrix_rank

from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse import vstack as svstack

from compas.numerical import connectivity_matrix
from compas.numerical import normrow

from compas.utilities import geometric_key

from compas_tno.algorithms.independents import find_independents
from compas_tno.algorithms.independents import check_independents
from compas_tno.algorithms.independents import check_horizontal

import time


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'initialise_problem',
    'initialise_form'
]


def initialise_problem(form, indset=None, printout=None, find_inds=True, tol=0.001):
    """ Initialise the problem for a given Form-Diagram and return the set of matrices and vectors to optimise.

    Parameters
    ----------
    form : obj
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
    i_k = form.index_key()
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
    for key in form.vertices_where({'is_fixed': False}):
        i = k_i[key]
        lbi = form.vertex_attribute(key, 'lb')
        ubi = form.vertex_attribute(key, 'ub')
        if lb is not None:
            lb.append(lbi)
            lb_ind.append(i)
        if ub is not None:
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
    print('Equilibrium Matrix Shape: ', E.shape)
    print('Rank Equilibrium Matrix: ', matrix_rank(E))

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

    if any(p) == True:
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


def initialise_form(form, indset=None, printout=False, find_inds=True):
    """ Initialise the problem for a Form-Diagram and return the FormDiagram with independent edges assigned

    Parameters
    ----------
    form : obj
        The FormDiagram.

    Returns
    -------
    form : obj
        The FormDiagram.

    """

    i_uv = form.index_uv()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
    # form.attributes['args'] = args
    q, ind, dep = args[:3]

    form.update_default_edge_attributes({'is_ind': False})
    gkeys = []
    for i in ind:
        u, v = i_uv[i]
        gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
        form.edge_attribute((u, v), 'is_ind', value=True)
    form.attributes['indset'] = gkeys

    return form


def initialise_problem_torch(form, indset=None, printout=None, find_inds=True, tol=0.001):
    """ Initialise the problem for a given Form-Diagram and return the set of matrices and vectors to optimise.

    Parameters
    ----------
    form : obj
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

    # import torch as th

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
    E = svstack((Citx.dot(U), City.dot(V))).toarray()
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

    if any(p) == True:
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
