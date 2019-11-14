
from compas_tna.diagrams import FormDiagram

from numpy import array
from numpy import dot
from numpy import zeros
from numpy import vstack
from numpy import hstack
from numpy import newaxis
from numpy.linalg import pinv

from scipy.sparse import csr_matrix
from scipy.sparse import diags

from compas.numerical import connectivity_matrix
from compas.numerical import equilibrium_matrix
from compas.numerical import normrow

from compas.utilities import geometric_key

from compas_thrust.algorithms.independents import find_independents
from compas_thrust.algorithms.independents import check_independents

import time


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'initialise_problem'
]

def initialise_problem(form, indset = None, printout = None, find_inds=True, tol = 0.001):

    """ Initialise the problem for a given Form-Diagram.

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

    k_i  = form.key_index()
    i_k  = form.index_key()
    i_uv = form.index_uv()
    uv_i = form.uv_index()

    # Vertices and edges

    n     = form.number_of_vertices()
    m     = form.number_of_edges()
    fixed = [k_i[key] for key in form.fixed()]
    rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
    sym   = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
    free  = list(set(range(n)) - set(fixed) - set(rol))

    # Constraints

    lb_ind = []
    ub_ind = []
    lb = []
    ub = []
    for key, vertex in form.vertex.items():
        if vertex.get('lb', None):
            lb_ind.append(k_i[key])
            lb.append(vertex['lb'])
        if vertex.get('ub', None):
            ub_ind.append(k_i[key])
            ub.append(vertex['ub'])
    lb = array(lb)
    ub = array(ub)
    lb.shape = (len(lb),1)
    ub.shape = (len(ub),1)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))
    s   = zeros((n, 1))
    w   = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i]  = vertex.get('x')
        y[i]  = vertex.get('y')
        z[i]  = vertex.get('z')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i]  = vertex.get('target', 0)
        w[i] = vertex.get('weight', 1.0) # weight used in case of fiting...

    Wfree = diags(w[free].flatten())
    xy = xyz[:, :2]

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cf  = C[:, fixed]
    Ct = C.transpose()
    Cit = Ci.transpose()
    E   = equilibrium_matrix(C, xy, free, 'csr').toarray()
    uvw = C.dot(xyz)
    U   = uvw[:, 0]
    V   = uvw[:, 1]

    start_time = time.time()

    # Independent and dependent branches

    if find_inds:

        if indset:
            ind = []
            for u, v in form.edges():
                if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
                    ind.append(uv_i[(u, v)])
        else:
            ind = find_independents(E)
        
        k   = len(ind)
        dep = list(set(range(m)) - set(ind))

        elapsed_time = time.time() - start_time

        if printout:
            print('Shape Equilibrium Matrix: ', E.shape)
            print('Found {0} independents'.format(k))
            print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

        for u, v in form.edges():
            form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

        Edinv  = -csr_matrix(pinv(E[:, dep]))
        Ei     = E[:, ind]

    else:
        k = m
        ind = list(set(range(m)))
        dep = []
        Edinv  = []
        Ei = []

    # Set-up

    try:
        t = form.attributes['offset']
    except:
        t = None

    lh     = normrow(C.dot(xy))**2
    p      = vstack([px[free], py[free]])
    q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

    args   = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y)

    if find_independents:
        checked = check_independents(args, tol = tol)
        if checked:
            pass
        else:
            print('Warning: independent edges not equilibrated')
    
    return args



