
from numpy import hstack
from numpy import array
from numpy import any
from numpy.linalg import matrix_rank
from numpy.random import rand
from math import sqrt
from copy import deepcopy


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'find_independents',
    'independents_exclude',
    'independents_include',
    'inds_incl_excl',
    'check_independents',
    'check_horizontal'
]


def find_independents(E):
    """ Find inndependent edges in a given equilibrium matrix.

    Parameters
    ----------
    array : E
        Equilibrium matrix.

    Returns
    -------
    list : ind
        Independent columns.

    """

    ind = []
    _, m = E.shape

    # Take as independent all edges connecting supported points.

    for i in range(m):
        if any(E[:, i]):
            pass
        else:
            ind.append(i)

    internal = list(set(range(m)) - set(ind))
    Etemp = E[:, [internal[0]]]

    # Proceed for the internal edges only.

    for i in internal[1:]:
        Etest = hstack([Etemp, E[:, [i]]])
        _, ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind


def independents_exclude(E, outs):
    """ Find inndependent edges with certain to exclude.

    Parameters
    ----------
    array : E
        Equilibrium matrix.

    list : outs
        List to exclude.

    Returns
    -------
    list : ind
        Independent columns.

    """

    _, m = E.shape
    possible = list(set(range(m)) - set(outs))
    Eouts = E[:, outs]
    if matrix_rank(Eouts) == Eouts.shape[1]:
        Etemp = Eouts
        ind = []
    else:
        print('Warning, could not exclude all')
        return find_independents(E)

    for i in possible:
        Etest = hstack([Etemp, E[:, [i]]])
        _, ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind


def independents_include(E, ins):
    """ Find inndependent edges with certain to include.

    Parameters
    ----------
    array : E
        Equilibrium matrix.

    list : ins
        List to include.

    Returns
    -------
    list : ind
        Independent columns.

    """

    n, m = E.shape
    if len(ins) > (m-n):
        print('Too many included edges - limit to number of independents: {0}'.format((m-n)))
        ins = ins[:(m-n)]
    not_in = list(set(range(m)) - set(ins))
    Ein = E[:, ins]
    while matrix_rank(Ein) < Ein.shape[1]:
        print('Warning, edges are dependent among each other')
        ins = ins[matrix_rank(Ein)]
    if len(ins) == (m-n):
        Enot_in = E[:, not_in]
        if matrix_rank(Enot_in) == Enot_in.shape[1]:
            return ins
        else:
            print('Warning, edges do not form an independent set')
            return find_independents(E)
    ind = ins
    Etemp = E[:, [not_in[0]]]
    Etemp.shape

    for i in not_in[1:]:
        Etest = hstack([Etemp, E[:, [i]]])
        _, ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind


def inds_incl_excl(E, ins, outs):
    """ Find inndependent edges with certain to exclude and include.

    Parameters
    ----------
    array : E
        Equilibrium matrix.

    list : outs
        List to exclude.

    list : ins
        List to include.

    Returns
    -------
    list : ind
        Independent columns.

    """

    n, m = E.shape
    if len(ins) > (m-n):
        print('Too many included edges - limit to number of independents: {0}'.format((m-n)))
        ins = ins[:(m-n)]
    not_in = list(set(range(m)) - set(ins))
    possible = list(set(not_in)-set(outs))
    Ein = E[:, ins]
    Eouts = E[:, outs]
    while matrix_rank(Ein) < Ein.shape[1]:
        print('Warning, included edges are dependent among each other')
        ins = ins[matrix_rank(Ein)]
    if len(ins) == (m-n):
        Enot_in = E[:, not_in]
        if matrix_rank(Enot_in) == Enot_in.shape[1]:
            return ins
        else:
            print('Warning, edges do not form an independent set')
            return find_independents(E)
    if matrix_rank(Eouts) == Eouts.shape[1]:
        Etemp = Eouts
        ind = ins
    else:
        print('Warning, could not exclude all')
        ind = ins
        Etemp = E[:, [not_in[0]]]
        possible = not_in[1:]

    for i in possible:
        Etest = hstack([Etemp, E[:, [i]]])
        _, ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind


def check_independents(args_inds, tol=0.001):
    """ Run checks to verify the independennts.

    Parameters
    ----------
    list : args_inds
        Arguments from the optimisation.

    tol : float (0.001)
        Allowed error.

    Returns
    -------
    bool : checked
        True if independents are checked.

    """

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Citx, City, Cf, U, V, p, px, py, pz, z, free_x, free_y, fixed, lh, sym, k = args_inds
    checked = True
    q_ = deepcopy(q)
    if tol > 0:
        for i in range(10**3):
            q_[ind, 0] = rand(k) * 10
            q_[dep] = Edinv.dot(- p + Ei.dot(q_[ind]))
            Rx = array(Citx.dot(U * q_.ravel()) - px[free_x].ravel())
            Ry = array(City.dot(V * q_.ravel()) - py[free_y].ravel())
            R = max(sqrt(max(Rx**2)), sqrt(max(Ry**2)))
            if R > tol:
                checked = False
                break

    return checked


def check_horizontal(E, p):
    """ Run checks to verify of a certain horizontal distribution can be taken by the fixed form diagram.

    Parameters
    ----------
    array : E
        Equilibrium matrix.
    array : p
        Vector of horizontal loads.
    Returns
    -------
    bool : checked
        True if loads can be taken.

    """

    r = matrix_rank(E)
    r_ = matrix_rank(hstack([E, p]))

    if r == r_:
        checked = True
    else:
        checked = False

    return checked
