
from numpy import vstack
from numpy import hstack
from numpy import array
from numpy.linalg import matrix_rank
from numpy.random import rand
from math import sqrt


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'find_independents',
    'independents_exclude',
    'independents_include',
    'inds_incl_excl',
    'check_independents',
]

def find_independents(E):

    _, m = E.shape
    Etemp = E[:,[0]]
    ind = []

    for i in range(1,m):
        Etest = hstack([Etemp,E[:,[i]]])
        _ , ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind

def independents_exclude(E, outs):

    _, m = E.shape
    possible = list(set(range(m)) - set(outs))
    Eouts = E[:,outs]
    if matrix_rank(Eouts) == Eouts.shape[1]:
        Etemp = Eouts
        ind = []
    else:
        print('Warning, could not exclude all')
        return find_independents(E)

    for i in possible:
        Etest = hstack([Etemp,E[:,[i]]])
        _ , ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind

def independents_include(E, ins):

    n, m = E.shape
    if len(ins) > (m-n):
        print('Too many included edges - limit to number of independents: {0}'.format((m-n)))
        ins = ins[:(m-n)]
    not_in = list(set(range(m)) - set(ins))
    Ein = E[:,ins]
    while matrix_rank(Ein) < Ein.shape[1]:
        print('Warning, edges are dependent among each other')
        ins = ins[matrix_rank(Ein)]
    if len(ins) == (m-n):
        Enot_in = E[:,not_in]
        if matrix_rank(Enot_in) == Enot_in.shape[1]:
            return ins
        else:
            print('Warning, edges do not form an independent set')
            return find_independents(E)
    ind = ins
    Etemp = E[:,[not_in[0]]]
    Etemp.shape

    for i in not_in[1:]:
        Etest = hstack([Etemp,E[:,[i]]])
        _ , ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind

def inds_incl_excl(E, ins, outs):

    n, m = E.shape
    if len(ins) > (m-n):
        print('Too many included edges - limit to number of independents: {0}'.format((m-n)))
        ins = ins[:(m-n)]
    not_in = list(set(range(m)) - set(ins))
    possible = list(set(not_in)-set(outs))
    Ein = E[:,ins]
    Eouts = E[:,outs]
    while matrix_rank(Ein) < Ein.shape[1]:
        print('Warning, included edges are dependent among each other')
        ins = ins[matrix_rank(Ein)]
    if len(ins) == (m-n):
        Enot_in = E[:,not_in]
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
        Etemp = E[:,[not_in[0]]]
        possible = not_in[1:]
    
    for i in possible:
        Etest = hstack([Etemp,E[:,[i]]])
        _ , ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind


def check_independents(args, tol = 0.001):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y = args
    checked = True
    if tol > 0:
        for i in range(10**3):
            q[ind, 0] = rand(k) * 10
            q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
            Rx = array(Cit.dot(U * q.ravel()) - px[free].ravel())
            Ry = array(Cit.dot(V * q.ravel()) - py[free].ravel())
            R  = sqrt(max(Rx**2 + Ry**2))
            if R > tol:
                checked = False
                break

    return checked