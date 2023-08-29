
from numpy import hstack
from numpy import array
from numpy import asarray
from numpy import any
from numpy.linalg import matrix_rank
from numpy.linalg import svd
from numpy.random import rand
from numpy import diag
from scipy.linalg import qr
from math import sqrt


def find_independents_forward(E, nind=None, tol=None):
    """ Find independent edges of the matrix E with the forward method.
    The matrix E is reconstructed column by column and the rank is computed at each step.
    Everytime that a column is added and the rank does not increase, the edge corresp. that column is selected as independent.

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    nind : ind
        Number of independent to stop the process, if known in advance, by default None.
    tol : float
        Tolerance for small Singular Values. Default is None.

    Returns
    -------
    ind : list
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
        if matrix_rank(Etest, tol=tol) < ncol:  # testing tol
            ind.append(i)
        else:
            Etemp = Etest
        if nind:
            if len(ind) == nind:
                break

    return ind


def find_independents_backward(E, nind=None, tol=None):
    """ Find independent edges of the matrix E with the backward method.
    The last columns of the matrix are removed sequentially and the rank is checked.
    Everytime thae rank does not decrease, the edge is selected as independent.

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    nind : ind
        Number of independent to stop the process, if known in advance, by default None.
    tol : float
        Tolerance for small singular values. Default is None.

    Returns
    -------
    ind : list
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
    rank = matrix_rank(E, tol=tol)

    for i in range(m):
        i_inv = m - i - 1
        indices_red = [j for j in internal if not j == i_inv]
        Ered = E[:, indices_red]
        ri = matrix_rank(Ered, tol=tol)
        if ri == rank:  # belongs to the nullspace
            ind.append(i_inv)
            internal.remove(i_inv)
            rank = ri
        if nind:
            if len(ind) == nind:
                break

    return ind


def find_independents_QR(E, tol=None):
    """ Find independent edges of the matrix E based on a permuted QR factorization.
    The matrix E [m x n] (n > m) is factored in Q, R, P matrices, where Q [m x m] is a square matrix, R [m x n] an upper triangular
    matrix and P is the column permutation applied.
    This method is faster, but it does not preserve the numbering structure of the form diagram.

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    tol : float
        Tolerance for small singular values. Default is None.

    Returns
    -------
    ind : list
        Independent columns.

    """

    if not tol:
        tol = 1e-12

    _, R, P = qr(E, pivoting=True)

    dr = abs(diag(R))
    rank_ = sum(dr > tol*dr[0])
    ind = sorted(list(P[rank_:]))

    return ind


def find_independents(E, method='SVD', tol=None):
    """ Overall method to find independent edges dependent on method

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    method : str
        Method to find independents. Default is 'SVD'. More to come.
    tol : float
        Tolerance for small singular values. Default is None.

    Returns
    -------
    ind : list
        Independent columns.

    """

    if method == 'SVD':
        ind = find_independents_forward(E, tol)
    elif method == 'QR':
        ind = find_independents_QR(E, tol)
    else:
        raise ValueError('Plese select a valid method to find the independent edges')

    return ind


def independents_exclude(E, outs):
    """ Find inndependent edges with certain to exclude.

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    outs : list
        List of edges (columns in the matrix) to exclude.

    Returns
    -------
    ind : list
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
        return find_independents_forward(E)

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
    E : array
        Equilibrium matrix.
    ins : list
        List of edges (columns in the matrix) to include.

    Returns
    -------
    ind : list
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
            return find_independents_forward(E)
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
    E : array
        Equilibrium matrix.
    ins : list
        List of edges (columns in the matrix) to include.
    outs : list
        List of edges (columns in the matrix) to exclude.

    Returns
    -------
    ind : list
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
            return find_independents_forward(E)
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


def check_independents(M, tol=0.001):
    """ Run checks to verify the independennts.

    Parameters
    ----------
    M : :class:`~compas_tno.problems.Problem`
        Arguments from the optimisation.
    tol : float, optional
        Allowed error.
        The default values is ``0.001``.

    Returns
    -------
    checked : bool
        True if independents are checked.

    """

    checked = True
    if tol > 0:
        for i in range(10**2):
            qid = array(rand(M.k) * 10).reshape(-1, 1)
            q_ = M.B @ qid + M.d
            res_x = M.Citx.dot(M.U * q_.ravel()) - M.ph[:len(M.free_x)].ravel()
            res_y = M.City.dot(M.V * q_.ravel()) - M.ph[:len(M.free_y)].ravel()
            R = max(sqrt(max(res_x**2)), sqrt(max(res_y**2)))
            if R > tol:
                checked = False
                print('deviation exceeded limit:', R)
                break

    # checked = True
    # q_ = deepcopy(M.q)
    # if tol > 0:
    #     for i in range(10**3):
    #         q_[M.ind, 0] = rand(M.k) * 10
    #         q_[M.dep] = M.Edinv.dot(- M.ph + M.Ei.dot(q_[M.ind]))
    #         Rx = array(M.Citx.dot(M.U * q_.ravel()) - M.ph[:len(M.free_x)].ravel())
    #         Ry = array(M.City.dot(M.V * q_.ravel()) - M.ph[len(M.free_y):].ravel())
    #         R = max(sqrt(max(Rx**2)), sqrt(max(Ry**2)))
    #         # print(R)
    #         if R > tol:
    #             checked = False
    #             break

    return checked


def check_horizontal_loads(E, p):
    """ Run checks to verify of a certain horizontal distribution can be taken by the fixed form diagram.

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    p : array
        Vector of horizontal loads.

    Returns
    -------
    checked : bool
        True if independents are checked.

    """

    r = matrix_rank(E)
    r_ = matrix_rank(hstack([E, p]))

    if r == r_:
        checked = True
    else:
        checked = False

    return checked


def count_inds(E, tol=None):
    """ Count the number of independent edges for a given singular value tolerance.

    Parameters
    ----------
    E : array
        Equilibrium matrix.
    tol : float
        Tolerance for small Singular Values. Default is None.

    Returns
    -------
    k : int
        number of independent edges
    rd : int
        Rank Deficiency of the Equilibrium Matrix.

    """

    nl, ncol = E.shape
    k = max(nl, ncol) - matrix_rank(E, tol=tol)
    rd = k - abs(nl - ncol)

    return k, rd


def matrix_svd(E):
    """ Return the singular values of the matrix E.

    Parameters
    ----------
    E : array
        Matrix to compute the singular values.

    Returns
    -------
    s : array
        vector with the singular values/
    """

    _, s, _ = svd(asarray(E))

    return s
