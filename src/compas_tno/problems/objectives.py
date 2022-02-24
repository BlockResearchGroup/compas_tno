from numpy import sum as npsum

from compas_tno.algorithms import xyz_from_q
from compas_tno.algorithms import q_from_variables
from compas.numerical import normrow

from scipy.sparse import diags


def f_min_thrust(variables, M):
    """Objective function to minimise the horizontal thrust

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k].reshape(-1, 1)
    M.q = q_from_variables(qid, M.B, M.d)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)
    xy = M.X[:, :2]
    P_xy_fixed = M.P[M.fixed][:, :2]

    CfQC = M.Cb.transpose() @ diags(M.q.flatten()) @ M.C
    Rh = CfQC @ xy - P_xy_fixed
    f = sum(normrow(Rh))

    return f


def f_max_thrust(variables, M):
    """Objective function to maximise the horizontal thrust

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    return -1 * f_min_thrust(variables, M)


def f_bestfit(variables, M):
    """Objective function to minimise the vertical squared distance to a given target

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k].reshape(-1, 1)
    M.q = q_from_variables(qid, M.B, M.d)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    f = sum((M.X[:, [2]] - M.s)**2)

    return f


def f_horprojection(variables, M):
    """Objective function to minimise the horizontal squared distance of the nodes on the form diagram to a given pattern

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k].reshape(-1, 1)
    M.q = q_from_variables(qid, M.B, M.d)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    f = sum((M.X[:, [0]] - M.x0)**2) + sum((M.X[:, [1]] - M.y0)**2)

    return f


def f_loadpath_general(variables, M):
    """Objective function to minimise the loadpath

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k].reshape(-1, 1)
    M.q = q_from_variables(qid, M.B, M.d)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    uvw = M.C @ M.X

    l2 = npsum(uvw**2, axis=1).reshape(-1, 1)

    f = (abs(M.q.transpose()).dot(l2)).item()

    return f


def f_complementary_energy(variables, M):
    """Objective function to minimise the complementary energy to a given applied foundation displacement

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k].reshape(-1, 1)
    M.q = q_from_variables(qid, M.B, M.d)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    CfQC = M.Cb.transpose().dot(diags(M.q.flatten())).dot(M.C)
    R = CfQC.dot(M.X) - M.P[M.fixed]
    f = -1 * npsum(R*M.dXb)

    return f


def f_complementary_energy_nonlinear(variables, M):
    """Objective function to minimise the nonlinear complementary energy to a given applied foundation displacement

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    flin = f_complementary_energy(variables, M)

    if M.Ecomp_method == 'simplified':
        fquad = npsum(M.stiff * M.q.reshape(-1, 1) ** 2)  # assuming area and lengths constant - computed in beginning
    if M.Ecomp_method == 'complete':
        fquad = f_loadpath_general(variables, M) * M.stiff

    f = flin + fquad

    return f


def f_max_section(variables, M):
    """Objective function to minimise additional thickness required to find a feasible thrust network

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    nb = len(M.fixed)
    n = M.X.shape[0]

    check = k

    f = 0.0

    if 'xyb' in M.variables:
        check = check + 2*nb
    if 'zb' in M.variables:
        check = check + nb
    if 'tub' in M.variables:
        tub = variables[check: check + n]
        f += sum(tub**2)
        M.tub = tub
        check = check + n
    if 'tlb' in M.variables:
        tlb = variables[check: check + n]
        f += sum(tlb**2)
        M.tlb = tlb
        check = check + n
    if 'tlb' in M.variables:
        tlb = variables[check: check + n]
        f += sum(tlb**2)
        M.tlb = tlb
        check = check + n
    if 'tub_reac' in M.variables:
        tub_reac = variables[check: check + 2*nb]
        M.tub_reac = tub_reac
        f += sum(tub_reac[:nb]**2 + tub_reac[nb:]**2)  # this is not really 100% ok for 3D. x^2 +y^2 != (x + y)^2
        check = check + 2*nb

    return f


def f_constant(variables, M):
    """Constant or feasible objective function f=1

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    f = 1.0

    return f


def f_reduce_thk(variables, M):
    """Objective function to reduce the thickness of the structure

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.

    Note
    ----
        Please respect the order of the variables.
    """

    return variables[-1]


def f_tight_crosssection(variables, M):
    """Objective function to tight the cross section using normal vectors

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    return -1 * variables[-1]
