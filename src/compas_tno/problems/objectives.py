from numpy import sum as npsum
from numpy.linalg import norm

from compas_tno.algorithms import xyz_from_q
from compas_tno.algorithms import q_from_variables

from scipy.sparse import diags
from scipy.sparse.linalg import splu

from compas_tno.algorithms import weights_from_xyz

from compas_tno.problems import gradient_fmin
from compas_tno.problems import gradient_fmax
from compas_tno.problems import gradient_bestfit
from compas_tno.problems import gradient_horprojection
from compas_tno.problems import gradient_loadpath
from compas_tno.problems import gradient_complementary_energy
from compas_tno.problems import gradient_complementary_energy_nonlinear
from compas_tno.problems import gradient_max_section
from compas_tno.problems import gradient_feasibility
from compas_tno.problems import gradient_reduce_thk
from compas_tno.problems import gradient_tight_crosssection


def objective_selector(objective):
    """Select objective callable and gradient vector based on the desired objective function.
    For the complete list of objectives that can be selected see [link].

    Parameters
    ----------
    objective : str
        The name of the objective. See ``Optimiser`` for a complete list of objectivees.

    Returns
    -------
    fobj : callable
        Callable to compute the value of the objective function in the point.
    fgrad : callable
        Callable to compute the gradient of the objective function in the point.

    """

    if objective == 'loadpath':
        fobj = f_loadpath_general
        fgrad = gradient_loadpath
    elif objective == 'target' or objective == 'bestfit':
        fobj = f_bestfit
        fgrad = gradient_bestfit
    elif objective == 'min':
        fobj = f_min_thrust
        fgrad = gradient_fmin
    elif objective == 'max':
        fobj = f_max_thrust
        fgrad = gradient_fmax
    elif objective == 'feasibility':
        fobj = f_constant
        fgrad = gradient_feasibility
    elif objective == 'hor_projection':
        fobj = f_horprojection
        fgrad = gradient_horprojection
    elif objective == 't':  # analytical reduce thickness
        fobj = f_reduce_thk
        fgrad = gradient_reduce_thk
    elif objective == 's':  # tight UB and LB 0 -> 1/2
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    elif objective == 'n':  # vector n offset the surfaces -> larger the better (higher GSF)
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    elif objective == 'Ecomp-linear':  # vector lambda as hor multiplier larger the better (higher GSF)
        fobj = f_complementary_energy
        fgrad = gradient_complementary_energy
    elif objective == 'Ecomp-nonlinear':
        fobj = f_complementary_energy_nonlinear
        fgrad = gradient_complementary_energy_nonlinear
    elif objective == 'max_section':
        fobj = f_max_section
        fgrad = gradient_max_section
    elif objective == 'max_load':
        fobj = f_tight_crosssection
        fgrad = gradient_tight_crosssection
    else:
        print('Please, provide a valid objective for the optimisation')
        raise NotImplementedError

    return fobj, fgrad


def f_min_thrust(variables, M):
    """Objective function to minimise the horizontal thrust

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : :class:`~compas_tno.problems.Problem`
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    if isinstance(M, list):
        M = M[0]

    update_geometry = False  # In the typical case the geometry does not need to be updated
    k = M.k
    nb = len(M.fixed)

    P_Xh_fixed = M.P[M.fixed][:, :2]  # Horizontal loads in the fixed vertices
    P_free = M.P[M.free]  # Loads in the free vertices

    X = M.X
    Xh = X[:, :2]  # xy of the thrust (i.e. horizontal projection)

    SPLU_D = None

    qid = variables[:k].reshape(-1, 1)
    q = q_from_variables(qid, M.B, M.d)
    Q = diags(q.flatten())

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
        update_geometry = True
    if 'zb' in M.variables:
        zb = variables[-nb:]
        X[M.fixed, 2] = zb.flatten()
        if 'fixed' not in M.features:
            update_geometry = True

    if update_geometry:
        if 'update-loads' in M.features:
            CitQCi = M.Cit @ Q @ M.Ci
            SPLU_D = splu(CitQCi)
            X[M.free] = xyz_from_q(q, P_free, X[M.fixed], M.Ci, M.Cit, M.Cb, SPLU_D=SPLU_D)
            pz = -1 * weights_from_xyz(X, M.F, M.V0, M.V1, M.V2, thk=M.thk, density=M.ro)
            P_free[:, 2] = pz[M.free]

        X[M.free] = xyz_from_q(q, P_free, X[M.fixed], M.Ci, M.Cit, M.Cb, SPLU_D=SPLU_D)
        Xh = X[:, :2]

    CfQC = M.Cb.transpose() @ Q @ M.C
    Rh = CfQC @ Xh - P_Xh_fixed
    f = sum(norm(Rh, axis=1))

    return f


def f_max_thrust(variables, M):
    """Objective function to maximise the horizontal thrust

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
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
    M : :class:`~compas_tno.problems.Problem`
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    f : float
        The value of the objective function.
    """

    return -1 * variables[-1]
