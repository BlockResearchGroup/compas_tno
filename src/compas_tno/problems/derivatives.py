from numpy import zeros
from numpy import identity
from numpy import hstack
from numpy import vstack
from numpy import sign
from numpy import divide
from numpy import sum as npsum
from numpy import multiply

from scipy.sparse.linalg import splu
from scipy.sparse import diags

from compas.numerical import normrow
from numpy import array

from compas_tno.algorithms import xyz_from_q


def d_fobj(fobj, x0, eps, *args):
    """Gradient approximated by hand using finite differences.

    Parameters
    ----------
    fconstr : callable
        Function with the constraints
    x0 : array
        Point to compute the pertubation
    eps : float
        Size of the pertubations

    Returns
    -------
    df0dx : array
        Gradient of the function computed using finite diferences
    """

    f0val = fobj(x0, *args)
    n = len(x0)
    df0dx = zeros((n, 1))
    for i in range(n):
        diff = zeros((n, 1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    return df0dx


def compute_dQ(q, ind, dep, Edinv, Ei):
    """Sensitivity of (all) the force densities with regards to the independent force densities.

    Parameters
    ----------
    q : array
        Force densities
    ind : list
        Indices of the independent edges
    dep : list
        Indices of the dependent edges
    Edinv : array-2d
        Equilibrium matrix of the dependents inverted
    Ei : array-2d
        Equilibrum matrix of the independents

    Returns
    -------
    dQ : array-2d
        Sensitivities of all q's with respect to the independents
    dQdep : array-2d
        Sensitivities of dependent q's with respect to the independents
    """

    dQdep = Edinv.dot(Ei)
    dQ = zeros((len(q), len(ind)))
    dQ[ind] = identity(len(ind))
    dQ[dep] = dQdep[:, :len(ind)]

    return dQ, dQdep


def gradient_feasibility(variables, M):
    """Sensitivity of the feasibility objective function, which returns a null vector.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """
    return zeros((len(variables), 1))


def gradient_reduce_thk(variables, M):
    """Sensitivity of the objective function to minimise the thickness.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """
    grad = zeros((len(variables), 1))
    grad[-1] = 1.0
    return grad


def gradient_tight_crosssection(variables, M):
    """Sensitivity of the objective function to tight the cross section.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """
    grad = zeros((len(variables), 1))
    grad[-1] = -1.0
    return grad


def gradient_fmin(variables, M):
    """Sensitivity of the objective function to minimise the thrust.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    if isinstance(M, list):
        M = M[0]

    n = M.n
    k = M.k
    is_xyb_var = False
    is_zb_var = False

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k]
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
        is_xyb_var = True
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()
        is_zb_var = True

    # update geometry
    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    M.U = diags(M.C.dot(M.X[:, 0]))  # U = diag(Cx)
    M.V = diags(M.C.dot(M.X[:, 1]))  # V = diag(Cy)
    M.W = diags(M.C.dot(M.X[:, 2]))  # W = diag(Cz)

    Q = diags(M.q.ravel())
    CitQCi = M.Cit.dot(Q).dot(M.Ci)
    SPLU_D = splu(CitQCi)

    dxidq = SPLU_D.solve((-M.Cit.dot(M.U)).toarray()).dot(M.B)
    dxdq = zeros((n, k))
    dxdq[M.free] = dxidq

    dyidq = SPLU_D.solve((-M.Cit.dot(M.V)).toarray()).dot(M.B)
    dydq = zeros((n, k))
    dydq[M.free] = dyidq

    CfU = M.Cb.transpose().dot(M.U)
    CfV = M.Cb.transpose().dot(M.V)
    dRxdq = CfU.dot(M.B) + M.Cb.transpose().dot(Q).dot(M.C).dot(dxdq)
    dRydq = CfV.dot(M.B) + M.Cb.transpose().dot(Q).dot(M.C).dot(dydq)

    # print(dRxdq.shape, dRydq.shape)

    Rx = CfU.dot(M.q).reshape(-1, 1) - M.P[M.fixed, 0].reshape(-1, 1)  # check this +/- business
    Ry = CfV.dot(M.q).reshape(-1, 1) - M.P[M.fixed, 1].reshape(-1, 1)
    R = normrow(hstack([Rx, Ry]))
    Rx_over_R = divide(Rx, R)
    Ry_over_R = divide(Ry, R)

    gradient = (Rx_over_R.transpose().dot(dRxdq) + Ry_over_R.transpose().dot(dRydq)).transpose()

    if is_xyb_var:
        dxdxb = zeros((n, nb))
        CitQCf = M.Cit.dot(Q).dot(M.Cb).toarray()
        dxdxb[M.free] = SPLU_D.solve(-CitQCf)
        dxdxb[M.fixed] = identity(nb)
        dydyb = dxdxb
        dRxdx = M.Cb.transpose().dot(Q).dot(M.C).dot(dxdxb)
        dRydy = M.Cb.transpose().dot(Q).dot(M.C).dot(dydyb)
        gradient_xb = (Rx_over_R.transpose().dot(dRxdx)).transpose()
        gradient_yb = (Ry_over_R.transpose().dot(dRydy)).transpose()
        gradient = vstack([gradient, gradient_xb, gradient_yb])
    if is_zb_var:
        gradient = vstack([gradient, zeros((nb, 1))])

    return array(gradient).flatten()


def gradient_fmax(variables, M):
    """Sensitivity of the objective function to maximise the thrust.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    return -1 * gradient_fmin(variables, M)


def gradient_bestfit(variables, M):
    """Sensitivity of the objective function to minimise the vertical squared distance to the target.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    if isinstance(M, list):
        M = M[0]

    n = M.n
    k = M.k

    k = M.k
    nb = len(M.fixed)
    n = M.n

    qid = variables[:k]
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    M.W = diags(M.C.dot(M.X[:, 2]))  # W = diag(Cz)

    f = 2*(M.X[:, [2]] - M.s)

    Q = diags(M.q.ravel())
    CitQCi = M.Cit.dot(Q).dot(M.Ci)
    SPLU_D = splu(CitQCi)

    dzdq = zeros((n, k))
    dzidq = SPLU_D.solve(-M.Cit.dot(M.W).toarray()).dot(M.B)
    dzdq[M.free] = dzidq

    gradient = (f.transpose().dot(dzdq)).transpose()

    if 'zb' in M.variables:
        dz_dzb = zeros((n, nb))
        CitQCf = M.Cit.dot(Q).dot(M.Cb).toarray()
        dz_dzb[M.free] = SPLU_D.solve(-CitQCf)
        dz_dzb[M.fixed] = identity(nb)

        gradient_zb = (f.transpose().dot(dz_dzb)).transpose()
        gradient = vstack([gradient, gradient_zb])

    return gradient


def gradient_horprojection(variables, M):
    """Sensitivity of the objective function to minimise the horizontal squared distance of the nodes on the form diagram to a given pattern.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    if isinstance(M, list):
        M = M[0]

    n = M.n
    k = M.k

    k = M.k
    nb = len(M.fixed)
    n = M.n

    qid = variables[:k]
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    M.U = diags(M.C.dot(M.X[:, 0]))  # U = diag(Cx)
    M.V = diags(M.C.dot(M.X[:, 1]))  # V = diag(Cy)
    # M.W = diags(M.C.dot(M.X[:, 2]))  # W = diag(Cz)

    fx = 2*(M.X[:, [0]] - M.x0)
    fy = 2*(M.X[:, [1]] - M.y0)

    Q = diags(M.q.ravel())
    CitQCi = M.Cit.dot(Q).dot(M.Ci)
    SPLU_D = splu(CitQCi)

    dxdq = zeros((n, k))
    dxidq = SPLU_D.solve(-M.Cit.dot(M.U).toarray()).dot(M.B)
    dxdq[M.free] = dxidq

    dydq = zeros((n, k))
    dyidq = SPLU_D.solve(-M.Cit.dot(M.V).toarray()).dot(M.B)
    dydq[M.free] = dyidq

    # dzdq = zeros((n, k))
    # dzidq = SPLU_D.solve(-M.Cit.dot(M.W).toarray()).dot(M.B)
    # dzdq[M.free] = dzidq

    gradient_x = (fx.transpose().dot(dxdq)).transpose()
    gradient_y = (fy.transpose().dot(dydq)).transpose()

    gradient = gradient_x + gradient_y

    if 'zb' in M.variables:
        gradient_zb = zeros((nb, 1))
        gradient = vstack([gradient, gradient_zb])

    return gradient


def gradient_complementary_energy(variables, M):
    """Sensitivity of the objective function to minimise the complementary energy.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    if isinstance(M, list):
        M = M[0]

    n = M.n
    k = M.k
    is_xyb_var = False
    is_zb_var = False

    k = M.k
    nb = len(M.fixed)

    qid = variables[:k]
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
        is_xyb_var = True
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()
        is_zb_var = True

    # update geometry
    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    M.U = diags(M.C.dot(M.X[:, 0]))  # U = diag(Cx)
    M.V = diags(M.C.dot(M.X[:, 1]))  # V = diag(Cy)
    M.W = diags(M.C.dot(M.X[:, 2]))  # W = diag(Cz)

    Q = diags(M.q.ravel())
    CitQCi = M.Cit.dot(Q).dot(M.Ci)
    SPLU_D = splu(CitQCi)

    dxidq = SPLU_D.solve((-M.Cit.dot(M.U)).toarray()).dot(M.B)
    dxdq = zeros((n, k))
    dxdq[M.free] = dxidq

    dyidq = SPLU_D.solve((-M.Cit.dot(M.V)).toarray()).dot(M.B)
    dydq = zeros((n, k))
    dydq[M.free] = dyidq

    dzidq = SPLU_D.solve(-M.Cit.dot(M.W).toarray()).dot(M.B)
    dzdq = zeros((n, k))
    dzdq[M.free] = dzidq

    CfU = M.Cb.transpose().dot(M.U)
    CfV = M.Cb.transpose().dot(M.V)
    CfW = M.Cb.transpose().dot(M.W)
    dRxdq = CfU.dot(M.B) + M.Cb.transpose().dot(Q).dot(M.C).dot(dxdq)
    dRydq = CfV.dot(M.B) + M.Cb.transpose().dot(Q).dot(M.C).dot(dydq)
    dRzdq = CfW.dot(M.B) + M.Cb.transpose().dot(Q).dot(M.C).dot(dzdq)

    gradient = (M.dXb[:, [0]].transpose().dot(dRxdq) + M.dXb[:, [1]].transpose().dot(dRydq) + M.dXb[:, [2]].transpose().dot(dRzdq)).transpose()

    CitQCf = M.Cit.dot(Q).dot(M.Cb).toarray()
    if is_xyb_var:
        dxdxb = zeros((n, nb))
        dxdxb[M.free] = SPLU_D.solve(-CitQCf)
        dxdxb[M.fixed] = identity(nb)
        dydyb = dxdxb
        dRxdxb = M.Cb.transpose().dot(Q).dot(M.C).dot(dxdxb)
        dRydyb = M.Cb.transpose().dot(Q).dot(M.C).dot(dydyb)
        gradient_xb = (M.dXb[:, [0]].transpose().dot(dRxdxb)).transpose()
        gradient_yb = (M.dXb[:, [1]].transpose().dot(dRydyb)).transpose()
        gradient = vstack([gradient, gradient_xb, gradient_yb])
    if is_zb_var:
        dzdzb = zeros((n, nb))
        dzdzb[M.free] = SPLU_D.solve(-CitQCf)
        dzdzb[M.fixed] = identity(nb)
        dRzdzb = M.Cb.transpose().dot(Q).dot(M.C).dot(dzdzb)
        gradient_zb = (M.dXb[:, [2]].transpose().dot(dRzdzb)).transpose()
        gradient = vstack([gradient, gradient_zb])

    return -1 * array(gradient).flatten()


def gradient_complementary_energy_nonlinear(variables, M):
    """Sensitivity of the objective function to minimise nonlinear complementary energy.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    grad_lin = gradient_complementary_energy(variables, M)

    if M.Ecomp_method == 'simplified':  # assuming area and lengths constant - computed in beginning

        nb = len(M.fixed)

        dEdq_vector = 2 * M.stiff * M.q.reshape(-1, 1)
        dEdq = diags(dEdq_vector.flatten())
        dEdqid = dEdq.dot(M.B)
        grad_quad = npsum(dEdqid, axis=0).reshape(-1, 1)

        if 'xyb' in M.variables:
            grad_quad = vstack([grad_quad, zeros((2*nb, 1))])
        if 'zb' in M.variables:
            grad_quad = vstack([grad_quad, zeros((nb, 1))])

    elif M.Ecomp_method == 'complete':
        grad_quad = gradient_loadpath(variables, M) * M.stiff

    fgrad = grad_lin + grad_quad.flatten()

    return fgrad


def gradient_loadpath(variables, M):
    """Sensitivity of the objective function to minimise the loadpath.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    if isinstance(M, list):
        M = M[0]

    n = M.n
    k = M.k

    k = M.k
    nb = len(M.fixed)
    n = M.n

    qid = variables[:k]
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[k:k + 2*nb]
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[-nb:]
        M.X[M.fixed, [2]] = zb.flatten()

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    uvw = M.C.dot(M.X)
    l2 = npsum(uvw**2, axis=1).reshape(-1, 1)

    M.U = diags(uvw[:, 0])  # U = diag(Cx)
    M.V = diags(uvw[:, 1])  # V = diag(Cy)
    M.W = diags(uvw[:, 2])  # W = diag(Cz)

    Q = diags(M.q.ravel())
    CitQCi = M.Cit.dot(Q).dot(M.Ci)
    SPLU_D = splu(CitQCi)

    dxdq = zeros((n, k))
    dydq = zeros((n, k))
    dzdq = zeros((n, k))

    dxidq = SPLU_D.solve(-M.Cit.dot(M.U).toarray()).dot(M.B)  # not needed if fixed
    dxdq[M.free] = dxidq

    dyidq = SPLU_D.solve(-M.Cit.dot(M.V).toarray()).dot(M.B)  # not needed if fixed
    dydq[M.free] = dyidq

    dzidq = SPLU_D.solve(-M.Cit.dot(M.W).toarray()).dot(M.B)
    dzdq[M.free] = dzidq

    # print(multiply(sign(M.q), l2).transpose().shape)
    # print(multiply(sign(M.q), l2).transpose().dot(M.B).shape)
    # print(M.U.dot(M.C.dot(dxdq)).shape)
    # print(M.V.dot(M.C.dot(dydq)).shape)
    # print(M.W.dot(M.C.dot(dzdq)).shape)
    # print(M.q.transpose().shape)
    # print(abs(M.q.transpose()) * (M.U.dot(M.C.dot(dxdq)) + M.V.dot(M.C.dot(dydq)) + M.W.dot(M.C.dot(dzdq))).shape)
    # print(M.C.dot(dzdq).shape)

    dudq = M.C.dot(dxdq)
    dvdq = M.C.dot(dydq)
    dwdq = M.C.dot(dzdq)

    dldq_1 = (multiply(sign(M.q).reshape(-1, 1), l2).transpose().dot(M.B)).flatten()

    # gradient = (multiply(sign(M.q), l2).transpose().dot(M.B) + 2*abs(M.q.transpose()).dot(M.U.dot(M.C.dot(dxdq)) + M.V.dot(M.C.dot(dydq)) + M.W.dot(M.C.dot(dzdq)))).transpose()
    gradient = zeros((k, 1))

    for i in range(k):
        dlidqi = multiply(uvw[:, [0]], dudq[:, [i]]) + multiply(uvw[:, [1]], dvdq[:, [i]]) + multiply(uvw[:, [2]], dwdq[:, [i]])
        # gradient[i] = sign(M.q[i]) * l2[i] + 2 * abs(M.q.transpose()).dot(dlidqi)
        gradient[i] = dldq_1[i] + 2*abs(M.q.transpose()).dot(dlidqi)

    if 'zb' in M.variables:
        dz_dzb = zeros((n, nb))
        CitQCf = M.Cit.dot(Q).dot(M.Cb).toarray()
        dz_dzb[M.free] = SPLU_D.solve(-CitQCf)
        dz_dzb[M.fixed] = identity(nb)

        gradient_zb = (2*abs(M.q.transpose()).dot(M.W.dot(M.C.dot(dz_dzb)))).transpose().reshape(-1, 1)
        gradient = vstack([gradient, gradient_zb])

    return gradient


def gradient_max_section(variables, M):
    """Sensitivity of the objective function to minimise additional thickness required to find a feasible thrust network.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    grad : array (k x 1)
        The gradient of the objective function in the point. Represents the sensitivity of each variable in the objective function value.
    """

    k = M.k
    nb = len(M.fixed)
    n = M.X.shape[0]

    check = k

    gradient = zeros((k, 1))

    if 'xyb' in M.variables:
        check = check + 2*nb
        gradient = vstack([gradient, zeros((2*nb, 1))])
    if 'zb' in M.variables:
        check = check + nb
        gradient = vstack([gradient, zeros((nb, 1))])
    if 'tub' in M.variables:
        tub = variables[check: check + n].reshape(-1, 1)
        gradient = vstack([gradient, 2 * tub])
        check = check + n
    if 'tlb' in M.variables:
        tlb = variables[check: check + n].reshape(-1, 1)
        gradient = vstack([gradient, 2 * tlb])
        check = check + n
    if 'tub_reac' in M.variables:
        tub_reac = variables[check: check + 2*nb].reshape(-1, 1)
        gradient = vstack([gradient, 2 * tub_reac])  # check if this changes because of the modulus
        check = check + 2*nb

    return gradient
