

from jax.numpy import dot
from jax.numpy import diag
from jax.numpy import zeros
from jax.numpy import array
from jax.numpy import vstack
from jax.numpy import hstack
from jax.numpy.linalg import norm
from jax import jacfwd

from jax.numpy import multiply
from jax.numpy import divide

from compas_tno.algorithms import q_from_variables

from compas_tno.autodiff.jax_equilibrium import xyz_from_q_jax
from compas_tno.autodiff.jax_equilibrium import q_from_variables_jax
from compas_tno.autodiff.jax_swt import weights_from_xyz_jax

from compas_tno.problems.bounds_update import ub_lb_update
from compas_tno.problems.bounds_update import b_update


def f_jacobian_jax():
    """Returns the wrapper of the constraints and Jacobian Matrix function using JAX,

    Returns
    -------
    fconstr : callable
        Callable returning the constraints evaluated in the point.
    fjac : callable
        Callable to compute the jacobian matrix of the constraints.

    """

    fjac = jacfwd(f_constraints_jax)

    return f_constraints_jax, fjac


def f_constraints_jax(variables, M):
    """Wrapper of the constraints assigned using JAX

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    constraints: array (m x 1)
        Vector with the value of the m constraints in the point.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    n = M.n
    nb = M.nb
    check = k
    qid = variables[:k].reshape(-1, 1)
    lambd = 1.0
    thk = M.thk
    t = M.shape.datashape['t']
    P = array(M.P)
    xyb = M.X[M.fixed, :2]
    zb = M.X[M.fixed, 2]

    q = q_from_variables_jax(qid, M.B, M.d)
    Q = diag(q.flatten())

    if 'xyb' in M.variables:
        xyb = variables[check:check + 2*nb]
        check = check + 2*nb
        # X.at[M.fixed, :2].set(xyb.reshape(-1, 2, order='F'))
    if 'zb' in M.variables:
        zb = variables[check: check + nb]
        check = check + nb
        # X.at[M.fixed, [2]].set(zb.flatten())
    if 't' in M.variables or 'n' in M.variables:
        thk = variables[check: check + 1]
        check = check + 1
    if 'lambdh' in M.variables:
        lambd = variables[check: check + 1]
        M.P[:, [0]] = lambd * M.px0
        M.P[:, [1]] = lambd * M.py0
        check = check + 1
    if 'tub' in M.variables:
        tub = variables[check: check + n].reshape(-1, 1)
        M.tub = tub
        check = check + n
    if 'tlb' in M.variables:
        tlb = variables[check: check + n].reshape(-1, 1)
        M.tlb = tlb
        check = check + n
    if 'tub_reac' in M.variables:
        tub_reac = variables[check: check + 2*nb].reshape(-1, 1)
        M.tub_reac = tub_reac
        check = check + 2*nb

    if 'update-loads' in M.features:
        pz_new = -1 * weights_from_xyz_jax(X, M.F, M.V0, M.V1, M.V2, thk=M.thk, density=M.ro)
        P.at[:, 2].set(pz_new)

    Xfixed = hstack([xyb, zb])
    Xfree = xyz_from_q_jax(q, P[array(M.free)], Xfixed, M.Ci, M.Cit, M.Cb)
    X = dot(M.Pmatrix, vstack([Xfree, Xfixed]))

    # X.at[array(M.free)].set(Xfree)
    # print(X)

    # if 'fixed' in M.features:
    #     X[M.free, 2] = X_free[:, 2]
    # else:
    #     X[M.free] = X_free

    constraints = zeros([0, 1])

    if 'funicular' in M.constraints:
        qmin = q.reshape(-1, 1) - M.qmin
        qmax = M.qmax - q.reshape(-1, 1)
        constraints = vstack([constraints, qmin, qmax])

    if 'envelopexy' in M.constraints:
        # constraints in x
        xmin = (X[:, 0] - M.xlimits[:, 0]).reshape(-1, 1)
        xmax = (M.xlimits[:, 1] - X[:, 0]).reshape(-1, 1)
        constraints = vstack([constraints, xmin, xmax])

        # constraints in y
        ymin = (X[:, 1] - M.ylimits[:, 0]).reshape(-1, 1)
        ymax = (M.ylimits[:, 1] - X[:, 1]).reshape(-1, 1)
        constraints = vstack([constraints, ymin, ymax])

    if 'envelope' in M.constraints:
        # constraints in z
        if 'update-envelope' in M.features:
            M.ub, M.lb = ub_lb_update(X[:, 0], X[:, 1], thk, t, M.shape, None, None, M.s, M.variables)
        elif 't' in M.variables or 'n' in M.variables:
            M.ub, M.lb = ub_lb_update(M.x0, M.y0, thk, t, M.shape, M.ub0, M.lb0, M.s, M.variables)
        else:
            pass
        zmin = (X[:, 2] - M.lb.flatten()).reshape(-1, 1)
        zmax = (M.ub.flatten() - X[:, 2]).reshape(-1, 1)
        if 'tub' in M.variables:
            zmax = zmax + tub
        if 'tlb' in M.variables:
            zmin = zmin + tlb
        constraints = vstack([constraints, zmin, zmax])

    if 'reac_bounds' in M.constraints:
        # constraints in compute_reactions
        if 't' in M.variables:
            M.b = b_update(M.x0, M.y0, thk, M.fixed, M.shape, M.b, M.variables)
        elif 'n' in M.variables:
            raise NotImplementedError
        else:
            pass

        CbQC = dot(dot(M.Cb.transpose(), Q), M.C)
        R = dot(CbQC, X) - M.P[M.fixed]
        Rx = abs(M.b[:, [0]].reshape(-1, 1)) - multiply(X[:, [2]][M.fixed] - M.s[M.fixed], abs(divide(R[:, [0]], R[:, [2]]).reshape(-1, 1)))  # >= 0
        Ry = abs(M.b[:, [1]].reshape(-1, 1)) - multiply(X[:, [2]][M.fixed] - M.s[M.fixed], abs(divide(R[:, [1]], R[:, [2]]).reshape(-1, 1)))  # >= 0

        if 'tub_reac' in M.variables:
            Rx = Rx + abs(tub_reac[:nb])
            Ry = Ry + abs(tub_reac[nb:])

        constraints = vstack([constraints, Rx, Ry])

    return constraints.flatten()
