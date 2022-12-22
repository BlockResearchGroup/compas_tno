from numpy import vstack
from numpy import multiply
from numpy import divide
from numpy import zeros
from scipy.sparse import diags

from compas_tno.problems.bounds_update import ub_lb_update
from compas_tno.problems.bounds_update import b_update

from compas_tno.algorithms import q_from_variables
from compas_tno.algorithms import xyz_from_q
from compas_tno.algorithms import weights_from_xyz


def constr_wrapper(variables, M):
    """Wrapper of the constraints assigned.

    Parameters
    ----------
    variables : array [(]k x 1]
        Variables to pass to the function.
    M : :class:`~compas_tno.problems.Problem`
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    constraints: array [m x 1]
        Vector with the value of the m constraints in the point.
    """

    if isinstance(M, list):
        M = M[0]

    k = M.k
    n = M.n
    nb = M.nb
    check = k
    qid = variables[:k].reshape(-1, 1)
    lambdh = 1.0
    delta = 0.0
    thk = M.thk
    t = M.shape.datashape['t']

    if 'xyb' in M.variables:
        xyb = variables[check:check + 2*nb]
        check = check + 2*nb
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
    if 'zb' in M.variables:
        zb = variables[check: check + nb]
        check = check + nb
        M.X[M.fixed, [2]] = zb.flatten()
    if 't' in M.variables or 'n' in M.variables:
        thk = variables[check: check + 1]
        check = check + 1
    if 'lambdh' in M.variables:
        lambdh = variables[check: check + 1]
        M.P[:, [0]] = lambdh * M.px0
        M.P[:, [1]] = lambdh * M.py0
        M.d = lambdh * M.d0
        # M.d[M.dep] = -M.Edinv.dot(vstack([M.P[M.free_x, 0].reshape(-1, 1), M.P[M.free_y, 1].reshape(-1, 1)]))  # make this line shorter
        check = check + 1
    if 'lambdv' in M.variables:
        lambdv = variables[check: check + 1]
        M.P[:, [2]] = lambdv * M.pzv + M.pz0
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
    if 'delta' in M.variables:
        delta = float(variables[check: check + 1])
        check = check + 1

    M.q = q_from_variables(qid, M.B, M.d)

    if 'update-loads' in M.features:
        M.P[:, 2] = -1 * weights_from_xyz(M.X, M.F, M.V0, M.V1, M.V2, thk=M.thk, density=M.ro)

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    # if 'fixed' in M.features:
    #     M.X[M.free, 2] = X_free[:, 2]
    # else:
    #     M.X[M.free] = X_free

    constraints = zeros([0, 1])

    if 'funicular' in M.constraints:
        qmin = M.q.reshape(-1, 1) - M.qmin
        qmax = M.qmax - M.q.reshape(-1, 1)
        constraints = vstack([constraints, qmin, qmax])

    if 'envelopexy' in M.constraints:
        # constraints in x
        xmin = (M.X[:, 0] - M.xlimits[:, 0]).reshape(-1, 1)
        xmax = (M.xlimits[:, 1] - M.X[:, 0]).reshape(-1, 1)
        constraints = vstack([constraints, xmin, xmax])

        # constraints in y
        ymin = (M.X[:, 1] - M.ylimits[:, 0]).reshape(-1, 1)
        ymax = (M.ylimits[:, 1] - M.X[:, 1]).reshape(-1, 1)
        constraints = vstack([constraints, ymin, ymax])

    if 'envelope' in M.constraints:
        # constraints in z
        if 'update-envelope' in M.features:
            M.ub, M.lb = ub_lb_update(M.X[:, 0], M.X[:, 1], thk, t, M.shape, None, None, M.s, M.variables)
        elif 't' in M.variables or 'n' in M.variables:
            M.ub, M.lb = ub_lb_update(M.x0, M.y0, thk, t, M.shape, M.ub0, M.lb0, M.s, M.variables)
        else:
            pass
        zmin = (M.X[:, 2] - M.lb.flatten()).reshape(-1, 1)
        zmax = (M.ub.flatten() - M.X[:, 2]).reshape(-1, 1)
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

        CbQC = M.Cb.transpose().dot(diags(M.q.flatten())).dot(M.C)
        R = CbQC.dot(M.X) - M.P[M.fixed]
        Rx = abs(M.b[:, [0]].reshape(-1, 1)) - multiply(M.X[:, [2]][M.fixed] - M.s[M.fixed], abs(divide(R[:, [0]], R[:, [2]]).reshape(-1, 1)))  # >= 0
        Ry = abs(M.b[:, [1]].reshape(-1, 1)) - multiply(M.X[:, [2]][M.fixed] - M.s[M.fixed], abs(divide(R[:, [1]], R[:, [2]]).reshape(-1, 1)))  # >= 0

        if 'tub_reac' in M.variables:
            Rx = Rx + abs(tub_reac[:nb])
            Ry = Ry + abs(tub_reac[nb:])

        constraints = vstack([constraints, Rx, Ry])

    if 'displ_map' in M.constraints:
        if delta:
            pos = (M.E @ M.q + delta * M.Ed @ M.q - M.ph).reshape(-1, 1)
        else:
            pos = (M.E @ M.q - M.ph).reshape(-1, 1)
        # pos = (M.E @ M.q - M.ph).reshape(-1, 1)
        neg = -pos

        constraints = vstack([constraints, pos, neg])

    return constraints.flatten()
