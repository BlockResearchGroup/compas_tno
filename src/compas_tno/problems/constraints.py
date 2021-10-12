from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import zeros
from numpy import transpose
from compas_tno.algorithms.equilibrium import zq_from_qid
from scipy.sparse import diags

from compas_tno.problems.bounds_update import ub_lb_update
from compas_tno.problems.bounds_update import b_update

from compas_tno.algorithms import xyz_from_q


__all__ = [
    'constr_wrapper',
    'constr_wrapper_inequalities',
    'constr_wrapper_general',
]


def constr_wrapper(xopt, *args):

    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]

    if 'ind' in variables:  # Not yet deal with all-q
        q[ind] = xopt[:k].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k+len(fixed)].reshape(-1, 1)
    if 't' in variables or 's' in variables or 'n' in variables:
        thk = xopt[-1].item()
        t = shape.datashape['t']

    args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
    z, q = zq_from_qid(q[ind], args)

    constraints = zeros([0, 1])

    if 'funicular' in dict_constr:
        constraints = vstack([constraints, (q[dep] - qmin).reshape(-1, 1)])  # >= 0
    if 'envelope' in dict_constr:
        if 't' in variables or 's' in variables or 'n' in variables:
            ub, lb = ub_lb_update(x, y, thk, t, shape, ub, lb, s, variables)
        constraints = vstack([constraints, ub[ub_ind] - z[ub_ind], z[lb_ind] - lb[lb_ind]])
    if 'reac_bounds' in dict_constr:
        if 't' in variables or 'n' in variables:
            b = b_update(x, y, thk, fixed, shape, b, variables)
        CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
        xyz = hstack([x, y, z])
        p_fixed = hstack([px, py, pz])[fixed]
        R = CfQC.dot(xyz) - p_fixed
        Rx = abs(b[:, 0].reshape(-1, 1)) - multiply(z[fixed] - s[fixed], abs(divide(R[:, 0], R[:, 2]).reshape(-1, 1)))  # >= 0 absoluute numpy
        Ry = abs(b[:, 1].reshape(-1, 1)) - multiply(z[fixed] - s[fixed], abs(divide(R[:, 1], R[:, 2]).reshape(-1, 1)))  # >= 0
        constraints = vstack([constraints, Rx, Ry])
    if 'cracks' in dict_constr:
        crack_tol = 10e-4
        constraints = vstack([constraints, (lb[cracks_lb] - z[cracks_lb]) + crack_tol, (z[cracks_ub] - ub[cracks_ub]) + crack_tol])
    if 'rollers' in dict_constr:
        rx_check = max_rol_rx - abs(Cftx.dot(U.dot(q)) - px[rol_x])
        ry_check = max_rol_ry - abs(Cfty.dot(V.dot(q)) - py[rol_y])
        constraints = vstack([constraints, rx_check, ry_check])
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        A_q = Asym.dot(vstack([q[ind], z[fixed]]))
        constraints = vstack([constraints, A_q])

    return transpose(constraints)[0]  # .reshape(-1, 1)  # Check if this transpose will be a problem for IPOPT


def constr_wrapper_general(variables, M):

    if isinstance(M, list):
        M = M[0]

    k = M.k
    n = M.n
    nb = M.nb
    qid = variables[:k]
    check = k
    M.q = M.B.dot(qid)
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
    if 'lambd' in M.variables:
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

    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)

    constraints = zeros([0, 1])  # missing compression only constraint

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
        if 'adapted-envelope' in M.features:
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
        # constraints in reactions
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

        constraints = vstack([constraints, Rx, Ry])

    return constraints.flatten()


def constr_wrapper_inequalities(xopt, *args):  # This considers a equality in Asym.
    """
    This computes the sensitivities considering only inequality constraints.
    """
    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]
    constraints = constr_wrapper(xopt, *args)
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        Aq = Asym.dot(vstack([q[ind], z[fixed]])).flatten()
        constraints = hstack([constraints, -Aq])

    return constraints
