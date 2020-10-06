from numpy import zeros
from numpy import identity
from numpy import hstack
from numpy import vstack

from numpy import divide

from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import splu
from scipy.sparse import diags

from compas.numerical import normrow


from compas_tno.problems.bounds_update import dub_dlb_update
from compas_tno.problems.bounds_update import db_update


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'd_fobj',
    'd_fconstr',
    'sensitivities_wrapper',
    'sensitivities_wrapper_inequalities',
    'gradient_fmin',
    'gradient_fmax',
    'gradient_feasibility',
    'gradient_reduce_thk',
    'gradient_tight_crosssection'
]


# Gradient "approximated by hand"
def d_fobj(fobj, x0, eps, *args):
    f0val = fobj(x0, *args)
    n = len(x0)
    df0dx = zeros((n, 1))
    for i in range(n):
        diff = zeros((n, 1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    return df0dx


# Jacobian "approximated by hand"
def d_fconstr(fconstr, x0, eps, *args):
    fval = fconstr(x0, *args).reshape(-1, 1)
    m = len(fval)
    n = len(x0)
    dfdx = zeros((m, n))
    for i in range(m):
        for j in range(n):
            diff = zeros((n, 1))
            diff[j] = eps
            dfdx[i, j] = (fconstr(x0 + diff, *args).reshape(-1, 1) - fval)[i]/diff[j]
    return dfdx


def d_min_thrust(fobj, x0, eps, *args):
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    f0val = fobj(x0, *args)
    n = len(x0)
    df0dx = zeros((n, 1))
    for i in range(n):
        diff = zeros((n, 1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    return df0dx


def d_f_ub_lb(fobj, x0, eps, *args):
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    f0val = fobj(x0, *args)
    n = len(x0)
    df0dx = zeros((n, 1))
    for i in range(n):
        diff = zeros((n, 1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    B = Edinv.dot(Ei)
    m_qpos = B.shape[0]
    qpos_contribution = B

    # transpose(vstack([qpos, upper_limit, lower_limit]))[0]

    return df0dx


# def sensitivities_wrapper(xopt, *args):

#     q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym = args[:48]
#     zf_include = False
#     if len(xopt) > k:
#         q[ind], z[fixed] = xopt[:k].reshape(-1, 1), xopt[k:].reshape(-1, 1)
#         zf_include = True
#     else:
#         q[ind] = xopt

#     q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
#     z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

#     deriv = zeros([0, len(xopt)])

#     if 'funicular' in dict_constr:
#         dQdep = Edinv.dot(Ei)
#         if zf_include:
#             dQdep = hstack([dQdep, zeros((len(dep), len(fixed)))])
#         deriv = vstack([deriv, dQdep])
#     if 'envelope' in dict_constr:
#         Q = diags(q.ravel())
#         CitQCi = Cit.dot(Q).dot(Ci)
#         SPLU_D = splu(CitQCi)
#         dQ = zeros((len(q), len(ind)))
#         dQ[ind] = identity(len(ind))
#         dQ[dep] = dQdep[:, :len(ind)]
#         Cz = diags((C.dot(z)).ravel())
#         B = - Cit.dot(Cz)
#         dz = SPLU_D.solve(B.dot(dQ))
#         if zf_include:
#             B = - Cit.dot(Q).dot(Cf).toarray()
#             dzi = SPLU_D.solve(B)
#             dz = hstack([dz, dzi])
#         deriv = vstack([deriv, -dz, dz])  # dz internal only -> Change!
#     if 'reac_bounds' in dict_constr:
#         CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
#         dRxdq = Cf.transpose().dot(U).dot(dQ)
#         dRydq = Cf.transpose().dot(V).dot(dQ)
#         dRzdq = Cf.transpose().dot(Cz).dot(dQ)
#         xyz = hstack([x, y, z])
#         p_fixed = hstack([px, py, pz])[fixed]
#         R = CfQC.dot(xyz) - p_fixed
#         dslope_dz = zeros((2 * len(fixed), len(fixed)))
#         dslope_dind = zeros((2 * len(fixed), len(ind)))
#         for i in range(len(fixed)):
#             i_ = len(fixed) + i
#             if R[i, 0]*R[i, 2] >= 0:
#                 signe = - 1.0
#             else:
#                 signe = + 1.0
#             dslope_dz[i, i] = signe * R[i, 0]/R[i, 2]
#             dslope_dind[i] = - signe * z[fixed][i]/R[i, 2]**2 * (abs(R[i, 2]) * dRxdq[i] + signe * abs(R[i, 0]) * dRzdq[i])
#             if R[i, 1]*R[i, 2] >= 0:
#                 signe = - 1.0
#             else:
#                 signe = + 1.0
#             dslope_dz[i_, i] = signe * R[i, 1]/R[i, 2]
#             dslope_dind[i_] = - signe * z[fixed][i]/R[i, 2]**2 * (abs(R[i, 2]) * dRydq[i] + signe * abs(R[i, 1]) * dRzdq[i])
#         dslope = hstack([dslope_dind, dslope_dz])
#         deriv = vstack([deriv, dslope])
#     if 'cracks' in dict_constr:
#         deriv = vstack([deriv, -dz[cracks_lb], dz[cracks_ub]])
#     if 'rollers' in dict_constr:
#         drolx = Cftx.dot(U).dot(dQ)
#         droly = Cfty.dot(V).dot(dQ)
#         rx = Cftx.dot(U.dot(q)) - px[rol_x]
#         ry = Cfty.dot(V.dot(q)) - py[rol_y]
#         drol = zeros((len(rol_x)+len(rol_y), len(xopt)))
#         for i in range(len(rol_x)):
#             if rx[i] < 0:
#                 drol[i, :len(ind)] = + drolx[i]
#             else:
#                 drol[i, :len(ind)] = - drolx[i]
#         for i in range(len(rol_y)):
#             if ry[i] < 0:
#                 drol[i + len(rol_x), :len(ind)] = + droly[i]
#             else:
#                 drol[i + len(rol_x), :len(ind)] = - droly[i]
#         deriv = vstack([deriv, drol])
#     if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
#         deriv = vstack([deriv, Asym])

#     return deriv


# def sensitivities_wrapper_inequalities(xopt, *args):
#     """
#     This computes the sensitivities considering only inequality constraints.
#     """
#     q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym = args[:48]
#     deriv = sensitivities_wrapper(xopt, *args)
#     if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
#         deriv = vstack([deriv, -Asym])

#     return deriv


def sensitivities_wrapper(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape = args[:50]

    if 'ind' in variables:  # Not yet deal with all-q
        q[ind] = xopt[:k].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k+len(fixed)].reshape(-1, 1)
    if 't' in variables or 's' in variables or 'n' in variables:
        thk = xopt[-1].item()
        t = shape.data['t']

    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

    deriv = zeros([0, len(xopt)])

    if 'funicular' in dict_constr:
        dQdep = Edinv.dot(Ei)
        if 'zb' in variables:
            dQdep = hstack([dQdep, zeros((len(dep), len(fixed)))])
        if 't' in variables or 's' in variables or 'n' in variables:
            dQdep = hstack([dQdep, zeros((len(dep), 1))])
        deriv = vstack([deriv, dQdep])
    if 'envelope' in dict_constr:
        Q = diags(q.ravel())
        CitQCi = Cit.dot(Q).dot(Ci)
        SPLU_D = splu(CitQCi)
        dQ = zeros((len(q), len(ind)))
        dQ[ind] = identity(len(ind))
        dQ[dep] = dQdep[:, :len(ind)]
        Cz = diags((C.dot(z)).ravel())
        B = - Cit.dot(Cz)
        dz = zeros((len(z), len(ind)))
        dz[free] = SPLU_D.solve(B.dot(dQ))
        if 'zb' in variables:
            B = - Cit.dot(Q).dot(Cf).toarray()
            dz_zb = zeros((len(z), len(fixed)))
            dz_zb[free] = SPLU_D.solve(B)
            dz_zb[fixed] = identity(len(fixed))
            dz = hstack([dz, dz_zb])
        if 't' in variables or 's' in variables or 'n' in variables:
            dub, dlb = dub_dlb_update(x, y, thk, t, shape, ub, lb, variables)
            dzub = hstack([-dz, dub])
            dzlb = hstack([dz, - dlb])
            deriv = vstack([deriv, dzub, dzlb])  # dz IN ALL Z's
        else:
            deriv = vstack([deriv, -dz, dz])  # dz IN ALL Z's
    if 'reac_bounds' in dict_constr:
        CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
        dRxdq = Cf.transpose().dot(U).dot(dQ)
        dRydq = Cf.transpose().dot(V).dot(dQ)
        dRzdq = Cf.transpose().dot(Cz).dot(dQ)
        xyz = hstack([x, y, z])
        p_fixed = hstack([px, py, pz])[fixed]
        R = CfQC.dot(xyz) - p_fixed
        dslope_dz = zeros((2 * len(fixed), len(fixed)))
        dslope_dind = zeros((2 * len(fixed), len(ind)))
        for i in range(len(fixed)):
            i_ = len(fixed) + i
            if R[i, 0]*R[i, 2] >= 0:
                signe = - 1.0
            else:
                signe = + 1.0
            dslope_dz[i, i] = signe * R[i, 0]/R[i, 2]
            dslope_dind[i] = - signe * z[fixed][i]/R[i, 2]**2 * (abs(R[i, 2]) * dRxdq[i] + signe * abs(R[i, 0]) * dRzdq[i])
            if R[i, 1]*R[i, 2] >= 0:
                signe = - 1.0
            else:
                signe = + 1.0
            dslope_dz[i_, i] = signe * R[i, 1]/R[i, 2]
            dslope_dind[i_] = - signe * z[fixed][i]/R[i, 2]**2 * (abs(R[i, 2]) * dRydq[i] + signe * abs(R[i, 1]) * dRzdq[i])
        dslope = hstack([dslope_dind, dslope_dz])
        if 't' in variables:
            db = db_update(x, y, thk, fixed, shape)
            db_column = vstack([db[:, 0].reshape(-1, 1), db[:, 1].reshape(-1, 1)])
            dslope = hstack([dslope, db_column])
        deriv = vstack([deriv, dslope])
    if 'cracks' in dict_constr:
        deriv = vstack([deriv, -dz[cracks_lb], dz[cracks_ub]])
    if 'rollers' in dict_constr:
        drolx = Cftx.dot(U).dot(dQ)
        droly = Cfty.dot(V).dot(dQ)
        rx = Cftx.dot(U.dot(q)) - px[rol_x]
        ry = Cfty.dot(V.dot(q)) - py[rol_y]
        drol = zeros((len(rol_x)+len(rol_y), len(xopt)))
        for i in range(len(rol_x)):
            if rx[i] < 0:
                drol[i, :len(ind)] = + drolx[i]
            else:
                drol[i, :len(ind)] = - drolx[i]
        for i in range(len(rol_y)):
            if ry[i] < 0:
                drol[i + len(rol_x), :len(ind)] = + droly[i]
            else:
                drol[i + len(rol_x), :len(ind)] = - droly[i]
        deriv = vstack([deriv, drol])
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        deriv = vstack([deriv, Asym])

    return deriv


def sensitivities_wrapper_inequalities(xopt, *args):
    """
    This computes the sensitivities considering only inequality constraints.
    """
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape = args[:50]
    deriv = sensitivities_wrapper(xopt, *args)
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        deriv = vstack([deriv, -Asym])

    return deriv


def compute_dQ(q, ind, dep, Edinv, Ei):

    dQdep = Edinv.dot(Ei)
    dQ = zeros((len(q), len(ind)))
    dQ[ind] = identity(len(ind))
    dQ[dep] = dQdep[:, :len(ind)]

    return dQ, dQdep


def gradient_fmin(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape_data = args[:50]

    if 'ind' in variables:
        q[ind] = xopt[:k].reshape(-1, 1)
    else:
        q = xopt[:len(q)].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k+len(fixed)].reshape(-1, 1)

    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

    dQ, _ = compute_dQ(q, ind, dep, Edinv, Ei)
    CfU = Cf.transpose().dot(U)
    CfV = Cf.transpose().dot(V)
    CfUdQ = CfU.dot(dQ)
    CfVdQ = CfV.dot(dQ)
    Rx = CfU.dot(q) - px[fixed]
    Ry = CfV.dot(q) - py[fixed]
    R = normrow(hstack([Rx, Ry]))
    Rx_over_R = divide(Rx, R)
    Ry_over_R = divide(Ry, R)

    gradient = (Rx_over_R.transpose().dot(CfUdQ) + Ry_over_R.transpose().dot(CfVdQ)).transpose()
    if 'zb' in variables:
        gradient = vstack([gradient, zeros((len(fixed), 1))])

    return gradient


def gradient_fmax(xopt, *args):
    return -1 * gradient_fmin(xopt, *args)


def gradient_feasibility(xopt, *args):
    return zeros((len(xopt), 1))


def gradient_reduce_thk(xopt, *args):
    grad = zeros((len(xopt), 1))
    grad[-1] = 1.0
    return grad


def gradient_tight_crosssection(xopt, *args):
    grad = zeros((len(xopt), 1))
    grad[-1] = -1.0
    return grad


def gradient_wrapper(xopt, *args):
    # WIP

    objective = 'min'

    if objective == 'min':
        grad = gradient_fmin
    elif objective == 'max':
        grad = gradient_fmax
    elif objective == 'feasibility':
        grad = gradient_feasibility
    else:
        raise NotImplementedError

    return grad
