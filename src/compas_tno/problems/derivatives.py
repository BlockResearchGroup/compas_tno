from numpy import zeros
from numpy import identity
from numpy import hstack
from numpy import vstack
from numpy import sign
from numpy import divide
from numpy import sum as npsum
from numpy import multiply

from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import splu
from scipy.sparse import diags

from compas.numerical import normrow
from numpy import array

from compas_tno.algorithms import xyz_from_q

from compas_tno.problems.bounds_update import dub_dlb_update
from compas_tno.problems.bounds_update import db_update


__all__ = [
    'd_fobj',
    'd_fconstr',
    'sensitivities_wrapper',
    'sensitivities_wrapper_inequalities',
    'sensitivities_wrapper_general',
    'gradient_fmin',
    'gradient_fmax',
    'gradient_feasibility',
    'gradient_reduce_thk',
    'gradient_bestfit',
    'gradient_loadpath',
    'gradient_tight_crosssection',
    'gradient_fmin_general',
    'gradient_fmax_general',
    'gradient_bestfit_general',
    'gradient_horprojection_general',
    'gradient_loadpath_general',
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
    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty) = args

    f0val = fobj(x0, *args)
    n = len(x0)
    df0dx = zeros((n, 1))
    for i in range(n):
        diff = zeros((n, 1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    return df0dx


def d_f_ub_lb(fobj, x0, eps, *args):
    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty) = args

    f0val = fobj(x0, *args)
    n = len(x0)
    df0dx = zeros((n, 1))
    for i in range(n):
        diff = zeros((n, 1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    # B = Edinv.dot(Ei)
    # m_qpos = B.shape[0]
    # qpos_contribution = B

    # transpose(vstack([qpos, upper_limit, lower_limit]))[0]

    return df0dx


def sensitivities_wrapper(xopt, *args):

    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]

    if 'ind' in variables:  # Not yet deal with all-q
        q[ind] = xopt[:k].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k+len(fixed)].reshape(-1, 1)
    if 't' in variables or 's' in variables or 'n' in variables:
        thk = xopt[-1].item()
        t = shape.datashape['t']

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
            dub, dlb = dub_dlb_update(x, y, thk, t, shape, ub, lb, s, variables)
            dzub = hstack([-1 * dz, dub])
            dzlb = hstack([dz, - dlb])
            deriv = vstack([deriv, dzub[ub_ind], dzlb[lb_ind]])  # dz IN ub_ind / lb_ind  # ind=indices
        else:
            deriv = vstack([deriv, -dz[ub_ind], dz[lb_ind]])  # dz IN ub_ind / lb_ind
    if 'reac_bounds' in dict_constr:
        CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
        dRxdq = Cf.transpose().dot(U).dot(dQ)
        dRydq = Cf.transpose().dot(V).dot(dQ)
        dRzdq = Cf.transpose().dot(Cz).dot(dQ)
        dRzdzb = Cf.transpose().dot(Q).dot(C).dot(dz_zb)  # new
        dRzdq_comp = Cf.transpose().dot(Q).dot(C).dot(dz[:, :len(ind)])  # new
        xyz = hstack([x, y, z])
        p_fixed = hstack([px, py, pz])[fixed]
        R = CfQC.dot(xyz) - p_fixed
        dslope_dz = zeros((2 * len(fixed), len(fixed)))
        dslope_dind = zeros((2 * len(fixed), len(ind)))
        # print(fixed)
        # print(R)
        for i in range(len(fixed)):
            i_ = len(fixed) + i
            zbi = z[fixed[i]]

            signe_x = 1.0
            signe_y = 1.0
            if R[i, 0] < 0:
                signe_x = -1.0
            if R[i, 1] < 0:
                signe_y = -1.0

            dslope_dz[i, i] = - 1 * abs(R[i, 0]/R[i, 2])
            dslope_dz[i] += - zbi * abs(R[i, 0])/R[i, 2]**2 * dRzdzb[:, i]  # new
            # dslope_dz[i] += + zbi * R[i, 2] * abs(R[i, 0])/abs(R[i, 2])**3 * dRzdzb[:, i]  # new
            # dslope_dind[i] = - signe * z[fixed][i]/R[i, 2]**2 * (abs(R[i, 2]) * dRxdq[i] + signe * abs(R[i, 0]) * (dRzdq[i] + dRzdq_comp[i]))
            dslope_dind[i] = - z[fixed][i]/R[i, 2]**2 * (signe_x * abs(R[i, 2]) * dRxdq[i] + abs(R[i, 0]) * (dRzdq[i] + dRzdq_comp[i]))

            dslope_dz[i_, i] = - 1 * abs(R[i, 1]/R[i, 2])
            dslope_dz[i_] += - zbi * abs(R[i, 1])/R[i, 2]**2 * dRzdzb[:, i]  # new
            # dslope_dz[i_] += + zbi * R[i, 2] * abs(R[i, 1])/abs(R[i, 2])**3 * dRzdzb[:, i]  # new
            # dslope_dind[i_] = - signe * z[fixed][i]/R[i, 2]**2 * (abs(R[i, 2]) * dRydq[i] + signe * abs(R[i, 1]) * (dRzdq[i] + dRzdq_comp[i]))
            dslope_dind[i_] = - z[fixed][i]/R[i, 2]**2 * (signe_y * abs(R[i, 2]) * dRydq[i] + abs(R[i, 1]) * (dRzdq[i] + dRzdq_comp[i]))
        dslope = hstack([dslope_dind, dslope_dz])
        if 't' in variables or 'n' in variables:
            db = db_update(x, y, thk, fixed, shape, b, variables)
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
        if len(variables) == 3:
            Asym = hstack([Asym, zeros((Asym.shape[0], 1))])
        deriv = vstack([deriv, Asym])

    return deriv


def sensitivities_wrapper_inequalities(xopt, *args):
    """
    This computes the sensitivities considering only inequality constraints.
    """
    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]
    deriv = sensitivities_wrapper(xopt, *args)
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        if len(variables) == 3:
            Asym = hstack([Asym, zeros((Asym.shape[0], 1))])
        deriv = vstack([deriv, -Asym])

    return deriv


def compute_dQ(q, ind, dep, Edinv, Ei):

    dQdep = Edinv.dot(Ei)
    dQ = zeros((len(q), len(ind)))
    dQ[ind] = identity(len(ind))
    dQ[dep] = dQdep[:, :len(ind)]

    return dQ, dQdep


def sensitivities_wrapper_general(variables, M):

    if isinstance(M, list):
        M = M[0]

    # variables
    k = M.k  # number of force variables
    n = M.n  # number of vertices
    m = M.m
    nb = len(M.fixed)  # number of fixed vertices
    nbz = 0
    nbxy = 0
    thk = M.thk  # Introduce this because it might be necessary
    t = M.shape.datashape['t']

    qid = variables[:k]
    check = k
    M.q = M.B.dot(qid)

    if 'xyb' in M.variables:
        xyb = variables[check:check + 2*nb]
        check = check + 2*nb
        M.X[M.fixed, :2] = xyb.reshape(-1, 2, order='F')
        nbxy = nb
    if 'zb' in M.variables:
        zb = variables[check: check + nb]
        check = check + nb
        M.X[M.fixed, [2]] = zb.flatten()
        nbz = nb
    if 't' in M.variables or 'n' in M.variables:
        thk = variables[check: check + 1]
        check = check + 1
    elif 'lambd' in M.variables:
        lambd = variables[check: check + 1]
        M.P[:, [0]] = lambd * M.px0
        M.P[:, [1]] = lambd * M.py0
        check = check + 1

    q = M.q
    Xfixed = M.X[M.fixed]

    # update geometry
    M.X[M.free] = xyz_from_q(q, M.P[M.free], Xfixed, M.Ci, M.Cit, M.Cb)

    M.U = diags(M.C.dot(M.X[:, 0]))  # U = diag(Cx)
    M.V = diags(M.C.dot(M.X[:, 1]))  # V = diag(Cy)
    M.W = diags(M.C.dot(M.X[:, 2]))  # W = diag(Cz)

    # initialize jac matrix
    deriv = zeros([0, M.k])

    Q = diags(q.ravel())
    CitQCi = M.Cit.dot(Q).dot(M.Ci)
    SPLU_D = splu(CitQCi)
    nlin_fun = 0
    nlin_limitxy = 0
    nlin_env = 0
    nlin_reacbounds = 0

    A = zeros((n, nb))
    CitQCf = M.Cit.dot(Q).dot(M.Cb).toarray()
    A[M.free] = SPLU_D.solve(-CitQCf)
    A[M.fixed] = identity(nb)

    dxdq = zeros((n, k))
    dydq = zeros((n, k))
    dzdq = zeros((n, k))
    db_column = zeros((nlin_reacbounds, 1))

    if 'funicular' in M.constraints:
        deriv = vstack([deriv, M.B, - M.B])
        nlin_fun = 2 * m

    if 'envelopexy' in M.constraints:
        # jacobian of in constraints on x
        dxidq = SPLU_D.solve(-M.Cit.dot(M.U).toarray()).dot(M.B)
        dxdq[M.free] = dxidq
        deriv = vstack([deriv, dxdq, - dxdq])

        # jacobian of in constraints on y
        dyidq = SPLU_D.solve(-M.Cit.dot(M.V).toarray()).dot(M.B)
        dydq[M.free] = dyidq
        deriv = vstack([deriv, dydq, - dydq])

        nlin_limitxy = 4 * n

    if 'envelope' in M.constraints:
        # jacobian of in constraints on z
        dzidq = SPLU_D.solve(-M.Cit.dot(M.W).toarray()).dot(M.B)
        dzdq[M.free] = dzidq

        if 'adapted-envelope' in M.features:
            dzmaxdt, dzmindt, dzmaxdx, dzmindx, dzmaxdy, dzmindy = dub_dlb_update(M.X[:, 0], M.X[:, 1], thk, t, M.shape, None, None, M.s, M.variables)
            dzmaxdq = dzmaxdx.dot(dxdq) + dzmaxdy.dot(dydq)
            dzmindq = dzmindx.dot(dxdq) + dzmindy.dot(dydq)
            deriv = vstack([deriv, dzdq - dzmindq, dzmaxdq - dzdq])
        else:
            deriv = vstack([deriv, dzdq, - dzdq])

        nlin_env = 2 * n

    if 'reac_bounds' in M.constraints:
        CbQC = M.Cb.transpose().dot(Q).dot(M.C)

        dRxdq = M.Cb.transpose().dot(M.U).dot(M.B) + CbQC.dot(dxdq)
        dRydq = M.Cb.transpose().dot(M.V).dot(M.B) + CbQC.dot(dydq)
        dRzdq = M.Cb.transpose().dot(M.W).dot(M.B) + CbQC.dot(dzdq)

        dRzdzb = CbQC.dot(A)

        R = CbQC.dot(M.X) - M.P[M.fixed]

        dslope_dind = zeros((2 * len(M.fixed), len(M.ind)))
        dslope_dzb = zeros((2 * len(M.fixed), len(M.fixed)))
        dslope_dlambd = zeros((2 * len(M.fixed), 1))

        for i in range(len(M.fixed)):
            i_ = len(M.fixed) + i
            zbi = M.X[M.fixed, 2][i]
            # px0i = M.P[M.fixed, 0][i]
            # py0i = M.P[M.fixed, 1][i]

            signe_x = 1.0
            signe_y = 1.0
            signe_z = 1.0
            if R[i, 0] < 0:
                signe_x = -1.0
            if R[i, 1] < 0:
                signe_y = -1.0
            if R[i, 2] < 0:
                signe_z = -1.0

            dslope_dzb[i, i] = - 1 * abs(R[i, 0]/R[i, 2])
            dslope_dzb[i] += signe_z * zbi * abs(R[i, 0])/R[i, 2]**2 * dRzdzb[:, i]

            # dslope_dlambd[i] = - zbi / abs(R[i, 2]) * signe_x * (- px0i)

            dslope_dind[i] = zbi * signe_x * (-R[i, 2] * dRxdq[i] + R[i, 0] * dRzdq[i]) / R[i, 2]**2 / signe_z

            dslope_dzb[i_, i] = - 1 * abs(R[i, 1]/R[i, 2])
            dslope_dzb[i_] += signe_z * zbi * abs(R[i, 1])/R[i, 2]**2 * dRzdzb[:, i]

            # dslope_dlambd[i_] = - zbi / abs(R[i, 2]) * signe_y * (- py0i)

            dslope_dind[i_] = zbi * signe_y * (-R[i, 2] * dRydq[i] + R[i, 1] * dRzdq[i]) / R[i, 2]**2 / signe_z

        deriv = vstack([deriv, dslope_dind])

        db = db_update(M.x0, M.y0, thk, M.fixed, M.shape, M.b, M.variables)
        db_column = vstack([db[:, 0].reshape(-1, 1), db[:, 1].reshape(-1, 1)])

        nlin_reacbounds = 2 * nb

    if nbxy or nbz:  # add a column to the derivatives to count the variables zb or xyb
        Anull = zeros((n, nb))
        if nbxy:
            deriv = hstack([deriv, vstack([zeros((nlin_fun, nb)), A, -A, Anull, -Anull, Anull, -Anull, zeros((nlin_reacbounds, nb))])])
            deriv = hstack([deriv, vstack([zeros((nlin_fun, nb)), Anull, -Anull, A, -A, Anull, -Anull, zeros((nlin_reacbounds, nb))])])
        if nbz:
            addcolumn = zeros((nlin_fun + nlin_limitxy, nb))
            if 'envelope' in M.constraints:
                addcolumn = vstack([addcolumn, A, -A])
            if 'reac_bounds' in M.constraints:
                addcolumn = vstack([addcolumn, dslope_dzb])
            deriv = hstack([deriv, addcolumn])

    if 't' in M.variables or 'n' in M.variables:  # add a column to the derivatives to count the variable t (thickness)
        if 'adapted-envelope' in M.features:
            pass
        else:
            dzmaxdt, dzmindt = dub_dlb_update(M.x0, M.y0, thk, t, M.shape, M.ub0, M.lb0, M.s, M.variables)[:2]

        dXdt = vstack([zeros((nlin_fun + nlin_limitxy, 1)), -dzmindt, +dzmaxdt, db_column])
        deriv = hstack([deriv, dXdt])

    if 'lambd' in M.variables:  # add a column to the derivatives to count the variable lambd (hor-multiplier)
        if 'fixed' in M.features:
            raise NotImplementedError

        dxdlambd = zeros((n, 1))
        dydlambd = zeros((n, 1))
        dxdlambd[M.free] = SPLU_D.solve(M.px0[M.free]).reshape(-1, 1)
        dydlambd[M.free] = SPLU_D.solve(M.py0[M.free]).reshape(-1, 1)

        if 'adapted-envelope' in M.features:
            dzmaxdlambd = dzmaxdx.dot(dxdlambd) + dzmaxdy.dot(dydlambd)
            dzmindlambd = dzmindx.dot(dxdlambd) + dzmindy.dot(dydlambd)
            dXdlambd = vstack([zeros((nlin_fun, 1)), dxdlambd, - dxdlambd, dydlambd, - dydlambd, - dzmindlambd, +dzmaxdlambd])
        else:
            dXdlambd = vstack([zeros((nlin_fun, 1)), dxdlambd, - dxdlambd, dydlambd, - dydlambd, zeros((nlin_env, 1))])

        if 'reac_bounds' in M.constraints:
            dXdlambd = vstack([dXdlambd, dslope_dlambd])

        deriv = hstack([deriv, dXdlambd])

    return deriv


def gradient_fmin(xopt, *args):

    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]

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
    elif objective == 'target' or objective == 'bestfit':
        grad = gradient_feasibility
    else:
        raise NotImplementedError

    return grad


def gradient_bestfit(xopt, *args):

    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]

    q[ind] = xopt[:k].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k + len(fixed)].reshape(-1, 1)

    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

    f = 2*(z - s)
    grad = zeros((len(xopt), 1))

    dz = zeros([len(z), len(xopt)])
    Q = diags(q.ravel())
    CitQCi = Cit.dot(Q).dot(Ci)
    SPLU_D = splu(CitQCi)
    dQ, dQdep = compute_dQ(q, ind, dep, Edinv, Ei)
    Cz = diags((C.dot(z)).ravel())
    B = - Cit.dot(Cz)
    dz[free, :k] = SPLU_D.solve(B.dot(dQ))

    if 'zb' in variables:
        B = - Cit.dot(Q).dot(Cf).toarray()
        dz[free, k:] = SPLU_D.solve(B)
        dz[fixed, k:] = identity(len(fixed))

    for j in range(len(xopt)):
        grad[j] = f.transpose().dot(dz[:, j].reshape(-1, 1))

    return grad


def gradient_loadpath(xopt, *args):

    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub,
     free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape) = args[:50]

    q[ind] = xopt[:k].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k + len(fixed)].reshape(-1, 1)

    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

    grad = zeros((len(xopt), 1))

    dZ = zeros((len(z), k))

    Q = diags(q.ravel())
    CitQCi = Cit.dot(Q).dot(Ci)
    SPLU_D = splu(CitQCi)
    dQ, dQdep = compute_dQ(q, ind, dep, Edinv, Ei)
    Cz = diags((C.dot(z)).ravel())
    B = - Cit.dot(Cz)
    dZi = SPLU_D.solve(B.dot(dQ))
    dZ[free, :] = dZi
    dW = C.dot(dZ)

    WdW = 2*Cz.dot(dW)
    l2 = lh + C.dot(z)**2

    for j in range(k):
        grad[j] = dQ[:, j].reshape(1, -1).dot(l2) + q.transpose().dot(WdW[:, j].reshape(-1, 1))

    if 'zb' in variables:  # This is correct but can lead to unexpected results for symmetric problems, try to constraint to symmetry
        B = - Cit.dot(Q).dot(Cf).toarray()
        dZ = zeros((len(z), len(fixed)))
        dZ[free] = SPLU_D.solve(B)
        dZ[fixed] = identity(len(fixed))
        dW = C.dot(dZ)
        WdW = 2*Cz.dot(dW)

        for i in range(len(fixed)):
            grad[j+i+1] = q.transpose().dot(WdW[:, i].reshape(-1, 1))

    return grad


def gradient_fmin_general(variables, M):

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


def gradient_fmax_general(variables, M):

    return -1 * gradient_fmin_general(variables, M)


def gradient_bestfit_general(variables, M):

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


def gradient_horprojection_general(variables, M):

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


def gradient_loadpath_general(variables, M):

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
    # print(sign(M.q).shape)
    # print(l2.shape)
    # print(multiply(sign(M.q), l2).transpose().shape)
    # print(M.B.shape)
    # print(dldq_1.shape)

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

        gradient_zb = 2*abs(M.q.transpose()) * (M.W.dot(M.C.dot(dz_dzb)))
        gradient = vstack([gradient, gradient_zb])

    return gradient
