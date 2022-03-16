from numpy import zeros
from numpy import identity
from numpy import hstack
from numpy import vstack

from scipy.sparse.linalg import splu
from scipy.sparse import diags

from compas_tno.problems.bounds_update import dub_dlb_update
from compas_tno.problems.bounds_update import db_update

from compas_tno.algorithms import q_from_variables
from compas_tno.algorithms import xyz_from_q


def d_fconstr(fconstr, x0, eps, *args):
    """Jacobian matrix approximated using finite differences.

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
    dfdx : array
        Jacobian of the constraints computed using finite diferences
    """

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


def sensitivities_wrapper(variables, M):
    """Jacobian matrix computed analytically based on the constraints and variables assigned.

    Parameters
    ----------
    variables : array (k x 1)
        Variables to pass to the function.
    M : Problem
        The class with necessary matrices, or arguments, to compute the objective function

    Returns
    -------
    [type]
        [description]

    Notes
    -----
        Observe the order in which variables are added.
    """

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
    lambd = 1.0

    qid = variables[:k].reshape(-1, 1)
    check = k

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
    if 'tub' in M.variables:
        tub = variables[check: check + n]
        M.tub = tub
        check = check + n
    if 'tlb' in M.variables:
        tlb = variables[check: check + n]
        M.tlb = tlb
        check = check + n
    if 'tub_reac' in M.variables:
        tub_reac = variables[check: check + 2*nb].reshape(-1, 1)
        M.tub_reac = tub_reac
        check = check + 2*nb
    if 'lambdv' in M.variables:
        lambdv = variables[check: check + 1]
        M.P[:, [2]] = lambdv * M.pzv + M.pz0
        check = check + 1

    M.q = q_from_variables(qid, M.B, M.d, lambd=lambd)

    # update geometry
    M.X[M.free] = xyz_from_q(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)
    # if 'fixed' in M.features:
    #     M.X[M.free, 2] = X_free[:, 2]
    # else:
    #     M.X[M.free] = X_free

    q = M.q

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

    if 'lambdv' in M.variables:  # add a column to the derivatives to count the variable lambdv (vertical load multiplier)

        dzdlambdv = zeros((n, 1))
        dzdlambdv[M.free] = SPLU_D.solve(M.pzv[M.free]).reshape(-1, 1)

        if 'fixed' in M.features:
            dXdlambd = vstack([zeros((nlin_fun, 1)), dzdlambdv, - dzdlambdv])
        else:
            raise NotImplementedError
            if 'update-envelope' in M.features:
                dzmaxdlambd = dzmaxdx.dot(dxdlambd) + dzmaxdy.dot(dydlambd)
                dzmindlambd = dzmindx.dot(dxdlambd) + dzmindy.dot(dydlambd)
                dXdlambd = vstack([zeros((nlin_fun, 1)), dxdlambd, - dxdlambd, dydlambd, - dydlambd, - dzmindlambd, +dzmaxdlambd])
            else:
                dXdlambd = vstack([zeros((nlin_fun, 1)), dxdlambd, - dxdlambd, dydlambd, - dydlambd, zeros((nlin_env, 1))])

        if 'reac_bounds' in M.constraints:
            dRzdlambdv = CbQC.dot(dzdlambdv) - M.pzv[M.fixed]
            dslope_dlambdv = zeros((2 * nb, 1))
            for i in range(nb):
                i_ = nb + i
                zbi = M.X[M.fixed, 2][i]
                signe_z = 1.0
                if R[i, 2] < 0:
                    signe_z = -1.0

                dslope_dlambdv[i] = signe_z * zbi * abs(R[i, 0])/R[i, 2]**2 * dRzdlambdv[i]
                dslope_dlambdv[i_] = signe_z * zbi * abs(R[i, 1])/R[i, 2]**2 * dRzdlambdv[i]
            dXdlambd = vstack([dXdlambd, dslope_dlambdv])


        deriv = hstack([deriv, dXdlambd])

    if 'tub' in M.variables:  # add a column to the derivatives to count the variable tub (max_section)

        nconst = deriv.shape[0]
        dXdtub = zeros((nconst, n))
        startline = nlin_fun + nlin_limitxy
        endline = startline + nlin_env

        In = identity(n)
        I0 = 0.0 * In
        Mt = vstack([I0, In])

        dXdtub[startline: endline, :] = Mt

        deriv = hstack([deriv, dXdtub])

    if 'tlb' in M.variables:  # add a column to the derivatives to count the variable tub (max_section)

        nconst = deriv.shape[0]
        dXdtub = zeros((nconst, n))
        startline = nlin_fun + nlin_limitxy
        endline = startline + nlin_env

        In = identity(n)
        I0 = 0.0 * In
        Mt = vstack([In, I0])

        dXdtub[startline: endline, :] = Mt

        deriv = hstack([deriv, dXdtub])

    if 'tub_reac' in M.variables:  # add a column to the derivatives to count the variable tub (max_section)

        nconst = deriv.shape[0]
        dXdtreacub = zeros((nconst, 2*nb))
        startline = nlin_fun + nlin_limitxy + nlin_env
        endline = startline + nlin_reacbounds

        Inb = identity(nb)  # check if the modulus need to be added here
        Inb0 = 0.0 * Inb
        Mt = hstack([vstack([Inb, Inb0]), vstack([Inb0, Inb])])

        dXdtreacub[startline: endline, :] = Mt

        deriv = hstack([deriv, dXdtreacub])

    return deriv
