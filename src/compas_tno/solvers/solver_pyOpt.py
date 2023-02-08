
from .post_process import post_process_general
from compas_tno.algorithms import xyz_from_q
from compas_tno.algorithms import q_from_qid
from numpy import hstack


__all__ = [
    'run_optimisation_pyOpt'
]


def run_optimisation_pyOpt(analysis):
    """Run Optimisation using PyOpt (requires installation).

    Parameters
    ----------
    analysis : :class:`~compas_tno.analysis.Analysis`
        The Analysis object to pass.

    Returns
    -------
    None
        Analysis is executed in place.

    Note
    -------
        Requires pyOpt installation.
    """

    import pyOpt

    form = analysis.form
    optimiser = analysis.optimiser
    solver = optimiser.settings['solver']
    objective = optimiser.settings['objective']
    variables = optimiser.settings['variables']
    fobj = optimiser.fobj
    fconstr = optimiser.fconstr
    args = optimiser.args
    x0 = optimiser.x0
    g0 = optimiser.g0
    f0 = optimiser.fo
    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints,
     cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty) = args
    i_uv = form.index_uv()
    k_i = form.key_index()
    bounds = optimiser.bounds
    x0 = optimiser.x0

    if solver.split('-')[0] == 'pyOpt':
        lower = [lw[0] for lw in bounds]
        upper = [up[1] for up in bounds]
        solver_pyOpt = solver.split('-')[1]
        title = 'Solver: ' + solver_pyOpt + 'Objective: ' + objective

        opt_prob = pyOpt.Optimization(title, pyOpt_wrapper)
        opt_prob.addObj('f', value=f0)
        opt_prob.addVarGroup('x', len(x0), type='c', value=x0, lower=lower, upper=upper)
        opt_prob.addConGroup('c', len(g0), type='i', value=g0)

        if solver_pyOpt == 'SLSQP':
            slv = pyOpt.SLSQP()
            slv.setOption('MAXIT', 1000)

        if solver_pyOpt == 'PSQP':
            slv = pyOpt.PSQP()

        if solver_pyOpt == 'CONMIN':
            slv = pyOpt.CONMIN()

        if solver_pyOpt == 'COBYLA':
            slv = pyOpt.COBYLA()

        if solver_pyOpt == 'SOLVOPT':
            slv = pyOpt.SOLVOPT()

        if solver_pyOpt == 'KSOPT':
            slv = pyOpt.KSOPT()

        if solver_pyOpt == 'NSGA2':
            slv = pyOpt.NSGA2()

        if solver_pyOpt == 'ALGENCAN':
            slv = pyOpt.ALGENCAN()

        if solver_pyOpt == 'FILTERSD':
            slv = pyOpt.FILTERSD()

        if solver_pyOpt == 'SDPEN':
            slv = pyOpt.SDPEN()

        if solver_pyOpt == 'ALPSO':
            slv = pyOpt.SDPEN()

        if solver_pyOpt == 'ALHSO':
            slv = pyOpt.SDPEN()

        if solver_pyOpt == 'MIDACO':
            slv = pyOpt.SDPEN()

        [fopt, xopt, info] = slv(opt_prob, sens_type='FD', args=args, objective=fobj, constraints=fconstr)
        # exitflag = info['text']
        if info['text'] == 'Optimization terminated successfully.':
            exitflag = 0
        else:
            exitflag = 1
        # exitflag = slv.informs
        # if printout:
        #     print(opt_prob.solution(0))
        print(info['text'])
        fopt = fopt.item()
        niter = 100
        if 'zb' in variables:
            q[ind] = xopt[:k].reshape(-1, 1)
            z[fixed] = xopt[k:].reshape(-1, 1)
        else:
            q[ind] = xopt[:k].reshape(-1, 1)

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i)
    q = q_from_qid(q, ind, Edinv, Ei, p)
    Pi = hstack([px, py, pz])
    Xb = hstack([x[fixed], y[fixed], z[fixed]])
    X = xyz_from_q(q, Pi, Xb, Ci, Cit)
    z = X[:, [2]]

    optimiser.niter = niter
    optimiser.exitflag = exitflag
    optimiser.fopt = fopt

    post_process_general(analysis)

    return analysis


def pyOpt_wrapper(x, **kwargs):
    args = kwargs['args']
    fobj = kwargs['objective']
    fconst = kwargs['constraints']
    # Temp
    f = fobj(x, *args)
    g = -1 * fconst(x, *args)  # pyOpt get constraints g1 <= 0 as default
    fail = 0
    return f, g, fail
