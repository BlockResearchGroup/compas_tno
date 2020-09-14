from scipy.optimize import fmin_slsqp
from scipy.optimize import shgo

from compas.numerical import devo_numpy
from compas.numerical import ga

from compas_tno.algorithms import reactions
from compas_tno.algorithms import zlq_from_qid

from compas.utilities import geometric_key

import time


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'run_optimisation_scipy'
]


def run_optimisation_scipy(analysis):
    """ Run nonlinear optimisation problem with scipy

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    obj : analysis
        Analysis object optimised.

    """

    form = analysis.form
    optimiser = analysis.optimiser
    solver = optimiser.data['solver']
    variables = optimiser.data['variables']
    fobj = optimiser.fobj
    fconstr = optimiser.fconstr
    fgrad = optimiser.fgrad
    fjac = optimiser.fjac
    args = optimiser.args
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints, max_rol_rx, max_rol_ry, Asym = args
    i_uv = form.index_uv()
    i_k = form.index_key()
    k_i = form.key_index()
    bounds = optimiser.bounds
    x0 = optimiser.x0
    plot = optimiser.data.get('plot', False)
    printout = optimiser.data.get('printout', True)
    summary = optimiser.data.get('summary', False)
    grad_choice = optimiser.data.get('gradient', False)
    jac_choice = optimiser.data.get('jacobian', False)

    if grad_choice is False:
        fgrad = None
    if jac_choice is False:
        fjac = None

    start_time = time.time()

    if solver == 'slsqp' or solver == 'SLSQP':
        fopt, xopt, exitflag, niter, message = _slsqp(fobj, x0, bounds, fgrad, fjac, printout, fconstr, args)
        while exitflag == 9:
            fopt, xopt, exitflag, niter, message = _slsqp(fobj, xopt, bounds, fgrad, fjac, printout, fconstr, args)
        if exitflag == 0:
            if 'zb' in variables:
                q[ind] = xopt[:k].reshape(-1, 1)
                z[fixed] = xopt[k:].reshape(-1, 1)
            else:
                q[ind] = xopt[:k].reshape(-1, 1)  # Code the option with only qinds WIP
        else:
            if printout:
                print(message)
    elif solver == 'shgo':
        dict_constr = []
        for i in range(len(fconstr(x0, *args))):
            args_constr = list(args)
            args_constr.append(i)
            args_constr.append(fconstr)
            dict_ = {
                'type': 'ineq',
                'fun': _shgo_constraint_wrapper,
                'args': args_constr,
            }
            dict_constr.append(dict_)
        args_constr[len(args_constr)-2] = 0
        result = _shgo(fobj, bounds, True, dict_constr, args)
        fopt = result['fun']
        xopt = result['x']
        sucess = result['success']
        message = result['message']
        if sucess is True:
            exitflag = 0
            q[ind] = xopt[:k].reshape(-1, 1)
            z[fixed] = xopt[k:].reshape(-1, 1)
        else:
            print(message)

    elapsed_time = time.time() - start_time
    if printout:
        print('Solving Time: {0:.1f} sec'.format(elapsed_time))
        optimiser.time = elapsed_time

    g_final = fconstr(xopt, *args)
    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i)
    z, _, q, q_ = zlq_from_qid(q[ind], args)

    gkeys = []
    for i in ind:
        u, v = i_uv[i]
        gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
    form.attributes['indset'] = gkeys

    for i in range(form.number_of_vertices()):
        key = i_k[i]
        form.vertex_attribute(key=key, name='z', value=float(z[i]))

    for c, qi in enumerate(list(q_.ravel())):
        u, v = i_uv[c]
        form.edge_attribute((u, v), 'q', float(qi))

    lp = 0
    for u, v in form.edges_where({'_is_edge': True}):
        if form.edge_attribute((u, v), 'is_symmetry') is False:
            qi = form.edge_attribute((u, v), 'q')
            li = form.edge_length(u, v)
            lp += abs(qi) * li**2
    form.attributes['loadpath'] = float(lp)

    # form.attributes['iter'] = niter
    optimiser.exitflag = exitflag
    optimiser.fopt = fopt
    analysis.form = form
    reactions(form, plot=plot)

    if printout or summary:
        print('\n' + '-' * 50)
        print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
        print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
        print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
        print('fopt      : {0:.3f}'.format(fopt))
        print('-' * 50 + '\n')

    return analysis


def _slsqp(fn, qid0, bounds, fprime, fprime_ieqcons, printout, fieq, args):
    pout = 2 if printout else 0
    opt = fmin_slsqp(fn, qid0, args=args, disp=pout, fprime=fprime, f_ieqcons=fieq, fprime_ieqcons=fprime_ieqcons, bounds=bounds, full_output=1, iter=500)

    return opt[1], opt[0], opt[3], opt[2], opt[4]


def _shgo(fn, bounds, printout, dict_constr, args):

    res = shgo(fn, bounds, args=args, constraints=dict_constr, n=3, iters=3, options={'disp': True})

    return res


def _shgo_constraint_wrapper(x, *args_constr):
    length = len(args_constr)
    args = args_constr[:length-2]
    i = args_constr[length-2]
    fconstr = args_constr[length-1]
    gi = fconstr(x, *args)[i]

    return gi


def _cobyla(fn, qid0, bounds, printout, fieq, args):

    # pout = 2 if printout else 0
    # opt  = fmin_cobyla(fn, qid0, args=args, disp=pout, bounds=bounds, full_output=1, iter=500, f_ieqcons=fieq)

    # W.I.P

    return None


def _diff_evo(fn, bounds, population, generations, printout, plot, frange, args):

    return devo_numpy(fn=fn, bounds=bounds, population=population, generations=generations, printout=printout,
                      plot=plot, frange=frange, args=args)


def _ga(fn, fit_type, num_var, boundaries, num_gen, num_pop, args):

    return ga(fit_function=fn, fit_type=fit_type, num_var=num_var, boundaries=boundaries, num_gen=num_gen, num_pop=num_pop, fargs=args)
