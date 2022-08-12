from scipy.optimize import fmin_slsqp
from scipy.optimize import shgo

from compas.numerical import devo_numpy
from compas.numerical import ga

from .post_process import post_process_general

import time


def run_optimisation_scipy(analysis):
    """ Run nonlinear optimisation problem with SciPy.

    Parameters
    ----------
    analysis : Analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    analysis : Analysis
        Analysis object optimised.

    """

    optimiser = analysis.optimiser
    solver = optimiser.settings['solver']
    fobj = optimiser.fobj
    fconstr = optimiser.fconstr
    fgrad = optimiser.fgrad
    fjac = optimiser.fjac
    args = [optimiser.M]
    bounds = optimiser.bounds
    x0 = optimiser.x0
    printout = optimiser.settings.get('printout', True)
    grad_choice = optimiser.settings.get('gradient', False)
    jac_choice = optimiser.settings.get('jacobian', False)
    max_iter = optimiser.settings.get('max_iter', 500)
    callback = optimiser.callback

    if grad_choice is False:
        fgrad = None
    if jac_choice is False:
        fjac = None

    start_time = time.time()

    if solver == 'slsqp' or solver == 'SLSQP':
        fopt, xopt, exitflag, niter, message = _slsqp(fobj, x0, bounds, fgrad, fjac, printout, fconstr, args, max_iter, callback)
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
        else:
            print(message)

    elapsed_time = time.time() - start_time
    if printout:
        print('Solving Time: {0:.1f} sec'.format(elapsed_time))

    # Store output info in optimiser

    optimiser.exitflag = exitflag
    optimiser.time = elapsed_time
    optimiser.fopt = float(fopt)
    optimiser.xopt = xopt
    optimiser.niter = niter
    optimiser.message = message

    post_process_general(analysis)

    return analysis


def _slsqp(fn, qid0, bounds, fprime, fprime_ieqcons, printout, fieq, args, iter, callback):
    pout = 2 if printout else 0
    opt = fmin_slsqp(fn, qid0, args=args, disp=pout, fprime=fprime, f_ieqcons=fieq, fprime_ieqcons=fprime_ieqcons, bounds=bounds, full_output=1, iter=iter, callback=callback)

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
