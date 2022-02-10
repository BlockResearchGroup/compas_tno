from compas_tno.solvers import mma_numpy
from numpy import hstack
from numpy import array
from numpy import zeros

from .post_process import post_process_general

import time

from compas_tno.problems import sensitivities_wrapper
from compas_tno.problems import constr_wrapper

try:
    from torch import tensor

    from compas_tno.autodiff.equilibrium_pytorch import f_constraints_pytorch_MMA
    from compas_tno.autodiff.equilibrium_pytorch import f_objective_pytorch
    from compas_tno.autodiff.equilibrium_pytorch import compute_autograd
    from compas_tno.autodiff.equilibrium_pytorch import compute_autograd_jacobian
except BaseException:
    pass


def run_optimisation_MMA(analysis):
    """ Run convex optimisation problem with MMA.

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    obj : analysis
        Analysis object optimised.

    """

    optimiser = analysis.optimiser
    solver = optimiser.settings['solver']
    constraints = optimiser.settings['constraints']
    objective = optimiser.settings['objective']
    gradient = optimiser.settings.get('gradient', False)
    jacobian = optimiser.settings.get('jacobian', False)
    fconstr = optimiser.fconstr
    args = optimiser.args

    bounds = optimiser.bounds
    x0 = optimiser.x0
    g0 = optimiser.g0

    if solver != 'MMA':
        print('Error, Only MMA solver is available for this library!')

    fobj = optimiser.fobj
    fgrad = optimiser.fgrad
    fjac = optimiser.fjac
    fconstr = optimiser.fconstr

    start_time = time.time()

    if gradient and jacobian:

        import nlopt

        # opt = nlopt.opt(nlopt.LD_MMA, len(x0))
        opt = nlopt.opt(nlopt.LD_SLSQP, len(x0))
        fobj_nlopt = create_nlopt_fobj(fobj, fgrad, *args)
        fconstr_nlopt = create_nlopt_fconstr(fconstr, fjac, *args)

        opt.set_min_objective(fobj_nlopt)
        lower = array([lw[0] for lw in bounds]).flatten()
        upper = array([up[1] for up in bounds]).flatten()
        opt.set_lower_bounds(lower)
        opt.set_upper_bounds(upper)
        opt.add_inequality_mconstraint(fconstr_nlopt, zeros((len(g0),)))
        opt.set_maxtime(1000)
        opt.set_xtol_rel(1e-8)

        xopt = opt.optimize(x0.flatten())
        fopt = opt.last_optimum_value()
        result = opt.last_optimize_result()
        if result > 0:
            exitflag = 0

        # args_MMA = list(args)
        # args_MMA.append([fobj, fgrad])
        # f_g_eval = analytical_f_g
        # f_g_df_dg_eval = analytical_f_g_df_dg
    else:

        (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints,
         cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints, max_rol_rx, max_rol_ry, Asym) = args[:48]

        args_MMA = list(args)
        args_MMA.append([fobj, fconstr])
        f_g_eval = pytorch_f_g
        f_g_df_dg_eval = pytorch_f_g_df_dg
        EdinvEi = Edinv*Ei
        Edinv_p = Edinv.dot(p)

        EdinvEi_th = tensor(EdinvEi)
        Edinv_p_th = tensor(Edinv_p)
        C_th = tensor(C.toarray())
        Ci_th = tensor(Ci.toarray())
        Cit_th = Ci_th.t()
        Cf_th = tensor(Cf.toarray())
        pzfree = tensor(pz[free])
        xyz = tensor(hstack([x, y, z]))
        xy = tensor(hstack([x, y]))
        pfixed = tensor(hstack([px, py, pz])[fixed])
        U_th = tensor(U.toarray())
        V_th = tensor(V.toarray())

        args_obj = (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, objective)
        args_constr = (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed,
                       ub, lb, ub_ind, lb_ind, b, constraints, max_rol_rx, max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th)

        args_MMA = [args_obj, args_constr]
        fopt, xopt, exitflag = mma_numpy(f_g_eval, f_g_df_dg_eval, x0, bounds, args_MMA, 10e-4, 100)

    elapsed_time = time.time() - start_time
    optimiser.exitflag = exitflag
    optimiser.time = elapsed_time
    optimiser.fopt = fopt
    optimiser.xopt = xopt
    optimiser.niter = None  # Did not find a way to display number of iteratinos
    if exitflag == 0:
        optimiser.message = 'Optimisation succesful terminated'  # Did not find a way to display message
    else:
        optimiser.message = 'Optimal solution not found'

    post_process_general(analysis)

    return analysis


# ---------------------------------------------------------------
# Calulation with pytorch
# ---------------------------------------------------------------


def pytorch_f_g(variables, *args):
    variables = tensor(variables, requires_grad=False)
    args_obj = args[0]
    args_constr = args[1]
    f0val = array(f_objective_pytorch(variables, *args_obj))
    fval = - 1 * array(f_constraints_pytorch_MMA(variables, *args_constr))
    return f0val, fval


def pytorch_f_g_df_dg(variables, *args):
    variables = tensor(variables, requires_grad=True)
    args_obj = args[0]
    args_constr = args[1]
    fval_tensor = -1 * f_constraints_pytorch_MMA(variables, *args_constr)
    dfdx = array(compute_autograd_jacobian(variables, fval_tensor))
    f0val_tensor = f_objective_pytorch(variables, *args_obj)
    df0dx = array(compute_autograd(variables, f0val_tensor))
    f0val = f0val_tensor.detach().numpy()
    fval = fval_tensor.detach().numpy()
    return f0val, df0dx, fval, dfdx


# ---------------------------------------------------------------
# Calculation analytical
# ---------------------------------------------------------------


def analytical_f_g(xopt, *args):
    [fobj, _] = args[-1]
    args = args[:-1]
    f0val = fobj(xopt, *args)
    fval = - 1 * constr_wrapper(xopt, *args)  # change args here
    return f0val, fval.reshape(-1, 1)


def analytical_f_g_df_dg(xopt, *args):
    [fobj, fgrad] = args[-1]
    args = args[:-1]
    f0val = fobj(xopt, *args)
    fval = - 1 * constr_wrapper(xopt, *args)  # change args here
    df0dx = fgrad(xopt, *args)
    dfdx = -1 * sensitivities_wrapper(xopt, *args)
    return f0val, df0dx, fval.reshape(-1, 1), dfdx


# ---------------------------------------------------------------
# Calculation with NLOPT
# ---------------------------------------------------------------


def create_nlopt_fobj(fobj, fgrad, *args):
    def fobj_nlopt(x, grad):
        f = fobj(x, *args)
        if grad.size > 0:
            grad[:] = array(fgrad(x, *args)).flatten()
        return f
    return fobj_nlopt


def create_nlopt_fconstr(fconstr, fjac, *args):
    def fconstr_nlopt(results, x, grad):
        results[:] = -1 * fconstr(x, *args).flatten()
        if grad.size > 0:
            grad[:] = -1 * fjac(x, *args)
        return
    return fconstr_nlopt
