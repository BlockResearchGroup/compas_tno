from compas_tno.solvers import mma_numpy
from numpy import hstack
from numpy import array
from numpy import zeros
from numpy import ones
import numpy as np

from compas_tno.algorithms import f_constraints_pytorch_MMA
from compas_tno.algorithms import f_objective_pytorch
from compas_tno.algorithms import compute_autograd
from compas_tno.algorithms import compute_autograd_jacobian
from compas_tno.algorithms import reactions
from compas_tno.algorithms import zlq_from_qid

from compas.utilities import geometric_key

from .post_process import post_process_general

import time

from compas_tno.problems import sensitivities_wrapper_inequalities
from compas_tno.problems import constr_wrapper
from compas_tno.problems import gradient_fmin
from compas_tno.problems import gradient_fmax
from compas_tno.problems import f_min_thrust
from compas_tno.problems import f_max_thrust

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'run_optimisation_MMA'
]


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

    form = analysis.form
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

        q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints, max_rol_rx, max_rol_ry, Asym = args[
        :48]

        from torch import tensor

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

    post_process_analysis(analysis)

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
    dfdx = array(compute_jacobian(variables, fval_tensor))
    f0val_tensor = f_objective_pytorch(variables, *args_obj)
    df0dx = array(compute_grad(variables, f0val_tensor))
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
    fval = - 1 * constr_wrapper(xopt, *args)
    return f0val, fval.reshape(-1, 1)


def analytical_f_g_df_dg(xopt, *args):
    [fobj, fgrad] = args[-1]
    args = args[:-1]
    f0val = fobj(xopt, *args)
    fval = - 1 * constr_wrapper(xopt, *args)
    df0dx = fgrad(xopt, *args)
    dfdx = -1 * sensitivities_wrapper_inequalities(xopt, *args)
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


# # ----------------------------------------------------------------
# # Helpers Calculation of Derivatives and function evaluation
# # ----------------------------------------------------------------

# def brute_f_g_eval(x0, *args):  # evaluation of objective and of constraints # Old beam1

#     fobj, fconstr = args[-1]
#     args = args[:-1]

#     f0val = fobj(x0, *args)
#     fval = -1 * fconstr(x0, *args).reshape(-1, 1)

#     return f0val, fval


# def brute_f_g_df_dg_eval(x0, *args):  # Derivatives of the objective and of constraints

#     fobj, fconstr = args[-1]
#     args = args[:-1]

#     eps = 1e-6

#     f0val = fobj(x0, *args)
#     fval = -1 * fconstr(x0, *args).reshape(-1, 1)

#     df0dx = d_fobj(fobj, x0, eps, *args)
#     dfdx = -1 * d_fconstr(fconstr, x0, eps, *args)

#     return f0val, df0dx, fval, dfdx

# # ----------------------------------------------------------------
# # Being more efficient in the gradient calculation
# # ----------------------------------------------------------------

# def reduced_f_g_eval(x0, *args):  # evaluation of objective and of constraints

#     variables, dict_constr, objective = args[-1]
#     args = args[:-1]
#     q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin = args[:44]

#     qid, z[fixed] = x0[:k], x0[k:].reshape(-1, 1)
#     q[ind] = array(qid).reshape(-1, 1)
#     q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
#     z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

#     CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
#     xy = hstack([x, y])
#     Rh = CfQC.dot(xy)
#     if objective == 'min':
#         f0val = sum(normrow(Rh))
#     elif objective == 'max':
#         f0val *= -1
#     else:
#         raise NotImplementedError

#     constraints = zeros([0,1])

#     if 'funicular' in dict_constr:
#         constraints = vstack([constraints, (q.ravel() - qmin).reshape(-1, 1)])  # >= 0
#     if 'envelope' in dict_constr:
#         constraints = vstack([constraints, ub - z[ub_ind], z[lb_ind] - lb])  # >= 0
#     if 'reac_bounds' in dict_constr:
#         CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
#         xyz = hstack([x, y, z])
#         p_fixed = hstack([px, py, pz])[fixed]
#         R = CfQC.dot(xyz) - p_fixed
#         Rx = abs(b[:, 0].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 0], R[:, 2]).reshape(-1, 1)))
#         Ry = abs(b[:, 1].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 1], R[:, 2]).reshape(-1, 1)))
#         constraints = vstack([constraints, Rx, Ry])
#     if 'cracks' in dict_constr:
#         crack_tol = 10e-4
#         constraints = vstack([constraints, (lb[cracks_lb] - z[cracks_lb]) + crack_tol, (z[cracks_ub] - ub[cracks_ub]) + crack_tol])

#     fval = -1 * transpose(constraints)[0].reshape(-1,1)

#     # fval = -1 * transpose(vstack([qpos, upper_limit, lower_limit]))[0].reshape(-1,1)

#     return f0val, fval

# def reduced_f_g_df_dg_eval(x0, *args): # Derivatives of the objective and of constraints

#     variables, dict_constr, objective = args[-1]
#     args = args[:-1]
#     q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin = args[:44]

#     qid, z[fixed] = x0[:k], x0[k:].reshape(-1, 1)
#     q[ind] = array(qid).reshape(-1, 1)
#     q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
#     z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

#     z0 = z.copy()
#     qd0 = q[dep].copy()

#     eps = 1e-6

#     CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
#     xy = hstack([x, y])
#     Rh = CfQC.dot(xy)
#     if objective == 'min':
#         f0val = sum(normrow(Rh))
#     elif objective == 'max':
#         f0val *= -1
#     else:
#         raise NotImplementedError

#     constraints = zeros([0,1])

#     if 'funicular' in dict_constr:
#         constraints = vstack([constraints, (q.ravel() - qmin).reshape(-1, 1)])  # >= 0
#     if 'envelope' in dict_constr:
#         constraints = vstack([constraints, ub - z[ub_ind], z[lb_ind] - lb])  # >= 0
#     if 'reac_bounds' in dict_constr:
#         CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
#         xyz = hstack([x, y, z])
#         p_fixed = hstack([px, py, pz])[fixed]
#         R = CfQC.dot(xyz) - p_fixed
#         Rx0 = abs(b[:, 0].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 0], R[:, 2]).reshape(-1, 1)))  # >= 0
#         Ry0 = abs(b[:, 1].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 1], R[:, 2]).reshape(-1, 1)))  # >= 0
#         constraints = vstack([constraints, Rx0, Ry0])
#     if 'cracks' in dict_constr:
#         crack_tol = 10e-4
#         constraints = vstack([constraints, (lb[cracks_lb] - z[cracks_lb]) + crack_tol, (z[cracks_ub] - ub[cracks_ub]) + crack_tol])

#     fval = -1 * transpose(constraints)[0].reshape(-1,1)

#     # fval = -1 * transpose(vstack([qpos, upper_limit, lower_limit]))[0].reshape(-1,1)

#     m = len(fval)
#     n = len(x0)

#     # -------- Derivatives

#     df0dx = zeros((n, 1))  # Derivatives of the objective function with regards to each variable in a vector
#     dfdx = zeros((m, n))  # Derivatives of the m constraints with regards to each variable in a matrix

#     m_dep = len(dep)
#     m_envelope = len(ub_ind)
#     i_zb = len(fixed)
#     p_fixed = hstack([px, py, pz])[fixed]

#     print('0 to ', m_dep, 'qpos //', m_dep, 'to', m_dep + 2*m_envelope, 'envelope //', 2*i_zb, 'additional')

#     for j in range(n):
#         diff = zeros((n, 1))
#         diff[j] = eps
#         x_ = x0 + diff

#         qid, z[fixed] = x_[:k], x_[k:].reshape(-1, 1)
#         q[ind] = array(qid).reshape(-1, 1)
#         q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
#         z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))
#         CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)

#         if 'funicular' in dict_constr:
#             dfdx[:m_dep, j] = -1 * ((q[dep] - qd0)/eps).reshape(-1,)
#         if 'envelope' in dict_constr:
#             dfdx[m_dep:m_dep + m_envelope, j] = ((z - z0)/eps).reshape(-1,)
#             dfdx[m_dep + m_envelope:m_dep + 2*m_envelope, j] = -1 * ((z - z0)/eps).reshape(-1,)
#         if 'reac_bounds' in dict_constr:
#             xyz = hstack([x, y, z])
#             R = CfQC.dot(xyz) - p_fixed
#             dRx = (abs(b[:, 0].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 0], R[:, 2]).reshape(-1, 1))) - Rx0)/eps
#             dRy = (abs(b[:, 1].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 1], R[:, 2]).reshape(-1, 1))) - Ry0)/eps
#             dfdx_R = -1 * vstack([dRx, dRy]).reshape(-1,)
#             dfdx[-2*i_zb:, j] = dfdx_R
#         if 'cracks' in dict_constr:
#             NotImplementedError

#         f = sum(normrow(CfQC.dot(xy)))
#         df0dx[j, 0] = (f - f0val)/eps

#     if objective == 'min':
#         pass
#     elif objective == 'max':
#         f0val *= -1
#         df0dx *= -1
#     else:
#         raise NotImplementedError

#     return f0val, df0dx, fval, dfdx

# # ----------------------------------------------------------------
# # Being analytic efficient in the gradient calculation
# # ----------------------------------------------------------------

# # def beam2(x0, fobj, fconstr, *args): # Derivatives of the objective and of constraints

# #     q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

# #     qid, z[fixed] = x0[:k], x0[k:].reshape(-1, 1)
# #     q[ind] = array(qid).reshape(-1, 1)
# #     q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
# #     z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

# #     z0 = z.copy()

# #     eps = 1e-6

# #     f0val = fobj(x0,*args)
# #     fval = -1 * fconstr(x0, *args).reshape(-1,1)

# #     m = len(fval)
# #     n = len(x0)

# #     # -------- Derivatives

# #     df0dx = d_fobj(fobj, x0, eps, *args)
# #     dfdx = zeros((m,n))

# #     m_dep = len(dep)
# #     B = Edinv.dot(Ei)
# #     B_ = zeros((m_dep,len(fixed)))
# #     qpos_contribution = hstack([B, B_])
# #     dfdx[list(range(m_dep))] = -1 * qpos_contribution

# #     m_envelope = len(ub_ind)

# #     for i in range(m_envelope):
# #         for j in range(n):
# #             diff = zeros((n, 1))
# #             diff[j] = eps
# #             x_ = x0 + diff

# #             qid, z[fixed] = x_[:k], x_[k:].reshape(-1, 1)
# #             q[ind] = array(qid).reshape(-1, 1)
# #             q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
# #             z[free, 0] = spsolve(Cit.dot(diags(q.flatten())).dot(Ci), pz[free] - Cit.dot(diags(q.flatten())).dot(Cf).dot(z[fixed]))

# #             dfdx[m_dep + i, j] =  -1 * (-z[i] - (-z0[i]))/diff[j]
# #             dfdx[m_dep + i + m_envelope, j] = - dfdx[i, j]

# #     # dfdx = -1 * d_fconstr(fconstr, x0, eps, *args)

# #     # upper_limit = z[ub_ind] - ub  # <= 0
# #     # lower_limit = lb - z[lb_ind]  # <= 0

#     return f0val, df0dx, fval, dfdx
