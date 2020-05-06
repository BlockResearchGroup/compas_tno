from compas_tno.solvers import mma_numpy
from torch import tensor

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
from compas.numerical import normrow

from numpy import identity
from numpy import hstack
from numpy import vstack
from numpy import transpose
from numpy import array
from numpy import zeros

from compas_tno.algorithms.equilibrium_pytorch import f_constraints_pytorch_MMA
from compas_tno.algorithms.equilibrium_pytorch import f_objective_pytorch
from compas_tno.algorithms.equilibrium_pytorch import compute_grad
from compas_tno.algorithms.equilibrium_pytorch import compute_jacobian

import logging

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import d_fobj
from compas_tno.algorithms import d_fconstr

from compas_tno.algorithms import zlq_from_qid
from compas.utilities import geometric_key

from numpy import multiply
from numpy import divide

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
    solver = optimiser.data['solver']
    variables = optimiser.data['variables']
    constraints = optimiser.data['constraints']
    objective = optimiser.data['objective']
    fobj = optimiser.fobj
    fconstr = optimiser.fconstr
    args = optimiser.args
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints, max_rol_rx, max_rol_ry, Asym = args[:48]
    i_uv = form.index_uv()
    i_k = form.index_key()
    k_i = form.key_index()
    bounds = optimiser.bounds
    x0 = optimiser.x0
    plot = optimiser.data['plot']

    if solver != 'MMA':
        print('Error, Only MMA solver is available for this library!')

    args_MMA = list(args)
    args_MMA.append([fobj, fconstr])
    f_g_eval = pytorch_f_g
    f_g_df_dg_eval = pytorch_f_g_df_dg

    # Tensor Transformation

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
    print('This is exitflag', exitflag)

    if 'ind' in variables and 'zb' in variables:
        q[ind] = xopt[:k].reshape(-1, 1)
        z[fixed] = xopt[k:].reshape(-1, 1)
    else:
        q[ind] = xopt[:k].reshape(-1, 1)

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
    form.attributes['loadpath'] = lp

    exitflag = 0
    optimiser.exitflag = exitflag
    optimiser.fopt = fopt
    analysis.form = form
    reactions(form)

    print('\n' + '-' * 50)
    print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
    print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
    print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
    print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
    print('fopt      : {0:.3f}'.format(fopt))
    print('-' * 50 + '\n')

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
