from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import array
from numpy import zeros
from numpy import transpose
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import distance_point_point_xy
from compas_tno.algorithms.equilibrium import zlq_from_qid
from compas_tno.algorithms.equilibrium import zq_from_qid
from compas_tno.algorithms.equilibrium import q_from_qid
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from compas_tno.problems.bounds_update import ub_lb_update
from compas_tno.problems.bounds_update import b_update


__all__ = [
    'constr_wrapper',
    'constr_wrapper_inequalities',
]


# def constr_wrapper(xopt, *args):

#     q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym = args

#     if len(xopt) > k:
#         qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
#         args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
#         z, q = zq_from_qid(qid, args)
#     else:
#         z, q = zq_from_qid(qid, args)

#     constraints = zeros([0, 1])

#     if 'funicular' in dict_constr:
#         constraints = vstack([constraints, (q[dep] - qmin).reshape(-1, 1)])  # >= 0
#     if 'envelope' in dict_constr:
#         constraints = vstack([constraints, ub - z[ub_ind], z[lb_ind] - lb])  # >= 0
#     if 'reac_bounds' in dict_constr:
#         CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
#         xyz = hstack([x, y, z])
#         p_fixed = hstack([px, py, pz])[fixed]
#         R = CfQC.dot(xyz) - p_fixed
#         Rx = abs(b[:, 0].reshape(-1, 1)) - multiply(z[fixed] - s[fixed], abs(divide(R[:, 0], R[:, 2]).reshape(-1, 1)))  # >= 0 absoluute numpy
#         Ry = abs(b[:, 1].reshape(-1, 1)) - multiply(z[fixed] - s[fixed], abs(divide(R[:, 1], R[:, 2]).reshape(-1, 1)))  # >= 0
#         constraints = vstack([constraints, Rx, Ry])
#     if 'cracks' in dict_constr:
#         crack_tol = 10e-4
#         constraints = vstack([constraints, (lb[cracks_lb] - z[cracks_lb]) + crack_tol, (z[cracks_ub] - ub[cracks_ub]) + crack_tol])
#     if 'rollers' in dict_constr:
#         rx_check = max_rol_rx - abs(Cftx.dot(U.dot(q)) - px[rol_x])
#         ry_check = max_rol_ry - abs(Cfty.dot(V.dot(q)) - py[rol_y])
#         constraints = vstack([constraints, rx_check, ry_check])
#     if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
#         A_q = Asym.dot(vstack([q[ind], z[fixed]]))
#         constraints = vstack([constraints, A_q, -1 * A_q])

#     return transpose(constraints)[0]  # .reshape(-1, 1)


def constr_wrapper(xopt, *args):  # change name shape data to shape

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape = args[:50]

    if 'ind' in variables:  # Not yet deal with all-q
        q[ind] = xopt[:k].reshape(-1, 1)
    if 'zb' in variables:
        z[fixed] = xopt[k:k+len(fixed)].reshape(-1, 1)
    if 't' in variables or 's' in variables or 'n' in variables:
        thk = xopt[-1].item()
        t = shape.data['t']

    args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
    z, q = zq_from_qid(q[ind], args)

    constraints = zeros([0, 1])

    if 'funicular' in dict_constr:
        constraints = vstack([constraints, (q[dep] - qmin).reshape(-1, 1)])  # >= 0
    if 'envelope' in dict_constr:
        if 't' in variables or 's' in variables or 'n' in variables:
            ub, lb = ub_lb_update(x, y, thk, t, shape, ub, lb, variables)
        constraints = vstack([constraints, ub[ub_ind] - z[ub_ind], z[lb_ind] - lb[lb_ind]])
    if 'reac_bounds' in dict_constr:
        if 't' in variables:
            b = b_update(x, y, thk, fixed, shape)
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


def constr_wrapper_inequalities(xopt, *args):  # This considers a equality in Asym.
    """
    This computes the sensitivities considering only inequality constraints.
    """
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, dict_constr, max_rol_rx, max_rol_ry, Asym, variables, shape = args[:50]
    constraints = constr_wrapper(xopt, *args)
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in dict_constr):
        constraints = vstack([constraints, -Asym])

    return constraints
