from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import array
from numpy import zeros
from numpy.linalg import matrix_rank
from numpy import transpose
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import distance_point_point_xy
from compas_tno.algorithms import zlq_from_qid
from compas_tno.algorithms import q_from_qid
from scipy.sparse import diags


__all__ = [
    'f_ub_lb',
    'f_compression',
    'f_joints',
    'f_cracks',
    'f_ub_lb_red'
]


def f_ub_lb(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Constraints on Heights
    upper_limit = ub - z[ub_ind]  # >= 0
    lower_limit = z[lb_ind] - lb  # >= 0
    tol = 1e-10

    # Constraints on Reactions
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x, y, z])
    p_fixed = hstack([px, py, pz])[fixed]
    R = CfQC.dot(xyz) - p_fixed
    Rx_angle = abs(b[:, 0].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 0], R[:, 2]).reshape(-1, 1))) # + tol >= 0
    Ry_angle = abs(b[:, 1].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 1], R[:, 2]).reshape(-1, 1))) # + tol >= 0

    # Positive Qs
    qpos = (q.ravel() + 10**(-5)).reshape(-1, 1)

    max_f_x = array([1.0]*len(rol_x))  # Change to be personalised for each node....
    max_f_y = array([1.0]*len(rol_y))

    # Reactions on the walls
    reac_rol_x = tol #(array(max_f_x) - abs(array(Cftx.dot(U * q.ravel()) - px[rol_x].ravel()))).reshape(-1, 1)
    reac_rol_y = tol# (array(max_f_y) - abs(array(Cfty.dot(V * q.ravel()) - px[rol_y].ravel()))).reshape(-1, 1)

    return transpose(vstack([qpos, upper_limit, lower_limit, Rx_angle, Ry_angle, reac_rol_x, reac_rol_y]))[0]

def f_ub_lb_red(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Constraints on Heights
    upper_limit = ub - z[ub_ind]  # >= 0
    lower_limit = z[lb_ind] - lb  # >= 0

    # Positive Qs
    qpos = (q[dep].ravel() + 10**(-5)).reshape(-1, 1)

    return transpose(vstack([qpos, upper_limit, lower_limit]))[0]

def f_compression(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
    else:
        qid = xopt

    q = q_from_qid(qid, args)

    return hstack([q.ravel() + 10**(-5)])


def f_joints(xopt, *args):
    # Boolean constraint - WIP to make it at least non-linear
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Reactions and additional edges
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x, y, z])
    R = CfQC.dot(xyz)
    zR = zeros((len(x), 1))
    xR = zeros((len(x), 1))
    xR[fixed] = x[fixed] + R[:, 0].reshape(-1, 1)
    zR[fixed] = z[fixed] + R[:, 2].reshape(-1, 1)

    constr = 0.0
    for i, elements in joints.items():
        p1, p2, edges = elements
        joint_int = - 50.0
        for k1, k2 in edges:
            if -k1 == k2:
                e1 = xR[-k1], zR[-k1]
            else:
                e1 = x[k1], z[k1]
            e2 = x[k2], z[k2]
            edge_2D = [e1, e2]
            joint_2D = [[p1[0], p1[2]], [p2[0], p2[2]]]
            if is_intersection_segment_segment_xy(joint_2D, edge_2D) == True:
                joint_int = 0.0
                break
            else:
                virtual_pt = intersection_line_segment_xy(joint_2D, edge_2D)
                if virtual_pt:
                    offset = min(distance_point_point_xy([virtual_pt[0], virtual_pt[2]], joint_2D[0]), distance_point_point_xy([virtual_pt[0], virtual_pt[2]], joint_2D[1]))
                    joint_int = - (offset + 1) ** 2
                    break
        constr += joint_int

    return constr


def f_cracks(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Compression only

    qpos = (q.ravel() + 10**(-5)).reshape(-1, 1)

    # Constraints on Heights
    upper_limit = ub - z[ub_ind]  # >= 0
    lower_limit = z[lb_ind] - lb  # >= 0

    tol = 1e-10

    # Constraints on Reactions
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x, y, z])
    R = CfQC.dot(xyz)
    Rx_angle = tol + abs(b[:, 0].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 0], R[:, 2]).reshape(-1, 1)))
    Ry_angle = tol + abs(b[:, 1].reshape(-1, 1)) - abs(multiply(z[fixed], divide(R[:, 1], R[:, 2]).reshape(-1, 1)))

    # Constraints on Cracks

    crack_tol = 0.030

    lower_cracks = (lb[cracks_lb] - z[cracks_lb]) + crack_tol
    upper_cracks = (z[cracks_ub] - ub[cracks_ub]) + crack_tol

    # Constraints on flying reactions

    max_f_x = array([5.0]*len(rol_x))
    max_f_y = array([5.0]*len(rol_y))

    reac_rol_x = (array(max_f_x) - abs(array(Cftx.dot(U * q.ravel()) - px[rol_x].ravel()))).reshape(-1, 1)
    reac_rol_y = (array(max_f_y) - abs(array(Cfty.dot(V * q.ravel()) - px[rol_y].ravel()))).reshape(-1, 1)

    return transpose(vstack([qpos, upper_limit, lower_limit, lower_cracks, upper_cracks, Rx_angle, Ry_angle, reac_rol_x, reac_rol_y]))[0]
