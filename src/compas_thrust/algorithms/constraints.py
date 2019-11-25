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
from compas_thrust.algorithms import zlq_from_qid
from compas_thrust.algorithms import q_from_qid
from scipy.sparse import diags


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'f_ub_lb',
    'f_compression',
    'f_joints'
]


def f_ub_lb(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Constraints on Heights
    upper_limit = ub - z[ub_ind] # >= 0
    lower_limit = z[lb_ind] - lb # >= 0

    tol = 1e-10

    # Constraints on Reactions
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x,y,z])
    R = CfQC.dot(xyz)
    Rx_angle = tol + abs(b[:,0].reshape(-1,1)) - abs(multiply(z[fixed],divide(R[:,0],R[:,2]).reshape(-1,1)))
    Ry_angle = tol + abs(b[:,1].reshape(-1,1)) - abs(multiply(z[fixed],divide(R[:,1],R[:,2]).reshape(-1,1)))

    return transpose(vstack([upper_limit, lower_limit, Rx_angle, Ry_angle]))[0]

def f_compression(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args
    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
    else:
        qid = xopt

    q = q_from_qid(qid, args)

    return hstack([q.ravel() + 10**(-5)])

def f_joints(xopt, *args):
    # Boolean constraint - WIP to make it at least non-linear
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Reactions and additional edges
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x,y,z])
    R = CfQC.dot(xyz)
    zR = zeros((len(x),1))
    xR = zeros((len(x),1))
    xR[fixed] = x[fixed] + R[:,0].reshape(-1,1)
    zR[fixed] = z[fixed] + R[:,2].reshape(-1,1)

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
            edge_2D = [e1,e2]
            joint_2D = [[p1[0],p1[2]],[p2[0],p2[2]]]
            if is_intersection_segment_segment_xy(joint_2D,edge_2D) == True:
                joint_int = 0.0
                break
            else:
                virtual_pt = intersection_line_segment_xy(joint_2D,edge_2D)
                if virtual_pt:
                    offset = min(distance_point_point_xy([virtual_pt[0],virtual_pt[2]],joint_2D[0]),distance_point_point_xy([virtual_pt[0],virtual_pt[2]],joint_2D[1]))
                    joint_int = - (offset +1) ** 2
                    break
        constr+= joint_int

    return constr


def f_cracks(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub = args

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    # Constraints on Heights
    upper_limit = 0 #ub - z[ub_ind] # >= 0
    lower_limit = 0 #z[lb_ind] - lb # >= 0

    tol = 1e-10

    # Constraints on Reactions
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x,y,z])
    R = CfQC.dot(xyz)
    Rx_angle = tol + abs(b[:,0].reshape(-1,1)) - abs(multiply(z[fixed],divide(R[:,0],R[:,2]).reshape(-1,1)))
    Ry_angle = tol + abs(b[:,1].reshape(-1,1)) - abs(multiply(z[fixed],divide(R[:,1],R[:,2]).reshape(-1,1)))

    # Constraints on Cracks

    crack_tol = 0.001

    lower_cracks = - abs(lb[cracks_lb] - z[cracks_lb]) + crack_tol
    upper_cracks = - abs(z[cracks_ub] - ub[cracks_ub]) + crack_tol

    # This can be simplified to not count with the abs... Will test in 3D later on... and see :)

    return transpose(vstack([upper_limit, lower_limit, lower_cracks, upper_cracks, Rx_angle, Ry_angle]))[0]