
from numpy import dot
from numpy import isnan
from numpy import newaxis
from numpy import hstack
from numpy import zeros
from numpy import array

from compas_thrust.algorithms import zlq_from_qid
from compas.numerical import normrow
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import distance_point_point_xy
from compas.geometry import norm_vector

from scipy.sparse import diags


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'f_min_loadpath',
    'f_min_loadpath_pen',
    'f_min_thrust',
    'f_min_thrust_pen',
    'f_max_thrust',
    'f_target',
    'f_constant',
]


def f_min_loadpath(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt.reshape(-1,1), args)
    f = dot(abs(q.transpose()), l2)

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f

def f_min_loadpath_pen(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)
    f = dot(abs(q.transpose()), l2)

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    # Reactions and additional edges

    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xyz = hstack([x,y,z])
    R = CfQC.dot(xyz)
    zR = zeros((len(x),1))
    xR = zeros((len(x),1))
    xR[fixed] = x[fixed] + R[:,0].reshape(-1,1)
    zR[fixed] = z[fixed] + R[:,2].reshape(-1,1)

    for i, elements in joints.items():
        p1, p2, edges = elements
        joint_int = 500.0
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
                    joint_int = offset * 10
                    break
        f+= joint_int

    return f

def f_min_thrust(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)

    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xy = hstack([x,y])
    Rh = CfQC.dot(xy)
    f = sum(normrow(Rh))   # Updated this, there was a mistake before

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f

def f_min_thrust_pen(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]
    xopt = array(xopt)
    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)
        # Verify if it is really necessary calculate z at this step. It may be necessary to only use the dependents equation
    
    # Reactions and additional edges

    CfQ = Cf.transpose().dot(diags(q.flatten()))
    f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])

    xyz = hstack([x,y,z])
    CfQC = CfQ.dot(C)
    R = CfQC.dot(xyz)
    # f = norm_vector(R[:,0]) + norm_vector(R[:,1])
    zR = zeros((len(x),1))
    xR = zeros((len(x),1))
    xR[fixed] = x[fixed] + R[:,0].reshape(-1,1)
    zR[fixed] = z[fixed] + R[:,2].reshape(-1,1)

    for i, elements in joints.items():
        p1, p2, edges = elements
        joint_int = 5000.0
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
                    joint_int = (1 +offset) **2 * 10
                    break
        f+= joint_int

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f

def f_max_thrust(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)
        # Verify if it is really necessary calculate z at this step. It may be necessary to only use the dependents equation
    
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xy = hstack([x,y])
    Rh = CfQC.dot(xy)
    f = (-1) * sum(normrow(Rh))

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f

def f_target(args):

    f = 0

    return f

def f_constant(xopt, *args):

    f = 1.0

    return f