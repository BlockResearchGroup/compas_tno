
from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import array
from numpy.linalg import matrix_rank
from numpy import transpose
from compas.geometry import is_intersection_segment_segment_xy
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

    q = q_from_qid(xopt, args)

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

    constr = 0.0
    for i, elements in joints.items():
        p1, p2, edges = elements
        joint_int = False
        for k1, k2 in edges:
            e1 = [x[k1], z[k1]]
            e2 = [x[k2], z[k2]]
            edge_2D = [e1,e2]
            joint_2D = [[p1[0],p1[2]],[p2[0],p2[2]]]
            if is_intersection_segment_segment_xy(joint_2D,edge_2D) == True:
                joint_int = True
                break
        constr+float(joint_int*1)-1
    constr = array(constr)
    # print(constr)

    if any(q) == False:
        return 10**10

    return constr.transpose()