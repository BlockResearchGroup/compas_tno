
from numpy import dot
from numpy import isnan
from numpy import newaxis

from compas_thrust.algorithms import zlq_from_qid

from scipy.sparse import diags


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'f_min_loadpath',
    'f_min_thrust',
    'f_max_thrust',
    'f_target',
    'f_constant',
]


def f_min_loadpath(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)
    f = dot(abs(q.transpose()), l2)

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f


def f_min_thrust(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)
        # Verify if it is really necessary calculate z at this step. It may be necessary to only use the dependents equation
    
    CfQ = Cf.transpose().dot(diags(q.flatten()))
    f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f

def f_max_thrust(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args

    if len(xopt)>k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1,1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt, args)
        # Verify if it is really necessary calculate z at this step. It may be necessary to only use the dependents equation
    
    CfQ = Cf.transpose().dot(diags(q.flatten()))
    f = -1* (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])

    if isnan(f) == True or any(xopt) == False:
        return 10**10

    return f

def f_target(args):

    f = 0

    return f

def f_constant(xopt, *args):

    f = 1

    return f