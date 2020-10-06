
from numpy import dot
from numpy import isnan
from numpy import hstack
from numpy import array

from compas_tno.algorithms.equilibrium import zlq_from_qid
from compas_tno.algorithms.equilibrium import zq_from_qid
from compas.numerical import normrow

from scipy.sparse import diags


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'f_min_loadpath',
    'f_min_thrust',
    'f_max_thrust',
    'f_target',
    'f_constant',
]


def f_min_loadpath(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, l2, q, _ = zlq_from_qid(qid, args)
    else:
        z, l2, q, _ = zlq_from_qid(xopt.reshape(-1, 1), args)
    f = dot(abs(q.transpose()), l2)

    if isnan(f) is True or any(xopt) is False:
        return 10**10

    return f


def f_min_thrust(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, q = zq_from_qid(qid, args)
    else:
        z, q = zq_from_qid(qid, args)

    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xy = hstack([x, y])
    Rh = CfQC.dot(xy) - hstack([px, py])[fixed]
    f = sum(normrow(Rh))

    if isnan(f) is True or any(xopt) is False:
        return 10**10

    return f


def f_max_thrust(xopt, *args):

    f = (-1) * f_min_thrust(xopt, *args)

    if isnan(f) is True or any(xopt) is False:
        return 10**10

    return f


def f_target(xopt, *args):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]

    if len(xopt) > k:
        qid, z[fixed] = xopt[:k], xopt[k:].reshape(-1, 1)
        args = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
        z, q = zq_from_qid(qid, args)
    else:
        z, q = zq_from_qid(qid, args)

    f = sum((z - s)**2)

    return f


def f_constant(xopt, *args):

    f = 1.0

    return f


def f_reduce_thk(xopt, *args):

    return xopt[-1]


def f_min_thrust_pytorch(xopt, *args):
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i = args[:35]
    qid = xopt[:k].reshape(-1, 1)
    q[ind] = array(qid).reshape(-1, 1)
    q[dep] = Edinv.dot(- p + Ei.dot(q[ind]))
    CfQC = Cf.transpose().dot(diags(q.flatten())).dot(C)
    xy = hstack([x, y])
    Rh = CfQC.dot(xy) - hstack([px, py])[fixed]
    f = sum(normrow(Rh))
    return f


def f_tight_crosssection(xopt, *args):

    return -1 * xopt[-1]
