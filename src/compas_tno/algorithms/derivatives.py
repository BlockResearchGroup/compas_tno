
from numpy import dot
from numpy import isnan
from numpy import newaxis
from numpy import hstack
from numpy import zeros
from numpy import array

from compas_tno.algorithms.equilibrium import zlq_from_qid
from compas.numerical import normrow
from compas.geometry import is_intersection_segment_segment_xy
from compas.geometry import intersection_line_segment_xy
from compas.geometry import distance_point_point_xy
from compas.geometry import norm_vector

from scipy.sparse import diags


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'd_fobj',
    'd_fconstr',
]


def d_fobj(fobj, x0, eps, *args):
    f0val = fobj(x0,*args)
    n = len(x0)
    df0dx = zeros((n,1))
    for i in range(n):
        diff = zeros((n,1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    return df0dx

def d_fconstr(fconstr, x0, eps, *args):
    fval = fconstr(x0, *args).reshape(-1,1)
    m = len(fval)
    n = len(x0)
    dfdx = zeros((m, n))
    for i in range(m):
        for j in range(n):
            diff = zeros((n, 1))
            diff[j] = eps
            dfdx[i, j] = (fconstr(x0 + diff, *args).reshape(-1,1) - fval)[i]/diff[j]
    return dfdx


def d_min_thrust(fobj, x0, eps, *args):
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    f0val = fobj(x0,*args)
    n = len(x0)
    df0dx = zeros((n,1))
    for i in range(n):
        diff = zeros((n,1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    return df0dx


def d_f_ub_lb(fobj, x0, eps, *args):
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    f0val = fobj(x0,*args)
    n = len(x0)
    df0dx = zeros((n,1))
    for i in range(n):
        diff = zeros((n,1))
        diff[i] = eps
        df0dx[i] = (fobj(x0 + diff, *args) - f0val)/diff[i]

    B = Edinv.dot(Ei)
    m_qpos = B.shape[0]
    qpos_contribution = B



    # transpose(vstack([qpos, upper_limit, lower_limit]))[0]


    return df0dx
