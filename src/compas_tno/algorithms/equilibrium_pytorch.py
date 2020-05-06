from torch import mm
from torch import diagflat
from torch import tensor
from torch import cat
from torch import mul
from torch import div
from torch import solve
from torch import norm
from torch import zeros as thzeros
from torch.autograd.gradcheck import zero_gradients
from torch import float64

from numpy import hstack
from numpy import zeros

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'q_from_variables_pytorch',
    'z_from_variables_pytorch',
    'reac_bound_variables_pytorch',
    'f_min_thrust_pytorch',
    'f_max_thrust_pytorch',
    'f_constraints_pytorch',
    'f_constraints_pytorch_MMA',
    'bounds_constraints_pytorch',
    'f_objective_pytorch',
    'compute_grad',
    'compute_jacobian'
]


def q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep):
    k = len(ind)
    m = k + len(dep)
    q = thzeros(m, 1, dtype=float64)
    q[ind] = variables[:k]
    q[dep] = - Edinv_p_th + mm(EdinvEi_th, q[ind])
    return q


def z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci, Cit, Cf, pzfree, free, fixed):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    zfixed = variables[-Cf.shape[1]:]
    Q = diagflat(q)
    Ai = mm(mm(Cit, Q), Ci)
    Af = mm(mm(Cit, Q), Cf)
    b = pzfree - mm(Af, zfixed)
    zi, LU = solve(b, Ai)
    z = tensor(zeros((len(free)+len(fixed), 1)))
    z[free] = zi
    z[fixed] = zfixed
    return z


def reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, pfixed, xyz):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    zfixed = variables[-Cf.shape[1]:]
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    R = mm(CfQC, xyz) - pfixed
    length_x = abs(mul(zfixed, div(R[:, 0], R[:, 2]).reshape(-1, 1)))
    length_y = abs(mul(zfixed, div(R[:, 1], R[:, 2]).reshape(-1, 1)))
    length = cat([length_x, length_y])
    return length


def f_min_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, xy):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    Rh = mm(CfQC, xy)
    R = norm(Rh, dim=1)
    f = sum(R)
    return f


def f_max_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C, Cf, xy):
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    Rh = mm(CfQC, xy)
    R = norm(Rh, dim=1)
    f = -1 * sum(R)
    return f


def f_loadpath_pytorch(variables):

    return


def f_target_pytorch(variables):

    return


def f_objective_pytorch(variables, *args):
    Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, objective = args
    if objective == 'min':
        f = f_min_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, xy)
    if objective == 'max':
        f = f_max_thrust_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, xy)
    if objective == 'loadpath':
        pass
    if objective == 'target':
        pass
    return f


def f_constraints_pytorch(variables, *args):
    Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr, max_rol_rx, max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th = args
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    if 'funicular' in dict_constr:
        constraints = q
    if 'envelope' in dict_constr:
        z = z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci_th, Cit_th, Cf_th, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1, 1)) - length
        constraints = cat([constraints, reac_bound])
    if 'cracks' in dict_constr:
        pass
    if 'rollers' in dict_constr:
        Cftx = Cf_th[:, rol_x].transpose()
        Cfty = Cf_th[:, rol_y].transpose()
        rx_check = max_rol_rx - abs(mm(Cftx, mm(U_th, q)) - px[rol_x])
        ry_check = max_rol_ry - abs(mm(Cfty, mm(V_th, q)) - py[rol_y])
        constraints = cat([constraints, rx_check, ry_check])
    if 'symmetry' in dict_constr:
        A_q = mm(tensor(Asym), variables)
        constraints = cat([constraints, A_q])
    return constraints


def f_constraints_pytorch_MMA(variables, *args):  # THe equality constraints are transformed in 2 inequalities (symmetry constraint)
    Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr, max_rol_rx, max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th = args
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    if 'funicular' in dict_constr:
        constraints = q
    if 'envelope' in dict_constr:
        z = z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci_th, Cit_th, Cf_th, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1, 1)) - length
        constraints = cat([constraints, reac_bound])
    if 'cracks' in dict_constr:
        pass
    if 'rollers' in dict_constr:
        Cftx = Cf_th[:, rol_x].transpose()
        Cfty = Cf_th[:, rol_y].transpose()
        rx_check = max_rol_rx - abs(mm(Cftx, mm(U_th, q)) - px[rol_x])
        ry_check = max_rol_ry - abs(mm(Cfty, mm(V_th, q)) - py[rol_y])
        constraints = cat([constraints, rx_check, ry_check])
    if 'symmetry' in dict_constr:
        A_q = mm(tensor(Asym), variables)
        constraints = cat([constraints, A_q, -1 * A_q])
    return constraints


def bounds_constraints_pytorch(variables, *args):
    Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr, max_rol_rx, max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th = args
    q = q_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep)
    if 'funicular' in dict_constr:
        constraints = q
    if 'envelope' in dict_constr:
        z = z_from_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, Ci_th, Cit_th, Cf_th, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_variables_pytorch(variables, Edinv_p_th, EdinvEi_th, ind, dep, C_th, Cf_th, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1, 1)) - length
        constraints = cat([constraints, reac_bound])
    if 'cracks' in dict_constr:
        pass
    if 'rollers' in dict_constr:
        Cftx = Cf_th[:, rol_x].transpose()
        Cfty = Cf_th[:, rol_y].transpose()
        rx_check = max_rol_rx - abs(mm(Cftx, mm(U_th, q)) - px[rol_x])
        ry_check = max_rol_ry - abs(mm(Cfty, mm(V_th, q)) - py[rol_y])
        constraints = cat([constraints, rx_check, ry_check])
    cu = [10e10]*len(constraints)
    cl = [0.0]*len(constraints)
    # End of inequality constraints
    if 'symmetry' in dict_constr:
        A_q = mm(tensor(Asym), variables)
        cu = hstack([cu, [0.0]*len(A_q)])
        cl = hstack([cl, [0.0]*len(A_q)])
    return cu, cl


def compute_grad(variables, f):
    f.backward(retain_graph=True)
    grad = variables.grad.data
    return grad


def compute_jacobian(inputs, outputs):
    d_otp = outputs.size()[0]
    d_inp = inputs.size()[0]
    jacobian = tensor(zeros((d_otp, d_inp)))
    grad_output = tensor(zeros((d_otp, 1)))
    for i in range(d_otp):
        zero_gradients(inputs)
        grad_output.zero_()
        grad_output[i] = 1
        outputs.backward(grad_output, retain_graph=True)
        jacobian[i] = inputs.grad.data.flatten()
    return jacobian
