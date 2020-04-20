# from torch import mm
# from torch import diagflat
# from torch import tensor
# from torch import cat
# from torch import mul
# from torch import div
# from torch import solve
# from torch import norm
# from torch.autograd.gradcheck import zero_gradients

from numpy import zeros

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'q_from_qid_matrix_pythorch',
    'z_from_variables_matrix_pythorch',
    'reac_bound_matrix_pytorch',
    'f_min_thrust_pytorch',
    'f_max_thrust_pytorch',
    'f_constraints_pytorch',
    'f_objective_pytorch',
    'compute_grad',
    'compute_jacobian'
]


def q_from_qid_matrix_pythorch(variables, p, A_Ik, B_Ik):
    k = B_Ik.shape[1]
    qid = variables[:k]
    q = mm(A_Ik, p) + mm(B_Ik, qid)
    return q


def z_from_variables_matrix_pythorch(variables, p, A_Ik, B_Ik, Ci, Cit, Cf, pzfree, free, fixed):
    k = B_Ik.shape[1]
    k_zb = Cf.shape[1]
    qid = variables[:k]
    zfixed = variables[-k_zb:]
    q = mm(A_Ik, p) + mm(B_Ik, qid)
    Q = diagflat(q)
    Ai = mm(mm(Cit, Q), Ci)
    Af = mm(mm(Cit, Q), Cf)
    b = pzfree - mm(Af, zfixed)
    X, LU = solve(b, Ai)
    zi = X
    z = tensor(zeros((len(free)+len(fixed), 1)))
    z[free] = zi
    z[fixed] = zfixed
    return z


def reac_bound_matrix_pytorch(variables, p, A_Ik, B_Ik, C, Cf, pfixed, xyz):
    k = B_Ik.shape[1]
    k_zb = Cf.shape[1]
    qid = variables[:k]
    zfixed = variables[-k_zb:]
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    R = mm(CfQC, xyz) - pfixed
    length_x = abs(mul(zfixed, div(R[:, 0], R[:, 2]).reshape(-1, 1)))
    length_y = abs(mul(zfixed, div(R[:, 1], R[:, 2]).reshape(-1, 1)))
    length = cat([length_x, length_y])
    return length


def f_min_thrust_pytorch(variables, p, A_Ik, B_Ik, C, Cf, xy):
    k = B_Ik.shape[1]
    qid = variables[:k]
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    Rh = mm(CfQC, xy)
    R = norm(Rh, dim=1)
    f = sum(R)
    return f


def f_max_thrust_pytorch(variables, p, A_Ik, B_Ik, C, Cf, xy):
    k = B_Ik.shape[1]
    qid = variables[:k]
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = diagflat(q)
    CfQC = mm(mm(Cf.t(), Q), C)
    Rh = mm(CfQC, xy)
    R = norm(Rh, dim=1)
    f = -1 * sum(R)
    return f


def f_objective_pytorch(variables, *args):
    A_Ik, B_Ik, p, C, Ci, Cit, Cf, pzfree, xyz, xy, pfixed, k, objective = args
    if objective == 'min':
        f = f_min_thrust_pytorch(variables, p, A_Ik, B_Ik, C, Cf, xy)
    if objective == 'max':
        f = f_max_thrust_pytorch(variables, p, A_Ik, B_Ik, C, Cf, xy)
    return f


def f_constraints_pytorch(variables, *args):
    A_Ik, B_Ik, p, C, Ci, Cit, Cf, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr = args
    if 'funicular' in dict_constr:
        constraints = q_from_qid_matrix_pythorch(variables, p, A_Ik, B_Ik)
    if 'envelope' in dict_constr:
        z = z_from_variables_matrix_pythorch(variables, p, A_Ik, B_Ik, Ci, Cit, Cf, pzfree, free, fixed)
        upper = tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - tensor(lb)
        constraints = cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_matrix_pytorch(variables, p, A_Ik, B_Ik, C, Cf, pfixed, xyz)
        reac_bound = abs(cat([tensor(b[:, 0]), tensor(b[:, 1])]).reshape(-1,1)) - length
        constraints = cat([constraints, reac_bound])
    return constraints


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
