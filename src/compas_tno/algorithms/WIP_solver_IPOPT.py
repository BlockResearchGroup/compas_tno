from scipy.optimize import fmin_slsqp
from scipy.optimize import shgo
from compas.numerical import devo_numpy
from compas.numerical import ga
from scipy.optimize import approx_fprime

from compas_tno.algorithms import d_fobj
from compas_tno.algorithms import d_fconstr

from ipopt import minimize_ipopt
import ipopt

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import zlq_from_qid
from compas.utilities import geometric_key

import torch as th
from torch.autograd.gradcheck import zero_gradients

from numpy import zeros
from numpy import identity
from numpy import vstack
from numpy import hstack
from numpy import multiply
from numpy import divide
from numpy import array


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'run_optimisation_scipy'
]


def q_from_qid_matrix_pythorch(variables, p, A_Ik, B_Ik):
    k = B_Ik.shape[1]
    qid = variables[:k]
    q = th.mm(A_Ik, p) + th.mm(B_Ik, qid)
    return q


def zi_from_qid_matrix_pytorch(variables, p, A_Ik, B_Ik, Ci, Cit, Cf, pzfree):
    k = B_Ik.shape[1]
    k_zb = Cf.shape[1]
    qid = variables[:k]
    zfixed = variables[-k_zb:]
    q = th.mm(A_Ik, p) + th.mm(B_Ik, qid)
    Q = th.diagflat(q)
    Ai = th.mm(th.mm(Cit, Q), Ci)
    Af = th.mm(th.mm(Cit, Q), Cf)
    b = pzfree - th.mm(Af, zfixed)
    X, LU = th.solve(b, Ai)
    zi = X
    return zi


def z_from_variables_matrix_pythorch(variables, p, A_Ik, B_Ik, Ci, Cit, Cf, pzfree, free, fixed):
    k = B_Ik.shape[1]
    k_zb = Cf.shape[1]
    qid = variables[:k]
    zfixed = variables[-k_zb:]
    q = th.mm(A_Ik, p) + th.mm(B_Ik, qid)
    Q = th.diagflat(q)
    Ai = th.mm(th.mm(Cit, Q), Ci)
    Af = th.mm(th.mm(Cit, Q), Cf)
    b = pzfree - th.mm(Af, zfixed)
    X, LU = th.solve(b, Ai)
    zi = X
    z = th.tensor(zeros((len(free)+len(fixed), 1)))
    z[free] = zi
    z[fixed] = zfixed
    return z


def reac_bound_matrix_pytorch(variables, p, A_Ik, B_Ik, C, Cf, pfixed, xyz):
    k = B_Ik.shape[1]
    k_zb = Cf.shape[1]
    qid = variables[:k]
    zfixed = variables[-k_zb:]
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = th.diagflat(q)
    CfQC = th.mm(th.mm(Cf.t(), Q), C)
    R = th.mm(CfQC, xyz) - pfixed
    length_x = abs(th.mul(zfixed, th.div(R[:, 0], R[:, 2]).reshape(-1, 1)))
    length_y = abs(th.mul(zfixed, th.div(R[:, 1], R[:, 2]).reshape(-1, 1)))
    length = th.cat([length_x, length_y])
    return length


def f_min_thrust_pytorch(variables, p, A_Ik, B_Ik, C, Cf, xy):
    k = B_Ik.shape[1]
    qid = variables[:k]
    q = q_from_qid_matrix_pythorch(qid, p, A_Ik, B_Ik)
    Q = th.diagflat(q)
    CfQC = th.mm(th.mm(Cf.t(), Q), C)
    Rh = th.mm(CfQC, xy)
    R = th.norm(Rh, dim=1)
    f = th.sum(R)
    return f


def f_objective(variables, *args):
    A_Ik, B_Ik, p, C, Ci, Cit, Cf, pzfree, xyz, xy, pfixed, k = args
    f = f_min_thrust_pytorch(variables, p, A_Ik, B_Ik, C, Cf, xy)
    return f


def f_constraints(variables, *args):
    A_Ik, B_Ik, p, C, Ci, Cit, Cf, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, dict_constr = args
    if 'funicular' in dict_constr:
        constraints = q_from_qid_matrix_pythorch(variables, p, A_Ik, B_Ik)
    if  'envelope' in dict_constr:
        z = z_from_variables_matrix_pythorch(variables, p, A_Ik, B_Ik, Ci, Cit, Cf, pzfree, free, fixed)
        upper = th.tensor(ub) - z[ub_ind]
        lower = z[lb_ind] - th.tensor(lb)
        constraints = th.cat([constraints, upper, lower])
    if 'reac_bounds' in dict_constr:
        length = reac_bound_matrix_pytorch(variables, p, A_Ik, B_Ik, C, Cf, pfixed, xyz)
        reac_bound = abs(th.cat([th.tensor(b[:, 0]), th.tensor(b[:, 1])]).reshape(-1,1)) - length
        constraints = th.cat([constraints, reac_bound])
    return constraints


def compute_grad(variables, f):
    f.backward(retain_graph=True)
    grad = variables.grad.data
    return grad


def compute_jacobian(inputs, outputs):
    d_otp = outputs.size()[0]
    d_inp = inputs.size()[0]
    jacobian = th.zeros(d_otp, d_inp)
    grad_output = th.zeros(d_otp, 1)
    for i in range(d_otp):
        zero_gradients(inputs)
        grad_output.zero_()
        grad_output[i] = 1
        outputs.backward(grad_output, retain_graph=True)
        jacobian[i] = inputs.grad.data.flatten()
    return jacobian.flatten()


def run_optimisation_ipopt(analysis):
    """ Run nonlinear optimisation problem with IPOPT

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    obj : analysis
        Analysis object optimised.

    """

    form = analysis.form
    optimiser = analysis.optimiser
    solver = optimiser.data['solver']
    fobj = optimiser.fobj
    fconstr = optimiser.fconstr
    args = optimiser.args
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints = args
    constraints = optimiser.data['constraints']
    i_uv = form.index_uv()
    i_k = form.index_key()
    k_i = form.key_index()
    bounds = optimiser.bounds
    x0 = optimiser.x0
    g0 = optimiser.g0
    plot = optimiser.data['plot']

    lower = [lw[0] for lw in bounds]
    upper = [up[1] for up in bounds]
    cu = [10e10]*len(g0)
    cl = [-10e-10]*len(g0)

    # Tensor modification

    m = len(q)
    B_Ik = zeros((m, len(ind)))
    A_Ik = zeros((m, len(p)))
    B_Ik[ind, :] = identity(len(ind))
    B_Ik[dep, :] = Edinv*Ei
    A_Ik[dep, :] = -1*Edinv.toarray()

    # Tensor Transformation

    A_Ik_th = th.tensor(A_Ik)
    B_Ik_th = th.tensor(B_Ik)
    p_th = th.tensor(p)
    C_th = th.tensor(C.toarray())
    Ci_th = th.tensor(Ci.toarray())
    Cit_th = Ci_th.t()
    Cf_th = th.tensor(Cf.toarray())
    pzfree = th.tensor(pz[free])
    xyz = th.tensor(hstack([x, y, z]))
    xy = th.tensor(hstack([x, y]))
    pfixed = th.tensor(hstack([px, py, pz])[fixed])

    args_obj = (A_Ik_th, B_Ik_th, p_th, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k)
    args_constr = (A_Ik_th, B_Ik_th, p_th, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, constraints)

    problem_obj = wrapper_ipopt()
    problem_obj.fobj = f_objective
    problem_obj.fconstr = f_constraints
    problem_obj.args_obj = args_obj
    problem_obj.args_constr = args_constr
    problem_obj.bounds = bounds
    problem_obj.x0 = x0
    print('shape variables:', x0.shape)
    print('shape inds:', len(ind))
    print('shape fixed:', len(pfixed))
    variables = th.tensor(x0, requires_grad=True)
    # qid = th.tensor(x0[:k], requires_grad=True)
    # zb = th.tensor(x0[-len(pfixed):], requires_grad=True)
    # print(qid.shape)
    print('variables tensor shape', variables.shape)
    g0 = f_constraints(variables, *args_constr)
    print('g0 shape', g0.shape)
    jac = compute_jacobian(variables, g0)
    print('jac.shape', jac.shape)
    f = f_objective(variables, *args_obj)
    print('f0: ', f)
    grad = compute_grad(variables, f)
    print('shape gradients', grad.shape)
    print('------------------ ALL WORk')

    nlp = ipopt.problem(
            n=len(x0),
            m=len(g0),
            problem_obj=problem_obj,
            lb=lower,
            ub=upper,
            cl=cl,
            cu=cu
            )

    nlp.addOption(b'hessian_approximation', b'limited-memory')
    # nlp.addOption('mu_strategy', 'adaptive')
    nlp.addOption('tol', 1e-6)
    # nlp.addOption('linear_solver', 'HSL_MA97')

    # Possible values:

    # ma27: use the Harwell routine MA27
    # ma57: use the Harwell routine MA57
    # ma77: use the Harwell routine HSL_MA77
    # ma86: use the Harwell routine HSL_MA86
    # ma97: use the Harwell routine HSL_MA97
    # pardiso: use the Pardiso package
    # wsmp: use WSMP package
    # mumps: use MUMPS package

    xopt, info = nlp.solve(x0)
    fopt = info['obj_val']
    exitflag = info['status']
    print(info['status_msg'])

    g_final = fconstr(xopt, *args)
    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i)
    z, _, q, q_ = zlq_from_qid(q[ind], args)

    gkeys = []
    for i in ind:
        u, v = i_uv[i]
        gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
    form.attributes['indset'] = gkeys

    for i in range(form.number_of_vertices()):
        key = i_k[i]
        form.vertex_attribute(key=key, name='z', value=float(z[i]))

    for c, qi in enumerate(list(q_.ravel())):
        u, v = i_uv[c]
        form.edge_attribute((u, v), 'q', float(qi))

    lp = 0
    for u, v in form.edges_where({'is_edge': True}):
        if form.edge_attribute((u, v), 'is_symmetry') is False:
            qi = form.edge_attribute((u, v), 'q')
            li = form.edge_length(u, v)
            lp += abs(qi) * li**2
    form.attributes['loadpath'] = float(lp)

    # form.attributes['iter'] = niter
    optimiser.exitflag = exitflag
    optimiser.fopt = fopt
    analysis.form = form
    reactions(form, plot=plot)

    print('\n' + '-' * 50)
    print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
    print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
    print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
    print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
    print('fopt      : {0:.3f}'.format(fopt))
    print('-' * 50 + '\n')

    return analysis


class wrapper_ipopt(object):
    def __init__(self):
        self.fobj = None
        self.fconstr = None
        self.args_obj = None
        self.args_constr = None
        self.bounds = None
        self.x0 = None
        self.eps = 1e-8
        pass
    def objective(self, x):
        #
        # The callback for calculating the objective
        #
        variables = th.tensor(x.reshape(-1,1))
        return array(self.fobj(variables, *self.args_obj))
    def gradient(self, x):
        #
        # The callback for calculating the gradient
        #
        variables = th.tensor(x.reshape(-1,1), requires_grad=True)
        f = self.fobj(variables, *self.args_obj)
        return array(compute_grad(variables, f))
    def constraints(self, x):
        #
        # The callback for calculating the constraints
        #
        variables = th.tensor(x.reshape(-1,1))
        return array(self.fconstr(variables, *self.args_constr))
    def jacobian(self, x):
        #
        # The callback for calculating the Jacobian
        #
        variables = th.tensor(x.reshape(-1,1), requires_grad=True)
        constraints = self.fconstr(variables, *self.args_constr)
        return array(compute_jacobian(variables, constraints))
    # def hessianstructure(self):
    #     #
    #     # The structure of the Hessian
    #     # Note:
    #     # The default hessian structure is of a lower triangular matrix. Therefore
    #     # this function is redundant. I include it as an example for structure
    #     # callback.
    #     #
    #     global hs

    #     hs = sps.coo_matrix(np.tril(np.ones((4, 4))))
    #     return (hs.col, hs.row)

    # def hessian(self, x, lagrange, obj_factor):
    #     #
    #     # The callback for calculating the Hessian
    #     #
    #     H = obj_factor*np.array((
    #             (2*x[3], 0, 0, 0),
    #             (x[3],   0, 0, 0),
    #             (x[3],   0, 0, 0),
    #             (2*x[0]+x[1]+x[2], x[0], x[0], 0)))

    #     H += lagrange[0]*np.array((
    #             (0, 0, 0, 0),
    #             (x[2]*x[3], 0, 0, 0),
    #             (x[1]*x[3], x[0]*x[3], 0, 0),
    #             (x[1]*x[2], x[0]*x[2], x[0]*x[1], 0)))

    #     H += lagrange[1]*2*np.eye(4)

    #     #
    #     # Note:
    #     #
    #     #
    #     return H[hs.row, hs.col]

    # def intermediate(
    #         self,
    #         alg_mod,
    #         iter_count,
    #         obj_value,
    #         inf_pr,
    #         inf_du,
    #         mu,
    #         d_norm,
    #         regularization_size,
    #         alpha_du,
    #         alpha_pr,
    #         ls_trials
    #         ):

    #     #
    #     # Example for the use of the intermediate callback.
    #     #
    #     print"Objective value at iteration #%d is - %g" % (iter_count, obj_value)
