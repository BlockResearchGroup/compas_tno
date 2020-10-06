import ipopt
from torch import tensor

from compas_tno.algorithms import f_constraints_pytorch
from compas_tno.algorithms import f_objective_pytorch
from compas_tno.algorithms import compute_autograd
from compas_tno.algorithms import compute_autograd_jacobian

from compas_tno.algorithms import reactions
from compas_tno.algorithms import zlq_from_qid

from compas.utilities import geometric_key

import time

from numpy import hstack
from numpy import array

from compas_tno.problems import sensitivities_wrapper
from compas_tno.problems import constr_wrapper
from compas_tno.problems import gradient_fmin
from compas_tno.problems import gradient_fmax
from compas_tno.problems import f_min_thrust
from compas_tno.problems import f_max_thrust

from .post_process import post_process_analysis

__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'run_optimisation_ipopt'
]


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
        variables = tensor(x.reshape(-1, 1))
        return array(self.fobj(variables, *self.args_obj))

    def gradient(self, x):
        #
        # The callback for calculating the gradient
        #
        variables = tensor(x.reshape(-1, 1), requires_grad=True)
        f = self.fobj(variables, *self.args_obj)
        return array(compute_autograd(variables, f))

    def constraints(self, x):
        #
        # The callback for calculating the constraints
        #
        variables = tensor(x.reshape(-1, 1))
        return array(self.fconstr(variables, *self.args_constr))

    def jacobian(self, x):
        #
        # The callback for calculating the Jacobian
        #
        variables = tensor(x.reshape(-1, 1), requires_grad=True)
        constraints = self.fconstr(variables, *self.args_constr)
        return array(compute_autograd_jacobian(variables, constraints)).flatten()


class wrapper_ipopt_analytical(object):
    def __init__(self):
        self.fobj = None
        self.fconstr = None
        self.args = None
        self.fgrad = None
        self.bounds = None
        self.x0 = None
        self.eps = 1e-8
        self.fjac = None
        pass

    def objective(self, x):
        #
        # The callback for calculating the objective
        #
        return self.fobj(x, *self.args)

    def gradient(self, x):
        #
        # The callback for calculating the gradient
        #
        return self.fgrad(x, *self.args)

    def constraints(self, x):
        #
        # The callback for calculating the constraints
        #
        return self.fconstr(x, *self.args).reshape(-1, 1)

    def jacobian(self, x):
        #
        # The callback for calculating the Jacobian
        #
        return self.fjac(x, *self.args).flatten()


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

    optimiser = analysis.optimiser

    constraints = optimiser.data['constraints']
    objective = optimiser.data['objective']
    printout = optimiser.data.get('printout', False)
    gradients = optimiser.data.get('gradient', False)

    bounds = optimiser.bounds
    x0 = optimiser.x0
    g0 = optimiser.g0
    args = optimiser.args

    lower = [lw[0] for lw in bounds]
    upper = [up[1] for up in bounds]

    # Tensor modification

    if not gradients:

        q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints, max_rol_rx, max_rol_ry, Asym = args[
        :48]

        EdinvEi = Edinv*Ei
        Edinv_p = Edinv.dot(p)

        EdinvEi_th = tensor(EdinvEi)
        Edinv_p_th = tensor(Edinv_p)
        C_th = tensor(C.toarray())
        Ci_th = tensor(Ci.toarray())
        Cit_th = Ci_th.t()
        Cf_th = tensor(Cf.toarray())
        pzfree = tensor(pz[free])
        xyz = tensor(hstack([x, y, z]))
        xy = tensor(hstack([x, y]))
        pfixed = tensor(hstack([px, py, pz])[fixed])
        U_th = tensor(U.toarray())
        V_th = tensor(V.toarray())

        args_obj = (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, objective)
        args_constr = (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed,
                       ub, lb, ub_ind, lb_ind, b, constraints, max_rol_rx, max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th)

        problem_obj = wrapper_ipopt()
        problem_obj.fobj = f_objective_pytorch
        problem_obj.fconstr = f_constraints_pytorch
        problem_obj.args_obj = args_obj
        problem_obj.args_constr = args_constr
        problem_obj.bounds = bounds
        problem_obj.x0 = x0
        problem_obj.args = args
        problem_obj.fgrad = gradient_fmin

        variables = tensor(x0, requires_grad=True)
        g0 = f_constraints_pytorch(variables, *args_constr)

    else:

        problem_obj = wrapper_ipopt_analytical()
        problem_obj.fobj = optimiser.fobj
        problem_obj.fconstr = optimiser.fconstr
        problem_obj.fjac = optimiser.fjac
        problem_obj.args = optimiser.args
        problem_obj.fgrad = optimiser.fgrad
        problem_obj.bounds = bounds
        problem_obj.x0 = x0

        if printout:
            g0 = optimiser.fconstr(x0, *args)

    cu = [10e10]*len(g0)
    cl = [0.0]*len(g0)
    if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in constraints):
        nsym = Asym.shape[0]
        cu[-nsym:] = [0.0]*nsym
        cl[-nsym:] = [0.0]*nsym

    nlp = ipopt.problem(
        n=len(x0),
        m=len(g0),
        problem_obj=problem_obj,
        lb=lower,
        ub=upper,
        cl=cl,
        cu=cu
    )

    # Set Options and Time
    nlp = _nlp_options(nlp, optimiser)
    start_time = time.time()

    # Solve
    xopt, info = nlp.solve(x0)
    fopt = info['obj_val']
    exitflag = info['status']
    if exitflag == 1:  # IPOPT consider solved = 1. TNO solved = 0
        exitflag = 0
    elif exitflag == 0:
        exitflag = 1
    if printout:
        print(info['status_msg'])

    elapsed_time = time.time() - start_time
    if printout:
        print('Solving Time: {0:.1f} sec'.format(elapsed_time))

    optimiser.exitflag = exitflag
    optimiser.time = elapsed_time
    optimiser.fopt = fopt
    optimiser.xopt = xopt
    optimiser.niter = None  # Did not find a way to display number of iteratinos
    optimiser.message = info['status_msg']

    post_process_analysis(analysis)

    # g_final = fconstr(xopt, *args)
    # args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i)
    # z, _, q, q_ = zlq_from_qid(q[ind], args)

    # gkeys = []
    # for i in ind:
    #     u, v = i_uv[i]
    #     gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
    # form.attributes['indset'] = gkeys

    # for i in range(form.number_of_vertices()):
    #     key = i_k[i]
    #     form.vertex_attribute(key=key, name='z', value=float(z[i]))

    # for c, qi in enumerate(list(q_.ravel())):
    #     u, v = i_uv[c]
    #     form.edge_attribute((u, v), 'q', float(qi))

    # lp = 0
    # for u, v in form.edges_where({'_is_edge': True}):
    #     if form.edge_attribute((u, v), 'is_symmetry') is False:
    #         qi = form.edge_attribute((u, v), 'q')
    #         li = form.edge_length(u, v)
    #         lp += abs(qi) * li**2
    # form.attributes['loadpath'] = float(lp)

    # # form.attributes['iter'] = niter
    # optimiser.exitflag = exitflag
    # optimiser.fopt = fopt
    # analysis.form = form
    # reactions(form, plot=plot)

    # summary = optimiser.data.get('summary', False)

    # if summary:
    #     print('\n' + '-' * 50)
    #     print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
    #     print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
    #     print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
    #     print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
    #     print('fopt      : {0:.3f}'.format(fopt))
    #     print('-' * 50 + '\n')

    return analysis


def _nlp_options(nlp, optimiser):

    # Link to instructions: https://coin-or.github.io/Ipopt/OPTIONS.html

    # nlp.addOption(b'hessian_approximation', b'limited-memory')
    # nlp.addOption('tol', 1e-4)              # Default 1e-8
    # nlp.addOption('max_iter', 500)        # Default 3000

    # nlp.addOption('dual_inf_tol', 100.0)  # Default 1.0
    # nlp.addOption('constr_viol_tol', 0.1) # Default 1e-4
    # nlp.addOption('compl_inf_tol', 0.1)   # Default 1e-4

    # nlp.addOption('acceptable_iter', 10)
    # nlp.addOption('acceptable_tol', 1e-3)
    # nlp.addOption('acceptable_constr_viol_tol', 1e-4)  # Default 1e-2
    # nlp.addOption('acceptable_dual_inf_tol', 10e10)  # Default 10e10
    # nlp.addOption('acceptable_compl_inf_tol', 1e-2)  # Default 1e-2
    # nlp.addOption('max_iter', 500)

    if not optimiser.data['printout']:
        nlp.addOption('print_level', 0)
    if optimiser.data.get('derivative_test', None):
        nlp.addOption('derivative_test', 'first-order')
        nlp.addOption('derivative_test_perturbation', 10e-8)
        # nlp.addOption('derivative_test_print_all', 'yes')

    return nlp
