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


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'run_optimisation_scipy'
]


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
    i_uv = form.index_uv()
    i_k = form.index_key()
    k_i = form.key_index()
    bounds = optimiser.bounds
    x0 = optimiser.x0
    g0 = optimiser.g0
    plot = optimiser.data['plot']

    lower = [lw[0] for lw in bounds]
    upper = [up[1] for up in bounds]
    cu = [10e20]*len(g0)
    cl = [0]*len(g0)

    problem_obj = wrapper_ipopt()
    problem_obj.fobj = fobj
    problem_obj.fconstr = fconstr
    problem_obj.args = args
    problem_obj.bounds = bounds
    problem_obj.x0 = x0

    nlp = ipopt.problem(
            n=len(x0),
            m=len(g0),
            problem_obj=problem_obj,
            lb=lower,
            ub=upper,
            cl=cl,
            cu=cu
            )

    nlp.addOption('mu_strategy', 'adaptive')
    nlp.addOption('tol', 1e-7)

    xopt, info = nlp.solve(x0)

    print(xopt)
    # print(info)

    # dict_constr = {'fun': fconstr}

    # exit_ = minimize_ipopt(fobj, x0, args=args, constraints=[dict_constr])

    # print(exit_)


    #     fopt, xopt, exitflag, niter, message = _slsqp(fobj, x0, bounds, True, fconstr, args)
    #     while exitflag == 9:
    #         fopt, xopt, exitflag, niter, message = _slsqp(fobj, xopt, bounds, True, fconstr, args)
    #     if exitflag == 0:
    #         q[ind] = xopt[:k].reshape(-1, 1)
    #         z[fixed] = xopt[k:].reshape(-1, 1)
    #         # else:
    #         #     q[ind] = xopt[:k].reshape(-1, 1) # Code the option with only qinds
    #     else:
    #         print(message)
    # elif solver == 'shgo':
    #     dict_constr = []
    #     for i in range(len(fconstr(x0, *args))):
    #         args_constr = list(args)
    #         args_constr.append(i)
    #         args_constr.append(fconstr)
    #         dict_ = {
    #             'type': 'ineq',
    #             'fun': _shgo_constraint_wrapper,
    #             'args': args_constr,
    #             }
    #         dict_constr.append(dict_)
    #     args_constr[len(args_constr)-2] = 0
    #     result = _shgo(fobj, bounds, True, dict_constr, args)
    #     fopt = result['fun']
    #     xopt = result['x']
    #     sucess = result['success']
    #     message = result['message']
    #     if sucess == True:
    #         exitflag = 0
    #         q[ind] = xopt[:k].reshape(-1, 1)
    #         z[fixed] = xopt[k:].reshape(-1, 1)
    #         # else:
    #         #     q[ind] = xopt[:k].reshape(-1, 1) # Code the option with only qinds
    #     else:
    #         print(message)

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
        self.args = None
        self.bounds = None
        self.x0 = None
        self.eps = 1e-8
        pass
    def objective(self, x):
        #
        # The callback for calculating the objective
        #
        return self.fobj(x.reshape(-1,1), *self.args)
    def gradient(self, x):
        #
        # The callback for calculating the gradient
        #
        return approx_fprime(x, self.fobj, self.eps, *self.args)
    def constraints(self, x):
        #
        # The callback for calculating the constraints
        #
        return self.fconstr(x.reshape(-1,1), *self.args)
    def jacobian(self, x):
        #
        # The callback for calculating the Jacobian
        #
        return d_fconstr(self.fobj, x.reshape(-1,1), self.eps, *self.args).flatten()
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
    #     print "Objective value at iteration #%d is - %g" % (iter_count, obj_value)
