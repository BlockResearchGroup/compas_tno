from scipy.optimize import fmin_slsqp
from compas.numerical import devo_numpy
from compas.numerical import ga

from compas_tno.algorithms.problems import initialise_problem

from compas_tno.algorithms.objectives import f_min_loadpath
from compas_tno.algorithms.objectives import f_min_loadpath_pen
from compas_tno.algorithms.objectives import f_min_thrust
from compas_tno.algorithms.objectives import f_min_thrust_pen
from compas_tno.algorithms.objectives import f_max_thrust
from compas_tno.algorithms.objectives import f_target
from compas_tno.algorithms.objectives import f_constant

from compas_tno.algorithms.constraints import f_compression
from compas_tno.algorithms.constraints import f_ub_lb
from compas_tno.algorithms.constraints import f_joints
from compas_tno.algorithms.constraints import f_cracks

from compas_tno.algorithms.equilibrium import reactions

from compas_tno.algorithms import zlq_from_qid
from compas.utilities import geometric_key

import pyOpt

from numpy.random import rand
from numpy import append
from numpy import array

from compas_tno.diagrams import FormDiagram


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'optimise_general'
]

def set_up_optimisation(analysis):

    form = analysis.form
    optimiser = analysis.optimiser

    i_k = form.index_key()

    args = initialise_problem(form) #, indset=indset, printout=printout, find_inds=find_inds, tol=tol)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    if 'reac_bounds' in optimiser.data['constraints']:
        b = set_b_constraint(form, True, True)
    else:
        b = None

    if 'cracks' in optimiser.data['constraints']:
        cracks_lb, cracks_ub = set_cracks_constraint(form, True, True)
    else:
        cracks_lb, cracks_ub = None, None

    if 'joints' in optimiser.data['constraints']:
        joints = set_joints_constraint(form, True)
    else:
        joints = None


    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind,
            ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty)

    objective = optimiser.data['objective']

    # Select Objetive

    if objective == 'loadpath':
        fobj, fconstr = f_min_loadpath, f_compression
    if objective == 'constr_lp':
        fobj, fconstr = f_min_loadpath, f_ub_lb
    if objective == 'lp_cracks':
        fobj, fconstr = f_min_loadpath, f_cracks
    if objective == 'lp_joints':
        fobj, fconstr = f_min_loadpath_pen, f_joints
    if objective == 'target':
        fobj, fconstr = f_target, f_compression
    if objective == 'min':
        fobj, fconstr = f_min_thrust, f_ub_lb
    if objective == 'min_joints':
        fobj, fconstr = f_min_thrust_pen, f_joints
    if objective == 'min_cracks':
        fobj, fconstr = f_min_thrust, f_cracks
    if objective == 'max':
        fobj, fconstr = f_max_thrust, f_ub_lb
    if objective == 'max_cracks':
        fobj, fconstr = f_max_thrust, f_cracks
    if objective == 'feasibility':
        fobj, fconstr = f_constant, f_ub_lb

    # Definition of the Variables and starting point

    variables = optimiser.data['variables']
    qmax = optimiser.data['qmax']

    if 'ind' in variables and 'zb' in variables:
        x0 = q[ind]
        zb_bounds = [[form.get_vertex_attribute(i_k[i], 'lb'), form.get_vertex_attribute(i_k[i], 'ub')] for i in fixed]
        bounds = [[-10e-6, qmax]] * k + zb_bounds
        x0 = append(x0, z[fixed]).reshape(-1, 1)
    else:
        x0 = q[ind]
        bounds = [[-10e-6, qmax]] * k

    print('Total of Independents:', len(ind))
    print('Number of Variables:', len(x0))
    f0 = fobj(x0, *args)
    g0 = fconstr(x0, *args)


    print('Non Linear Optimisation - Initial Objective Value: {0}'.format(f0))
    print('Non Linear Optimisation - Initial Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))

    optimiser.fobj = fobj
    optimiser.fconstr = fconstr
    optimiser.args = args
    optimiser.x0 = x0
    optimiser.bounds = bounds

    analysis.form = form
    analysis.optimiser = optimiser

    return analysis


def run_optimisation(analysis):

    form = analysis.form
    optimiser = analysis.optimiser
    solver = optimiser.data['solver']
    fobj = optimiser.fobj
    fconstr = optimiser.fconstr
    args = optimiser.args
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args
    i_uv = form.index_uv()
    i_k = form.index_key()
    k_i = form.key_index()
    bounds = optimiser.bounds
    x0 = optimiser.x0

    if solver == 'slsqp':
        fopt, xopt, exitflag, niter, message = _slsqp(fobj, x0, bounds, True, fconstr, args)
        while exitflag == 9:
            fopt, xopt, exitflag, niter, message = _slsqp(fobj, xopt, bounds, True, fconstr, args)
        if exitflag == 0:
            q[ind] = xopt[:k].reshape(-1, 1)
            z[fixed] = xopt[k:].reshape(-1, 1)
            # else:
            #     q[ind] = xopt[:k].reshape(-1, 1) # Code the option with only qinds
        else:
            print(message)

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
        form.set_vertex_attribute(key=key, name='z', value=float(z[i]))

    for c, qi in enumerate(list(q_.ravel())):
        u, v = i_uv[c]
        form.set_edge_attribute((u, v), 'q', float(qi))

    lp = 0
    for u, v in form.edges_where({'is_edge': True}):
        if form.get_edge_attribute((u, v), 'is_symmetry') is False:
            qi = form.get_edge_attribute((u, v), 'q')
            li = form.edge_length(u, v)
            lp += abs(qi) * li**2
    form.attributes['loadpath'] = lp

    # form.attributes['iter'] = niter
    optimiser.exitflag = exitflag
    optimiser.fopt = fopt

    reactions(form, plot=False)

    print('\n' + '-' * 50)
    print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
    print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
    print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
    print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
    print('fopt      : {0:.3f}'.format(fopt))
    print('-' * 50 + '\n')

    return analysis


def optimise_general(form, solver='slsqp', qmin=1e-6, qmax=10, find_inds=True, tol=0.001,
                     printout=10, plot=False, indset=None, tension=False, planar=False,
                     translation=None, use_bounds=None, bounds_width=5.0,
                     summary=True, objective='loadpath', bmax=False, cracks=False, rollers=False):

    # Mapping

    k_i = form.key_index()
    i_k = form.index_key()
    i_uv = form.index_uv()

    # Set-up of the problem

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds, tol=tol)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

    # Problem specific constraint setting

    b = set_b_constraint(form, bmax, printout)
    joints = set_joints_constraint(form, printout)
    cracks_lb, cracks_ub = set_cracks_constraint(form, cracks, printout)

    args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind,
            ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty)

    # Select Objetive and Constraints

    if objective == 'loadpath':
        fobj, fconstr = f_min_loadpath, f_compression
    if objective == 'constr_lp':
        fobj, fconstr = f_min_loadpath, f_ub_lb
    if objective == 'lp_cracks':
        fobj, fconstr = f_min_loadpath, f_cracks
    if objective == 'lp_joints':
        fobj, fconstr = f_min_loadpath_pen, f_joints
    if objective == 'target':
        fobj, fconstr = f_target, f_compression
    if objective == 'min':
        fobj, fconstr = f_min_thrust, f_ub_lb
    if objective == 'min_joints':
        fobj, fconstr = f_min_thrust_pen, f_joints
    if objective == 'min_cracks':
        fobj, fconstr = f_min_thrust, f_cracks
    if objective == 'max':
        fobj, fconstr = f_max_thrust, f_ub_lb
    if objective == 'max_cracks':
        fobj, fconstr = f_max_thrust, f_cracks
    if objective == 'feasibility':
        fobj, fconstr = f_constant, f_ub_lb

    # Definition of the Variables and starting point

    if translation:
        x0 = q[ind]
        zb_bounds = [[form.get_vertex_attribute(i_k[i], 'lb'), form.get_vertex_attribute(i_k[i], 'ub')] for i in fixed]
        bounds = [[qmin, qmax]] * k + zb_bounds
        x0 = append(x0, z[fixed]).reshape(-1, 1)
    else:
        x0 = q[ind]
        bounds = [[qmin, qmax]] * k

    print('Total of Independents:', len(ind))
    print('Number of Variables:', len(x0))
    f0 = fobj(x0, *args)
    g0 = fconstr(x0, *args)

    if printout:
        print('Non Linear Optimisation - Initial Objective Value: {0}'.format(f0))
        print('Non Linear Optimisation - Initial Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))
        # print('Initial step:', x0)

    # Solve with the appropriate solver

    if solver == 'slsqp':
        fopt, xopt, exitflag, niter, message = _slsqp(fobj, x0, bounds, printout, fconstr, args)
        while exitflag == 9:
            fopt, xopt, exitflag, niter, message = _slsqp(fobj, xopt, bounds, printout, fconstr, args)
        if exitflag == 0:
            if translation:
                q[ind] = xopt[:k].reshape(-1, 1)
                z[fixed] = xopt[k:].reshape(-1, 1)
            else:
                q[ind] = xopt[:k].reshape(-1, 1)
        else:
            if printout:
                print(message)

    if solver == 'devo':
        fopt, xopt = _diff_evo(fn=fobj, bounds=bounds, population=500, generations=100, printout=10, plot=False, frange=[], args=args)
        xopt = array(xopt).reshape(len(xopt), 1)
        niter = 100
        exitflag = 0
        if translation:
            q[ind] = xopt[:k]
            z[fixed] = xopt[k:].reshape(-1, 1)
        else:
            q[ind] = xopt[:k]

    if solver == 'ga':
        GA = _ga(fn=fobj, fit_type='min', num_var=len(x0), boundaries=bounds, num_pop=500, num_gen=300, args=args)
        fopt = float(GA.best_fit)
        best_i = GA.best_individual_index
        xopt = array(GA.current_pop['scaled'][best_i]).reshape(-1, 1)
        niter = 100
        exitflag = 0
        if translation:
            q[ind] = xopt[:k]
            z[fixed] = xopt[k:].reshape(-1, 1)
        else:
            q[ind] = xopt[:k]

    if solver.split('-')[0] == 'pyOpt':
        lower = [lw[0] for lw in bounds]
        upper = [up[1] for up in bounds]
        solver_pyOpt = solver.split('-')[1]
        title = 'Solver: ' + solver_pyOpt + 'Objective: ' + objective

        opt_prob = pyOpt.Optimization(title, pyOpt_wrapper)
        opt_prob.addObj('f', value=f0)
        opt_prob.addVarGroup('x', len(x0), type='c', value=x0, lower=lower, upper=upper)
        opt_prob.addConGroup('c', len(g0), type='i', value=g0)

        if solver_pyOpt == 'SLSQP':
            slv = pyOpt.SLSQP()
            slv.setOption('MAXIT', 1000)

        if solver_pyOpt == 'PSQP':
            slv = pyOpt.PSQP()

        if solver_pyOpt == 'CONMIN':
            slv = pyOpt.CONMIN()

        if solver_pyOpt == 'COBYLA':
            slv = pyOpt.COBYLA()

        if solver_pyOpt == 'SOLVOPT':
            slv = pyOpt.SOLVOPT()

        if solver_pyOpt == 'KSOPT':
            slv = pyOpt.KSOPT()

        if solver_pyOpt == 'NSGA2':
            slv = pyOpt.NSGA2()

        if solver_pyOpt == 'ALGENCAN':
            slv = pyOpt.ALGENCAN()

        if solver_pyOpt == 'FILTERSD':
            slv = pyOpt.FILTERSD()

        if solver_pyOpt == 'SDPEN':
            slv = pyOpt.SDPEN()

        if solver_pyOpt == 'ALPSO':
            slv = pyOpt.SDPEN()

        if solver_pyOpt == 'ALHSO':
            slv = pyOpt.SDPEN()

        if solver_pyOpt == 'MIDACO':
            slv = pyOpt.SDPEN()

        [fopt, xopt, info] = slv(opt_prob, sens_type='FD', args=args, objective=fobj, constraints=fconstr)
        # exitflag = info['text']
        if info['text'] == 'Optimization terminated successfully.':
            exitflag = 0
        else:
            exitflag = 1
        # exitflag = slv.informs
        # if printout:
        #     print(opt_prob.solution(0))
        print(info['text'])
        fopt = fopt.item()
        niter = 100
        if translation:
            q[ind] = xopt[:k].reshape(-1, 1)
            z[fixed] = xopt[k:].reshape(-1, 1)
        else:
            q[ind] = xopt[:k].reshape(-1, 1)

    # Recalculate equilibrium and update form-diagram

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
        form.set_vertex_attribute(key=key, name='z', value=float(z[i]))

    for c, qi in enumerate(list(q_.ravel())):
        u, v = i_uv[c]
        form.set_edge_attribute((u, v), 'q', float(qi))

    lp = 0
    for u, v in form.edges_where({'is_edge': True}):
        if form.get_edge_attribute((u, v), 'is_symmetry') is False:
            qi = form.get_edge_attribute((u, v), 'q')
            li = form.edge_length(u, v)
            lp += abs(qi) * li**2
    form.attributes['loadpath'] = lp

    # form.attributes['iter'] = niter
    form.attributes['exitflag'] = exitflag
    form.attributes['fopt'] = fopt
    form.attributes['objective'] = objective
    reactions(form, plot=plot)

    if summary:
        print('\n' + '-' * 50)
        print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
        print('zb range  : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
        print('constr    : {0:.3f} : {1:.3f}'.format(min(g_final), max(g_final)))
        print('fopt      : {0:.3f}'.format(fopt))
        print('-' * 50 + '\n')

    return fopt, q[ind], z[fixed], exitflag


def feasibility(form):

    # W.I.P.

    return


def _slsqp(fn, qid0, bounds, printout, fieq, args):

    pout = 2 if printout else 0
    opt = fmin_slsqp(fn, qid0, args=args, disp=pout, bounds=bounds, full_output=1, iter=500, f_ieqcons=fieq)

    return opt[1], opt[0], opt[3], opt[2], opt[4]


def _cobyla(fn, qid0, bounds, printout, fieq, args):

    # pout = 2 if printout else 0
    # opt  = fmin_cobyla(fn, qid0, args=args, disp=pout, bounds=bounds, full_output=1, iter=500, f_ieqcons=fieq)

    # W.I.P

    return None


def _diff_evo(fn, bounds, population, generations, printout, plot, frange, args):

    return devo_numpy(fn=fn, bounds=bounds, population=population, generations=generations, printout=printout,
                      plot=plot, frange=frange, args=args)


def _ga(fn, fit_type, num_var, boundaries, num_gen, num_pop, args):

    return ga(fit_function=fn, fit_type=fit_type, num_var=num_var, boundaries=boundaries, num_gen=num_gen, num_pop=num_pop, fargs=args)


def pyOpt_wrapper(x, **kwargs):
    args = kwargs['args']
    fobj = kwargs['objective']
    fconst = kwargs['constraints']
    # Temp
    f = fobj(x, *args)
    g = -1 * fconst(x, *args)  # pyOpt get constraints g1 <= 0 as default
    print(max(g))
    fail = 0
    return f, g, fail


def set_b_constraint(form, bmax, printout):
    if bmax:
        b = []
        for key in form.vertices_where({'is_fixed': True}):
            try:
                [b_] = form.get_vertex_attributes(key, 'b')
                b.append(b_)
            except:
                pass
        b = array(b)
    else:
        b = None
    if printout and bmax:
        print('Reaction bounds active in : {0} joints'.format(len(b)))
    return b


def set_joints_constraint(form, printout):
    try:
        joints = form.attributes['joints']
    except:
        joints = None
    if printout and joints:
        print('Constraints on the Joints set for {0} contacts.'.format(len(joints)))
    return joints


def set_cracks_constraint(form, cracks, printout):
    if cracks:
        try:
            cracks_lb, cracks_ub = form.attributes['cracks']
            if printout:
                print('Cracks Definition')
                print(cracks_lb, cracks_ub)
        except:
            cracks_lb = []
            cracks_ub = []
    else:
        cracks_lb = []
        cracks_ub = []
    print('Constraints on cracks activated in {0} lb and {1} ub.'.format(len(cracks_lb), len(cracks_ub)))
    return cracks_lb, cracks_ub
