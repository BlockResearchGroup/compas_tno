from scipy.optimize import fmin_slsqp
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
    indset = form.attributes['indset']
    find_inds = optimiser.data['find_inds']
    printout = optimiser.data['printout']

    i_k = form.index_key()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
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
        zb_bounds = [[form.vertex_attribute(i_k[i], 'lb'), form.vertex_attribute(i_k[i], 'ub')] for i in fixed]
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
    plot = optimiser.data['plot']

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
                [b_] = form.vertex_attributes(key, 'b')
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
