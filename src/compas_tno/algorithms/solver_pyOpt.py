# import pyOpt

from compas_tno.algorithms.equilibrium import reactions
from compas_tno.algorithms import zlq_from_qid
from compas.utilities import geometric_key

from numpy import array


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'run_optimisation_pyOpt'
]

def run_optimisation_pyOpt(analysis):

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
