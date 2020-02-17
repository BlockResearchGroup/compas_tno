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

from compas_tna.diagrams import FormDiagram


__author__ = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__ = 'MIT License'
__email__ = 'mricardo@ethz.ch'


__all__ = [
    'MMA_solver'
]


def MMA_solver(form, solver='slsqp', qmin=1e-6, qmax=10, find_inds=True, tol=0.001,
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

    # Solve with MMA
    # Start Time
    start_time = time.time()

    # Logger
    path = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(path, "GCMMA_TEST.log")
    logger = setup_logger(file)
    logger.info("Started\n")
    # Set numpy print options
    np.set_printoptions(precision=4, formatter={'float': '{: 0.4f}'.format})
    # Beam initial settings
    m = len(g0)
    n = len(x0)
    epsimin = 0.0000001
    # epsimin = 0.1
    eeen = np.ones((n, 1))
    eeem = np.ones((m, 1))
    zeron = np.zeros((n, 1))
    zerom = np.zeros((m, 1))
    xval = x0
    xold1 = xval.copy()
    xold2 = xval.copy()
    xmin = lower  # eeen.copy()
    xmax = upper  # 10*eeen
    low = xmin.copy()
    upp = xmax.copy()
    c = 1000*eeem
    d = eeem.copy()
    a0 = 1
    a = zerom.copy()
    raa0 = 0.01
    raa = 0.01*eeem
    raa0eps = 0.000001
    raaeps = 0.000001*eeem
    outeriter = 0
    maxoutit = 50
    kkttol = 1e-6
    # Calculate function values and gradients of the objective and constraints functions
    if outeriter == 0:
        f0val, df0dx, fval, dfdx = beam2(xval, fobj, fconstr, *args)
        innerit = 0
        outvector1 = np.array([outeriter, innerit, f0val, fval])
        outvector2 = xval.flatten()
        # Log
        # logger.info("outvector1 = {}".format(outvector1))
        logger.info("outvector2 = {}\n".format(outvector2))
    # The iterations starts
    kktnorm = kkttol + 10
    f0val_old = f0val + 10
    outit = 0
    while (kktnorm > kkttol) and (outit < maxoutit): # and (abs(f0val_old - f0val) > 10e-4):
        logger.info("iteration = {}".format(outit))
        f0val_old = f0val
        outit += 1
        outeriter += 1
        # The parameters low, upp, raa0 and raa are calculated:
        low, upp, raa0, raa = \
            asymp(outeriter, n, xval, xold1, xold2, xmin, xmax,
                  low, upp, raa0, raa, raa0eps, raaeps, df0dx, dfdx)
        # The MMA subproblem is solved at the point xval:
        xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, f0app, fapp = \
            gcmmasub(m, n, iter, epsimin, xval, xmin, xmax, low, upp,
                     raa0, raa, f0val, df0dx, fval, dfdx, a0, a, c, d)
        # The user should now calculate function values (no gradients) of the objective- and constraint
        # functions at the point xmma ( = the optimal solution of the subproblem).
        f0valnew, fvalnew = beam1(xmma, fobj, fconstr, *args)
        # It is checked if the approximations are conservative:
        conserv = concheck(m, epsimin, f0app, f0valnew, fapp, fvalnew)
        # While the approximations are non-conservative (conserv=0), repeated inner iterations are made:
        innerit = 0
        if conserv == 0:
            while conserv == 0 and innerit <= 15:
                innerit += 1
                # New values on the parameters raa0 and raa are calculated:
                raa0, raa = raaupdate(xmma, xval, xmin, xmax, low, upp, f0valnew, fvalnew, f0app, fapp, raa0,
                                      raa, raa0eps, raaeps, epsimin)
                # The GCMMA subproblem is solved with these new raa0 and raa:
                xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, f0app, fapp = gcmmasub(m, n, iter, epsimin, xval, xmin,
                                                                                    xmax, low, upp, raa0, raa, f0val, df0dx, fval, dfdx, a0, a, c, d)
                # The user should now calculate function values (no gradients) of the objective- and
                # constraint functions at the point xmma ( = the optimal solution of the subproblem).
                f0valnew, fvalnew = beam1(xmma, fobj, fconstr, *args)
                # It is checked if the approximations have become conservative:
                conserv = concheck(m, epsimin, f0app, f0valnew, fapp, fvalnew)
        # Some vectors are updated:
        xold2 = xold1.copy()
        xold1 = xval.copy()
        xval = xmma.copy()
        # Re-calculate function values and gradients of the objective and constraints functions
        f0val, df0dx, fval, dfdx = beam2(xval, fobj, fconstr, *args)
        # The residual vector of the KKT conditions is calculated
        residu, kktnorm, residumax = \
            kktcheck(m, n, xmma, ymma, zmma, lam, xsi, eta, mu, zet,
                     s, xmin, xmax, df0dx, fval, dfdx, a0, a, c, d)
        outvector1 = np.array([outeriter, innerit, f0val, fval])
        outvector2 = xval.flatten()
        # Log
        # logger.info("outvector1 = {}".format(outvector1))
        # logger.info("outvector2 = {}".format(outvector2))
        logger.info("kktnorm    = {}\n".format(kktnorm))
        logger.info("Objective Iteration: {}".format(f0val))
    # Final log
    elapsed_time = time.time() - start_time
    logger.info("Finished")
    logger.info("Objective Value: {}".format(f0val))
    print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    fopt = f0val.item()
    xopt = xval

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
    # form.attributes['exitflag'] = exitflag
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
        print('Reaction spread in : {0} joints'.format(len(b)))
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
