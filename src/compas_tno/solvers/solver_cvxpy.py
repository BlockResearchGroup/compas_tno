from compas_tno.algorithms import xyz_from_q

from numpy import array

from compas_tno.algorithms import compute_reactions

from compas_tno.problems import initialise_problem_general
from compas_tno.problems import adapt_problem_to_fixed_diagram

import cvxpy as cp
from cvxpy import diag
from cvxpy import matrix_frac
from cvxpy import Minimize
from cvxpy import Problem

# from numpy import hstack

# from compas_tno.algorithms import equilibrium_fdm
# from compas_tno.problems import initialise_problem

# from compas.numerical import equilibrium_matrix

# from numpy import zeros
# from numpy import ones
# from numpy import array

# import matlab.engine

# from compas.utilities import geometric_key
# from compas_tno.algorithms.equilibrium import compute_reactions


__all__ = [
    'optimise_convex',
    'call_cvx',
    'call_cvx_ind',
    'min_loadpath',
    'min_thrust',
    'min_thrust',
    'feasibility'
]


def run_loadpath_from_form_CVXPY(form, problem=None, find_inds=False, printout=False):
    """Run convex optimisation problem with CVXPY directly from the Form Diagram
        OBS: Requires installation of CVXPY and ....

    Parameters
    ----------
    form : FormDiagram
        FormDiagram object with form containing full information.
    problem : Problem, optional
        The problem with matrices of interest, by default None
    find_inds : bool, optional
        Whether or not independents must be computed before the analysis, by default False
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    obj : analysis
        Analysis object optimised.
    """

    if not problem:
        problem = initialise_problem_general(form)

    if find_inds:
        adapt_problem_to_fixed_diagram(problem, form)

    output = call_and_output_CVXPY(form, problem, printout=printout)

    return output


def call_and_output_CVXPY(form, problem, printout=False):
    """Call and output the loadpath optimisation with CVXPY

    Parameters
    ----------
    form : ::class:: FormDiagram
        The form Diagram of the analysis
    problem : ::class:: Problem
        The Problem with relevant matrices and vectors
    eng : matlab.engine
        The matlab engine initiated to call the analysis
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    output: dict
        Dictionary with results.
    """

    # if len(problem.ind) < problem.m:
    #     print('Calling LP-Optimisation via CVXPY with independents')
    #     fopt, qopt, exitflag, niter, status, sol_time = call_cvx_ind(problem, eng, printout=printout)
    # else:
    print('Calling LP-Optimisation via CVXPY with NO independents')
    # fopt, qopt, exitflag, niter, status, sol_time = call_cvxpy(problem, printout=printout)
    fopt = call_cvxpy(problem, printout=printout)

    problem.q = qopt
    Xfinal = xyz_from_q(problem.q, problem.P[problem.free], problem.X[problem.fixed], problem.Ci, problem.Cit, problem.Cb)
    problem.X[problem.free, 2] = Xfinal[:, 2]

    # Update form diagram optimised

    i_uv = problem.i_uv

    i = 0
    for key in form.vertices():
        # form.vertex_attribute(key, 'x', problem.X[i, 0])
        # form.vertex_attribute(key, 'y', problem.X[i, 1])
        form.vertex_attribute(key, 'z', problem.X[i, 2])
        i = i + 1

    for c, qi in enumerate(list(problem.q.ravel())):
        u, v = i_uv[c]
        li = form.edge_length(u, v)
        form.edge_attribute((u, v), 'q', float(qi))
        form.edge_attribute((u, v), 'f', float(qi*li))

    form.attributes['loadpath'] = form.loadpath()
    compute_reactions(form)

    summary = True

    # Output dictionary

    output = {}
    output['fopt'] = fopt
    output['exitflag'] = exitflag
    output['status'] = status
    output['niter'] = niter
    output['sol_time'] = sol_time

    if printout or summary:
        print('\n' + '-' * 50)
        print('LOADPATH OPTIMISATION WITH CVX (MATLAB)')
        print('status    :', status)
        print('fopt (lp) : {0:.3f}'.format(fopt))
        print('n-iter    : {0}'.format(niter))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(qopt)[0], max(qopt)[0]))
        print('sol. time : {0:.3f} sec'.format(sol_time))
        print('-' * 50 + '\n')

    return output


def call_cvxpy(problem, printout=False):

    q = problem.q
    E = problem.E
    C = problem.C
    Ci = problem.Ci
    Cit = problem.Cit
    Cf = problem.Cb
    pz = problem.P[:, 2].reshape(-1, 1)
    free = problem.free
    fixed = problem.fixed
    x = problem.x0
    y = problem.y0
    m = problem.m

    q = cp.Variable(m)

    fobj = matrix_frac(pz[free], Cit*cp.diag(q)*Ci) + x.T*C.T*diag(q)*Cf*x[fixed] + y.T*C.T*diag(q)*Cf*y[fixed]
    objective = Minimize(fobj)

    horz = E*q == 0
    pos = q >= 0
    maxq = q <= 2000.0

    # zmin = zeros(len(free)).reshape(len(free), 1)
    # zmax = 2.0*ones(len(free)).reshape(len(free), 1)
    # lower = pz[free]-Cit*diag(q)*Ci*zmin >= 0
    # upper = pz[free]-Cit*diag(q)*Ci*zmax <= 0
    # constraints = [horz, pos, maxq, upper, lower]
    constraints = [horz, pos, maxq]

    prob = Problem(objective, constraints)
    fopt = prob.solve(verbose=True)

    return fopt


# -----------------------
# -----------------------
# ALL OF THE BELOW IMPLEMENT THE CONVEX OPTIMISATION WITH CVX_PY -> THERE'S NO SOLVER SDPT3 ON THAT -> VERIFY IF SHOULD DELETE
# -----------------------
# -----------------------


def min_loadpath(form, args, printout=False):

    uv_i = form.uv_index()
    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target,
     s, Wfree, anchors, x, y, b) = args
    m = form.number_of_edges()
    q = cp.Variable(m)

    # -------- OBJECTIVE FUNCTION

    fobj = matrix_frac(pz[free], Cit*cp.diag(q)*Ci) + x.T*C.T*diag(q)*Cf*x[fixed] + y.T*C.T*diag(q)*Cf*y[fixed]
    objective = Minimize(fobj)

    # -------- CONSTRAINTS

    horz = E*q == 0
    # horz = q[dep] == -Edinv*(p - Ei*q[ind])
    pos = q >= 0
    maxq = q <= 2000.0

    zmin = zeros(len(free)).reshape(len(free), 1)
    zmax = 2.0*ones(len(free)).reshape(len(free), 1)
    lower = pz[free]-Cit*diag(q)*Ci*zmin >= 0
    upper = pz[free]-Cit*diag(q)*Ci*zmax <= 0
    constraints = [horz, pos, maxq, upper, lower]

    # ---------- ASSEMBLE PROBLEM

    prob = Problem(objective, constraints)

    # import sdpt3glue
    # matfile_target = '/Users/mricardo/Documents/MATLAB/optimisation/test_SDPT3_in.mat'
    # output_target = '/Users/mricardo/Documents/MATLAB/optimisation/test_SDPT3_out.mat'

    # sdpt3glue.sdpt3_solve_problem(prob, sdpt3glue.MATLAB, matfile_target,
    #                                    output_target=output_target)

    # -------- SOLVE

    form.attributes['loadpath'] = prob.solve(verbose=printout)
    form.attributes['solve'] = prob.solver_stats

    for u, v in form.edges():
        i = uv_i[(u, v)]
        qi = q.value[i]
        form.edge_attribute((u, v), 'q', qi)

    form = equilibrium_fdm(form)

    return form


def min_thrust(form, args, zmin, zmax, printout=False):

    uv_i = form.uv_index()

    (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target,
     s, Wfree, anchors, x, y, b) = args
    m = form.number_of_edges()
    q = cp.Variable(m)

    # -------- OBJECTIVE FUNCTION

    fobj = x.T*C.T*diag(q)*Cf*x[fixed] + y.T*C.T*diag(q)*Cf*y[fixed]
    objective = cp.Minimize(fobj)

    # -------- CONSTRAINTS

    horz = E*q == 0
    pos = q >= 0
    # Add constraints on Heights...
    zmin = zeros(len(free)).reshape(len(free), 1)
    zmax = 2.0*ones(len(free)).reshape(len(free), 1)
    upper = pz[free]-Cit*diag(q)*Ci*zmin >= 0
    lower = pz[free]-Cit*diag(q)*Ci*zmax <= 0
    constraints = [horz, pos, upper, lower]

    # ---------- ASSEMBLE PROBLEM

    prob = cp.Problem(objective, constraints)

    # -------- SOLVE

    lp = prob.solve(verbose=printout)
    form.attributes['solve'] = prob.solver_stats
    form.attributes['loadpath'] = lp

    for u, v in form.edges():
        i = uv_i[(u, v)]
        qi = q.value[i]
        form.edge_attribute((u, v), 'q', qi)

    form = equilibrium_fdm(form)

    return form


def max_thrust(form, args):

    (q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target,
     s, Wfree, anchors, x, y, b) = args
    xy = hstack([x, y])
    E = equilibrium_matrix(C, xy, free, 'csr')
    m = form.number_of_edges()
    q = cp.Variable(m)

    # -------- OBJECTIVE FUNCTION

    fobj = x.T*C.T*diag(q)*Cf*x[fixed] + y.T*C.T*diag(q)*Cf*y[fixed]
    objective = cp.Minimize(fobj)

    # -------- CONSTRAINTS

    horz = E*q == 0
    pos = q >= 0
    # Add constraints on Heights...
    constraints = [horz, pos]

    # ---------- ASSEMBLE PROBLEM

    prob = cp.Problem(objective, constraints)

    # -------- SOLVE

    lp = prob.solve()
    form.attributes['loadpath'] = lp

    return q.value


def feasibility(form, args):

    (q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target,
     s, Wfree, anchors, x, y, b) = args
    xy = hstack([x, y])
    E = equilibrium_matrix(C, xy, free, 'csr')
    m = form.number_of_edges()
    q = cp.Variable(m)

    # -------- OBJECTIVE FUNCTION

    fobj = 1.0
    objective = cp.Minimize(fobj)

    # -------- CONSTRAINTS

    horz = E*q == 0
    pos = q >= 0
    constraints = [horz, pos]

    # ---------- ASSEMBLE PROBLEM

    prob = cp.Problem(objective, constraints)

    # -------- SOLVE

    lp = prob.solve()
    form.attributes['loadpath'] = lp

    return q.value
