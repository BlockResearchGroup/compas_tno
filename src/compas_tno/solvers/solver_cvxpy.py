from compas_tno.algorithms import xyz_from_q
from compas_tno.algorithms import compute_reactions

from compas_tno.problems import initialise_problem_general
from compas_tno.problems import adapt_problem_to_fixed_diagram

import cvxpy as cp
from cvxpy import diag
from cvxpy import matrix_frac
from cvxpy import Minimize
from cvxpy import Problem


def run_optimisation_CVXPY(analysis):
    """ Run convex optimisation problem with CVXPY after going through the optimisation set up.

    Parameters
    ----------
    analysis : Analysis
        Analysis object with information about optimiser, form and shape.

    """

    form = analysis.form
    problem = analysis.optimiser.M
    find_inds = analysis.optimiser.settings.get('find_inds', False)
    printout = analysis.optimiser.settings.get('printout', False)

    problem = run_loadpath_from_form_CVXPY(form, problem=problem, find_inds=find_inds, printout=printout)

    return problem


def run_loadpath_from_form_CVXPY(form, problem=None, find_inds=False, printout=False):
    """ Run convex optimisation problem with CVXPY directly from the Form Diagram
    OBS: Requires installation of CVXPY and MOSEK

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

    problem = call_and_output_CVXPY(form, problem, printout=printout)

    return problem


def call_and_output_CVXPY(form, problem, printout=False):
    """Call and output the loadpath optimisation with CVXPY

    Parameters
    ----------
    form : ::class:: FormDiagram
        The form Diagram of the analysis
    problem : ::class:: Problem
        The Problem with relevant matrices and vectors`
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    output: dict
        Dictionary with results.
    """

    if len(problem.ind) < problem.m:
        print('Calling LP-Optimisation via CVXPY with independents')
        fopt, qopt, exitflag, niter, status, sol_time = call_cvxpy_ind(problem, printout=printout)
    else:
        print('Calling LP-Optimisation via CVXPY with NO independents')
        fopt, qopt, exitflag, niter, status, sol_time = call_cvxpy(problem, printout=printout)

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
        print('LOADPATH OPTIMISATION WITH CVXPY')
        print('status    :', status)
        print('fopt (lp) : {0:.3f}'.format(fopt))
        print('n-iter    : {0}'.format(niter))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(qopt), max(qopt)))
        print('sol. time : {0:.3f} sec'.format(sol_time))
        print('-' * 50 + '\n')

    return problem


def call_cvxpy(problem, printout=False):
    """Call and output the loadpath optimisation with CVXPY

    Parameters
    ----------
    problem : ::class:: Problem
        The Problem with relevant matrices and vectors`
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    fopt : float
        Objective function value. Loadpath.
    qopt : array
        Independent Force densities in the optimum.
    exitflag : int
        Whether or not optimisation worked.
    niter : int
        Number of iterations.
    status : str
        Message with statuss.
    sol_time : dict
        Time to solve optimisation.
    """

    # Retrieve important vector and matrices
    q = problem.q
    E = problem.E
    C = problem.C
    Ci = problem.Ci
    Cit = problem.Cit
    Cb = problem.Cb
    pz = problem.P[:, 2].reshape(-1, 1)
    ph = problem.ph
    free = problem.free
    fixed = problem.fixed
    qmin = problem.qmin
    qmax = problem.qmax
    x = problem.x0
    y = problem.y0
    m = problem.m

    # # If need to store the matrices/vectors in
    # dict_parameters = {}
    # dict_parameters['q'] = q.tolist()
    # dict_parameters['qmin'] = qmin.tolist()
    # dict_parameters['qmax'] = qmax.tolist()
    # dict_parameters['E'] = E.tolist()
    # dict_parameters['C'] = C.todense().tolist()
    # dict_parameters['Ci'] = Ci.todense().tolist()
    # dict_parameters['Cit'] = Cit.todense().tolist()
    # dict_parameters['Cb'] = Cb.todense().tolist()
    # dict_parameters['pz'] = pz.tolist()
    # dict_parameters['free'] = free
    # dict_parameters['fixed'] = fixed
    # dict_parameters['x'] = x.tolist()
    # dict_parameters['y'] = y.tolist()
    # dict_parameters['m'] = m

    # import json

    # with open('/Users/mricardo/compas_dev/compas_tno/data/data.json', 'w') as outfile:
    #     json.dump(dict_parameters, outfile)

    q = cp.Variable(m)
    # fobj = matrix_frac(pz[free], Cit@cp.diag(q)@Ci) + x.T@C.T@diag(q)@Cb@x[fixed] + y.T@C.T@diag(q)@Cb@y[fixed]  # for q positive
    fobj = matrix_frac(pz[free], - Cit@cp.diag(q)@Ci) - x.T@C.T@diag(q)@Cb@x[fixed] - y.T@C.T@diag(q)@Cb@y[fixed]  # for q negative
    objective = Minimize(fobj)

    horz = E@q == ph.flatten()
    pos = q >= qmin.flatten()
    maxq = q <= qmax.flatten()

    constraints = [horz, pos, maxq]

    prob = Problem(objective, constraints)
    prob.solve(solver='MOSEK', verbose=printout)
    # prob.solve(solver='MOSEK', verbose=True)

    # save output
    fopt = prob.value
    qopt = q.value
    status = prob.status
    niter = prob.solver_stats.num_iters
    sol_time = prob.solver_stats.solve_time

    if status not in ["infeasible", "unbounded"]:
        exitflag = 1
    else:
        exitflag = 0

    return fopt, qopt, exitflag, niter, status, sol_time


def call_cvxpy_ind(problem, printout=False):
    """Call and output the loadpath optimisation with CVXPY using independents

    Parameters
    ----------
    problem : ::class:: Problem
        The Problem with relevant matrices and vectors`
    printout : bool, optional
        Whether or not print results, by default False

    Returns
    -------
    fopt : float
        Objective function value. Loadpath.
    qopt : array
        Independent Force densities in the optimum.
    exitflag : int
        Whether or not optimisation worked.
    niter : int
        Number of iterations.
    status : str
        Message with statuss.
    sol_time : dict
        Time to solve optimisation.
    """

    # Retrieve important vector and matrices
    q = problem.q
    Edinv = problem.Edinv
    Ei = problem.Ei
    ph = problem.ph
    C = problem.C
    Ci = problem.Ci
    Cit = problem.Cit
    Cb = problem.Cb
    pz = problem.P[:, 2].reshape(-1, 1)
    free = problem.free
    fixed = problem.fixed
    dep = problem.dep
    ind = problem.ind
    qmin = problem.qmin
    qmax = problem.qmax
    x = problem.x0
    y = problem.y0
    m = problem.m

    q = cp.Variable(m)

    fobj = matrix_frac(pz[free], -Cit@cp.diag(q)@Ci) - x.T@C.T@diag(q)@Cb@x[fixed] - y.T@C.T@diag(q)@Cb@y[fixed]
    objective = Minimize(fobj)

    horz = q[dep] == Edinv@(Ei@q[ind] - ph.flatten())
    pos = q >= qmin.flatten()
    maxq = q <= qmax.flatten()

    constraints = [horz, pos, maxq]

    prob = Problem(objective, constraints)
    prob.solve(solver='MOSEK', verbose=printout)

    # save output
    fopt = prob.value
    qopt = q.value
    status = prob.status
    niter = prob.solver_stats.num_iters
    sol_time = prob.solver_stats.solve_time

    if status not in ["infeasible", "unbounded"]:
        exitflag = 1
    else:
        exitflag = 0

    return fopt, qopt, exitflag, niter, status, sol_time
