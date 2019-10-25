import cvxpy as cp
from cvxpy import diag
from cvxpy import matrix_frac
from cvxpy import Minimize
from cvxpy import Maximize
from cvxpy import Problem
from compas_tna.diagrams import FormDiagram
from compas_thrust.algorithms import z_from_form
from numpy import zeros
from numpy import ones



__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'



__all__ = [
    'min_loadpath',
    'min_thrust',
    'min_thrust',
    'feasability'
]

def min_loadpath(form, args, printout=False):

    uv_i = form.uv_index()

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    m = form.number_of_edges()
    q = cp.Variable(m)

    # -------- OBJECTIVE FUNCTION

    fobj = matrix_frac(pz[free],Cit*cp.diag(q)*Ci) + x.T*C.T*diag(q)*Cf*x[fixed] + y.T*C.T*diag(q)*Cf*y[fixed]
    objective = Minimize(fobj)

    # -------- CONSTRAINTS

    horz = E*q == 0
    # horz = q[dep] == -Edinv*(p - Ei*q[ind])
    pos = q >= 0
    maxq = q <= 2000.0


    zmin = zeros(len(free)).reshape(len(free),1)
    zmax = 2.0*ones(len(free)).reshape(len(free),1)
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

    for u,v in form.edges():
        i = uv_i[(u,v)]
        qi= q.value[i]
        form.set_edge_attribute((u,v), 'q', qi)

    form = z_from_form(form)

    return form

def min_thrust(form, args, zmin, zmax, printout=False):

    uv_i = form.uv_index()
    
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    m = form.number_of_edges()
    q = cp.Variable(m)
    
    # -------- OBJECTIVE FUNCTION
    
    fobj = x.T*C.T*diag(q)*Cf*x[fixed] + y.T*C.T*diag(q)*Cf*y[fixed]
    objective = cp.Minimize(fobj)

    # -------- CONSTRAINTS
    
    horz = E*q == 0
    pos = q >= 0
    # Add constraints on Heights...
    zmin = zeros(len(free)).reshape(len(free),1)
    zmax = 2.0*ones(len(free)).reshape(len(free),1)
    upper = pz[free]-Cit*diag(q)*Ci*zmin >= 0
    lower = pz[free]-Cit*diag(q)*Ci*zmax <= 0
    constraints = [horz, pos, upper, lower]

    # ---------- ASSEMBLE PROBLEM

    prob = cp.Problem(objective, constraints)

    # -------- SOLVE

    lp = prob.solve(verbose=printout)
    form.attributes['solve'] = prob.solver_stats

    for u,v in form.edges():
        i = uv_i[(u,v)]
        qi = q.value[i]
        form.set_edge_attribute((u,v), 'q', qi)
    
    form = z_from_form(form)

    return form

def max_thrust(form, args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    xy = hstack([x,y])
    E   = equilibrium_matrix(C, xy, free, 'csr')
    m = form.number_of_edges()
    q = cp.Variable(m)
    
    # -------- OBJECTIVE FUNCTION
    
    fobj = x.T*C.T*diag(q)*Cf*xf + y.T*C.T*diag(q)*Cf*yf
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
    q_sol = q.value

    return q_sol

def feasability(form, args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    xy = hstack([x,y])
    E   = equilibrium_matrix(C, xy, free, 'csr')
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
    q_sol = q.value

    return q_sol
