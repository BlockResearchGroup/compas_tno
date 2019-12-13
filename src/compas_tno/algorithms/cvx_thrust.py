import cvxpy as cp
from cvxpy import diag
from cvxpy import matrix_frac
from cvxpy import Minimize
from cvxpy import Maximize
from cvxpy import Problem
from compas_tna.diagrams import FormDiagram
from compas_tno.algorithms import z_from_form
from compas_tno.algorithms.problems import initialise_problem
from compas_tno.algorithms import zlq_from_qid
from compas_tno.algorithms import q_from_qid
from compas_tno.algorithms import zlq_from_q
from numpy import zeros
from numpy import ones
import matlab.engine

from numpy import array

from compas.utilities import geometric_key
from compas_tno.algorithms.equilibrium import reactions


__all__ = [
    'optimise_convex',
    'call_cvx',
    'call_cvx_ind',
    'min_loadpath',
    'min_thrust',
    'min_thrust',
    'feasibility'
]


def optimise_convex(form, qmin=1e-6, qmax=10, find_inds=True, tol=0.001,
                    printout=1, plot=False, indset=None, tension=False, planar=False,
                    translation=None, summary=True, objective='loadpath'):

    # Mapping

    k_i = form.key_index()
    i_k = form.index_key()
    i_uv = form.index_uv()
    uv_i = form.uv_index()

    # Set-up of the problem and start matlab engine

    future = matlab.engine.start_matlab(background=True)
    eng = future.result()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds, tol=tol)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cfx, Cfy = args

    # Problem specifities

    args_cvx = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, eng, qmax, i_uv, k_i)

    # Select Objetive and Constraints

    if find_inds or indset:
        print('Calling CVX WITH independents')
        fopt, qopt, exitflag, niter = call_cvx_ind(objective, args_cvx)
        z, _, q, q_ = zlq_from_qid(qopt[ind], args)
    else:
        print('Calling CVX with NO independents')
        fopt, qopt, exitflag, niter = call_cvx(objective, args_cvx)
        z, _, q, q_ = zlq_from_q(qopt, args)

    # Recalculate equilibrium and update form-diagram

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
    for u, v in form.edges():
        if form.get_edge_attribute((u, v), 'is_symmetry') is False:
            qi = form.get_edge_attribute((u, v), 'q')
            li = form.edge_length(u, v)
            lp += abs(qi) * li**2
    form.attributes['loadpath'] = lp

    form.attributes['iter'] = niter
    form.attributes['exitflag'] = exitflag
    reactions(form, plot=plot)

    if summary:
        print('\n' + '-' * 50)
        print('qid range : {0:.3f} : {1:.3f}'.format(min(q[ind])[0], max(q[ind])[0]))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(q)[0], max(q)[0]))
        print('zb range : {0:.3f} : {1:.3f}'.format(min(z[fixed])[0], max(z[fixed])[0]))
        print('fopt      : {0:.3f}'.format(fopt))
        print('-' * 50 + '\n')

    return fopt, q[ind], z[fixed], exitflag


def call_cvx(objective, args_cvx):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, eng, qmax, i_uv, k_i = args_cvx

    eng.workspace['m'] = int(len(q))
    eng.workspace['C'] = matlab.double(C.toarray().tolist())
    eng.workspace['Ci'] = matlab.double(Ci.toarray().tolist())
    eng.workspace['Cf'] = matlab.double(Cf.toarray().tolist())
    eng.workspace['Cit'] = matlab.double(Cit.toarray().tolist())
    eng.workspace['xt'] = matlab.double(x.transpose().tolist())
    eng.workspace['yt'] = matlab.double(y.transpose().tolist())
    eng.workspace['xb'] = matlab.double(x[fixed].tolist())
    eng.workspace['yb'] = matlab.double(y[fixed].tolist())
    eng.workspace['pz'] = matlab.double(pz[free].tolist())
    eng.workspace['qmax'] = float(qmax)
    eng.workspace['E'] = matlab.double(E.tolist())

    eng.save('/Users/mricardo/Documents/MATLAB/optimisation/discretize/form2.mat', nargout=0)

    print('Solving with SDPT3')
    eng.cvx_begin(nargout=0)
    eng.variable('q(double(m))', nargout=0)
    if objective == 'loadpath':
        eng.minimize('matrix_frac(pz,(Cit*diag(q)*Ci)) + xt*transpose(C)*diag(q)*Cf*xb + yt*transpose(C)*diag(q)*Cf*yb', nargout=0)
    if objective == 'feasibility':
        eng.minimize('1', nargout=0)
    eng.eval('q >= -0.001', nargout=0)
    eng.eval('q <= qmax', nargout=0)
    eng.eval('E * q == 0.0', nargout=0)
    eng.cvx_end(nargout=0)

    fopt = eng.workspace['cvx_optval']
    qopt = array(eng.workspace['q'])
    status = eng.workspace['cvx_status']
    niter = eng.workspace['cvx_slvitr']

    if status is not 'Infeasible':
        exitflag = 0
    else:
        exitflag = 1

    return fopt, qopt, exitflag, niter


def call_cvx_ind(fobj, args_cvx):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, eng, qmax, i_uv, k_i = args_cvx
    ind_ = [x+1 for x in ind]
    dep_ = [x+1 for x in dep]

    eng.workspace['m'] = int(len(q))
    eng.workspace['C'] = matlab.double(C.toarray().tolist())
    eng.workspace['Ci'] = matlab.double(Ci.toarray().tolist())
    eng.workspace['Cf'] = matlab.double(Cf.toarray().tolist())
    eng.workspace['dep'] = matlab.double(dep_)
    eng.workspace['ind'] = matlab.double(ind_)
    eng.workspace['Edinv'] = matlab.double(Edinv.toarray().tolist())
    eng.workspace['Ei'] = matlab.double(Ei.tolist())
    eng.workspace['p'] = matlab.double(p.tolist())
    eng.workspace['Cit'] = matlab.double(Cit.toarray().tolist())
    eng.workspace['xt'] = matlab.double(x.transpose().tolist())
    eng.workspace['yt'] = matlab.double(y.transpose().tolist())
    eng.workspace['xb'] = matlab.double(x[fixed].tolist())
    eng.workspace['yb'] = matlab.double(y[fixed].tolist())
    eng.workspace['pz'] = matlab.double(pz[free].tolist())
    eng.workspace['qmax'] = float(qmax)
    eng.workspace['E'] = matlab.double(E.tolist())

    eng.cvx_begin(nargout=0)
    eng.variable('q(double(m))', nargout=0)
    eng.minimize('matrix_frac(pz,(Cit*diag(q)*Ci)) + xt*transpose(C)*diag(q)*Cf*xb + yt*transpose(C)*diag(q)*Cf*yb', nargout=0)
    eng.eval('q >= 0.0', nargout=0)
    eng.eval('q <= qmax', nargout=0)
    eng.eval('q(dep) == - Edinv*(p - Ei*q(ind))', nargout=0)
    eng.cvx_end(nargout=0)

    fopt = eng.workspace['cvx_optval']
    status = eng.workspace['cvx_status']
    qopt = array(eng.workspace['q'])
    niter = eng.workspace['cvx_slvitr']

    if status is not 'Infeasible':
        exitflag = 0
    else:
        exitflag = 1

    return fopt, qopt, exitflag, niter


# -----------------------
# -----------------------
# ALL OF THE BELOW IMPLEMENT THE CONVEX OPTIMISATION WITH CVX_PY -> THERE'S NO SOLVER SDPT3 ON THAT -> VERIFY IF SHOULD DELETE
# -----------------------
# -----------------------


def min_loadpath(form, args, printout=False):

    uv_i = form.uv_index()

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
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
        form.set_edge_attribute((u, v), 'q', qi)

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

    for u, v in form.edges():
        i = uv_i[(u, v)]
        qi = q.value[i]
        form.set_edge_attribute((u, v), 'q', qi)

    form = z_from_form(form)

    return form


def max_thrust(form, args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
    xy = hstack([x, y])
    E = equilibrium_matrix(C, xy, free, 'csr')
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


def feasibility(form, args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args
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
    q_sol = q.value

    return q_sol
