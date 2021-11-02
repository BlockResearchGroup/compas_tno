from compas_tno.algorithms import xyz_from_q

import matlab.engine
import time

from numpy import array

from compas_tno.algorithms import reactions

from compas_tno.problems import initialise_problem_general
from compas_tno.problems import adapt_problem_to_fixed_diagram

__all__ = [
    'run_optimisation_MATLAB',
    'run_loadpath_from_form_MATLAB'
]


def run_optimisation_MATLAB(analysis):
    """ Run convex optimisation problem with MATLAB.

    Parameters
    ----------
    obj : analysis
        Analysis object with information about optimiser, form and shape.

    Returns
    -------
    obj : analysis
        Analysis object optimised.

    """

    # Initiate Matlab Engine

    future = matlab.engine.start_matlab(background=True)
    eng = future.result()

    form = analysis.form
    optimiser = analysis.optimiser
    args_cvx = optimiser.args
    indset = form.attributes['indset']
    find_inds = optimiser.settings['find_inds']
    objective = optimiser.settings['objective']
    printout = optimiser.settings['printout']
    plot = optimiser.settings['plot']

    i_k = form.index_key()
    i_uv = form.index_uv()

    q, ind = args_cvx[:2]
    args_cvx_ = list(args_cvx)
    args_cvx_.append(eng)
    args_cvx = tuple(args_cvx_)
    args = args_cvx[:22]

    fopt, exitflag = call_and_output_CVX_MATLAB(form, find_inds, indset, printout, plot, objective, args_cvx, args, ind, i_uv, i_k)

    optimiser.exitflag = exitflag
    optimiser.fopt = fopt
    analysis.form = form
    reactions(form, plot=plot)

    return analysis


def run_loadpath_from_form_MATLAB(form, problem=None, find_inds=False, printout=False):
    """ Run convex optimisation problem with MATLAB directly from the Form Diagram
        OBS: Requires installation of CVX and MATLAB.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram object with form containing full information.

    Returns
    -------
    obj : analysis
        Analysis object optimised.

    """

    future = matlab.engine.start_matlab(background=True)

    eng = future.result()

    if not problem:
        problem = initialise_problem_general(form)

    if find_inds:
        adapt_problem_to_fixed_diagram(problem, form)

    output = call_and_output_CVX_MATLAB(form, problem, eng)

    return output


def call_and_output_CVX_MATLAB(form, problem, eng, printout=False):

    if len(problem.ind) < problem.m:
        print('Calling LP-Optimisation via CVX (MATLAB) with independents')
        fopt, qopt, exitflag, niter, status, sol_time = call_cvx_ind(problem, eng, printout=printout)
    else:
        print('Calling LP-Optimisation via CVX (MATLAB) with NO independents')
        fopt, qopt, exitflag, niter, status, sol_time = call_cvx(problem, eng, printout=printout)

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

    # print(problem.q)
    # from compas_plotters import MeshPlotter
    # plotter = MeshPlotter(form)
    # plotter.draw_edges(text={edge: round(form.edge_attribute(edge, 'q'), 2) for edge in form.edges()})
    # plotter.show()

    form.attributes['loadpath'] = form.loadpath()
    reactions(form)

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


def call_cvx(problem, eng, printout=False):

    # q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i, eng = args_cvx

    start_time1 = time.time()

    eng.workspace['m'] = problem.m

    eng.workspace['C'] = matlab.double(problem.C.toarray().tolist())
    eng.workspace['Ci'] = matlab.double(problem.Ci.toarray().tolist())
    eng.workspace['Cb'] = matlab.double(problem.Cb.toarray().tolist())
    eng.workspace['E'] = matlab.double(problem.E.tolist())

    # load_csr_matrix_to_matlab(problem.E, 'E', eng)

    # print(problem.qmin.shape, problem.qmax.shape, problem.q.shape)

    eng.workspace['xt'] = matlab.double(problem.X[:, 0].reshape(-1, 1).transpose().tolist())
    eng.workspace['yt'] = matlab.double(problem.X[:, 1].reshape(-1, 1).transpose().tolist())
    eng.workspace['xb'] = matlab.double(problem.X[:, 0][problem.fixed].reshape(-1, 1).tolist())
    eng.workspace['yb'] = matlab.double(problem.X[:, 1][problem.fixed].reshape(-1, 1).tolist())
    eng.workspace['pz'] = matlab.double(problem.P[:, 2][problem.free].reshape(-1, 1).tolist())
    eng.workspace['p'] = matlab.double(problem.P[problem.free, :2].flatten('F').reshape((-1, 1)).tolist())

    eng.workspace['qmax'] = matlab.double(problem.qmax.reshape(-1, 1).tolist())
    eng.workspace['qmin'] = matlab.double(problem.qmin.reshape(-1, 1).tolist())

    # If intended to save the MATLAB engine
    # eng.save('/Users/mricardo/Documents/MATLAB/data.mat', nargout=0)

    start_time = time.time()

    if printout:
        eng.cvx_begin(nargout=0)
    else:
        eng.cvx_begin('quiet', nargout=0)
    eng.variable('q(double(m))', nargout=0)
    # if objective == 'loadpath':
    eng.minimize('matrix_frac(pz, -(transpose(Ci)*diag(q)*Ci)) - xt*transpose(C)*diag(q)*Cb*xb - yt*transpose(C)*diag(q)*Cb*yb', nargout=0)
    # if objective == 'feasibility':
    # eng.minimize('1', nargout=0)
    eng.eval('q >= qmin', nargout=0)
    eng.eval('q <= qmax', nargout=0)
    eng.eval('E * q == p', nargout=0)
    eng.cvx_end(nargout=0)

    sol_time = time.time() - start_time

    fopt = eng.workspace['cvx_optval']
    qopt = array(eng.workspace['q'])
    status = eng.workspace['cvx_status']
    niter = eng.workspace['cvx_slvitr']

    elapsed_time = time.time() - start_time1
    print('Elapsed time on LP:', elapsed_time)

    if status != 'Infeasible':
        exitflag = 0
    else:
        exitflag = 1

    return fopt, qopt, exitflag, niter, status, sol_time


def call_cvx_ind(problem, eng, printout=True):

    # q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i, eng = args_cvx

    ind_ = [x+1 for x in problem.ind]
    dep_ = [x+1 for x in problem.dep]

    start_time = time.time()

    # load_csr_matrix_to_matlab(problem.C, 'C', eng)
    # load_csr_matrix_to_matlab(problem.Ci, 'Ci', eng)
    # load_csr_matrix_to_matlab(problem.Cb, 'Cb', eng)
    # load_csr_matrix_to_matlab(problem.Edinv, 'Edinv', eng)
    # load_csr_matrix_to_matlab(problem.Ei, 'Ei', eng)

    eng.workspace['C'] = matlab.double(problem.C.toarray().tolist())
    eng.workspace['Ci'] = matlab.double(problem.Ci.toarray().tolist())
    eng.workspace['Cb'] = matlab.double(problem.Cb.toarray().tolist())
    eng.workspace['Edinv'] = matlab.double(problem.Edinv.toarray().tolist())
    eng.workspace['Ei'] = matlab.double(problem.Ei.toarray().tolist())

    eng.workspace['m'] = problem.m

    eng.workspace['dep'] = matlab.double(dep_)
    eng.workspace['ind'] = matlab.double(ind_)
    eng.workspace['p'] = matlab.double(problem.P[problem.free, :2].flatten('F').reshape((-1, 1)).tolist())
    eng.workspace['xt'] = matlab.double(problem.X[:, 0].reshape(-1, 1).transpose().tolist())
    eng.workspace['yt'] = matlab.double(problem.X[:, 1].reshape(-1, 1).transpose().tolist())
    eng.workspace['xb'] = matlab.double(problem.X[:, 0][problem.fixed].reshape(-1, 1).tolist())
    eng.workspace['yb'] = matlab.double(problem.X[:, 1][problem.fixed].reshape(-1, 1).tolist())
    eng.workspace['pz'] = matlab.double(problem.P[:, 2][problem.free].reshape(-1, 1).tolist())
    eng.workspace['qmax'] = matlab.double(problem.qmax.reshape(-1, 1).tolist())
    eng.workspace['qmin'] = matlab.double(problem.qmin.reshape(-1, 1).tolist())

    sol_time0 = time.time()

    if printout:
        eng.cvx_begin(nargout=0)
    else:
        eng.cvx_begin('quiet', nargout=0)
    eng.variable('q(double(m))', nargout=0)
    eng.minimize('matrix_frac(pz, -transpose(Ci)*diag(q)*Ci) - xt*transpose(C)*diag(q)*Cb*xb - yt*transpose(C)*diag(q)*Cb*yb', nargout=0)  # with compression negative
    eng.eval('q >= qmin', nargout=0)
    eng.eval('q <= qmax', nargout=0)
    eng.eval('q(dep) == Edinv*(Ei*q(ind) - p)', nargout=0)
    eng.cvx_end(nargout=0)

    sol_time = time.time() - sol_time0

    fopt = eng.workspace['cvx_optval']
    status = eng.workspace['cvx_status']
    qopt = array(eng.workspace['q'])
    niter = int(eng.workspace['cvx_slvitr'])

    elapsed_time = time.time() - start_time
    print('Elapsed time on LP:', elapsed_time)

    # If intended to save the MATLAB engine
    eng.save('/Users/mricardo/Documents/MATLAB/data.mat', nargout=0)

    if status != 'Infeasible':
        exitflag = 0
    else:
        exitflag = 1

    return fopt, qopt, exitflag, niter, status, sol_time


def load_csr_matrix_to_matlab(mat, name, eng):

    (m, n) = mat.shape
    s = mat.data
    i = mat.tocoo().row
    j = mat.indices

    i_ = [x+1 for x in i]
    j_ = [x+1 for x in j]

    eng.workspace['m_'] = matlab.double([m])
    eng.workspace['n_'] = matlab.double([n])
    eng.workspace['s'] = matlab.double(s.tolist())
    eng.workspace['i_'] = matlab.double(i_)
    eng.workspace['j_'] = matlab.double(j_)

    eng.eval(name + ' = sparse(i_, j_, s, m_, n_, m_*n_);', nargout=0)

    return
