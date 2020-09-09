
from compas_tno.algorithms import zlq_from_qid
from compas_tno.algorithms import zlq_from_q

try:
    import matlab.engine
except:
    print('Error: MATLAB Python Engine Not Available!')

from numpy import array

import time

from compas.utilities import geometric_key
from compas_tno.algorithms import reactions

from compas_tno.problems import initialise_problem

__all__ = [
    'run_optimisation_MATLAB',
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
    find_inds = optimiser.data['find_inds']
    objective = optimiser.data['objective']
    printout = optimiser.data['printout']
    plot = optimiser.data['plot']

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


def run_loadpath_from_form_MATLAB(form, find_inds=True, qmax=3000, printout=False):
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

    # Initiate Matlab Engine

    future = matlab.engine.start_matlab(background=True)
    eng = future.result()

    objective = 'loadpath'
    indset = None
    find_inds = True
    plot = printout

    k_i = form.key_index()
    i_uv = form.index_uv()

    args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y = args[:31]
    args_cvx = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i)
    indset = form.attributes['indset']

    i_k = form.index_key()
    i_uv = form.index_uv()

    q, ind = args_cvx[:2]
    args_cvx_ = list(args_cvx)
    args_cvx_.append(eng)
    args_cvx = tuple(args_cvx_)
    args = args_cvx[:22]

    call_and_output_CVX_MATLAB(form, find_inds, indset, printout, plot, objective, args_cvx, args, ind, i_uv, i_k)

    return form


def call_and_output_CVX_MATLAB(form, find_inds, indset, printout, plot, objective, args_cvx, args, ind, i_uv, i_k):

    if find_inds or indset:
        print('Calling LP-Optimisation via CVX (MATLAB) with independents')
        fopt, qopt, exitflag, niter, status, sol_time = call_cvx_ind(objective, args_cvx, printout=printout)
        z, _, q, q_ = zlq_from_qid(qopt[ind], args)
    else:
        print('Calling LP-Optimisation via CVX (MATLAB) with NO independents')
        fopt, qopt, exitflag, niter, status, sol_time = call_cvx(objective, args_cvx, printout=printout)
        z, _, q, q_ = zlq_from_q(qopt, args)

    # Update form diagram optimised

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
    for u, v in form.edges_where({'_is_edge': True}):
        if form.edge_attribute((u, v), 'is_symmetry') is False:
            qi = form.edge_attribute((u, v), 'q')
            li = form.edge_length(u, v)
            lp += abs(qi) * li**2
    form.attributes['loadpath'] = float(lp)
    reactions(form, plot=plot)

    summary = True

    if printout or summary:
        print('\n' + '-' * 50)
        print('LOADPATH OPTIMISATION WITH CVX (MATLAB)')
        print('status    :', status)
        print('fopt (lp) : {0:.3f}'.format(fopt))
        print('n-iter    : {0}'.format(niter))
        print('q range   : {0:.3f} : {1:.3f}'.format(min(qopt)[0], max(qopt)[0]))
        print('sol. time : {0:.3f} sec'.format(sol_time))
        print('-' * 50 + '\n')

    return fopt, exitflag


def call_cvx(objective, args_cvx, printout=False):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i, eng = args_cvx

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

    # If intended to save the MATLAB engine
    # eng.save('/Users/mricardo/Documents/MATLAB/data.mat', nargout=0)

    start_time = time.time()

    if printout:
        eng.cvx_begin(nargout=0)
    else:
        eng.cvx_begin('quiet', nargout=0)
    eng.variable('q(double(m))', nargout=0)
    if objective == 'loadpath':
        eng.minimize('matrix_frac(pz,(Cit*diag(q)*Ci)) + xt*transpose(C)*diag(q)*Cf*xb + yt*transpose(C)*diag(q)*Cf*yb', nargout=0)
    if objective == 'feasibility':
        eng.minimize('1', nargout=0)
    eng.eval('q >= -0.001', nargout=0)
    eng.eval('q <= qmax', nargout=0)
    eng.eval('E * q == 0.0', nargout=0)
    eng.cvx_end(nargout=0)

    sol_time = time.time() - start_time

    fopt = eng.workspace['cvx_optval']
    qopt = array(eng.workspace['q'])
    status = eng.workspace['cvx_status']
    niter = eng.workspace['cvx_slvitr']

    if status != 'Infeasible':
        exitflag = 0
    else:
        exitflag = 1

    return fopt, qopt, exitflag, niter, status, sol_time


def call_cvx_ind(fobj, args_cvx, printout=False):

    q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, qmax, i_uv, k_i, eng = args_cvx

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

    start_time = time.time()

    if printout:
        eng.cvx_begin(nargout=0)
    else:
        eng.cvx_begin('quiet', nargout=0)
    eng.variable('q(double(m))', nargout=0)
    eng.minimize('matrix_frac(pz,(Cit*diag(q)*Ci)) + xt*transpose(C)*diag(q)*Cf*xb + yt*transpose(C)*diag(q)*Cf*yb', nargout=0)
    eng.eval('q >= 0.0', nargout=0)
    eng.eval('q <= qmax', nargout=0)
    eng.eval('q(dep) == - Edinv*(p - Ei*q(ind))', nargout=0)
    eng.cvx_end(nargout=0)

    sol_time = time.time() - start_time

    fopt = eng.workspace['cvx_optval']
    status = eng.workspace['cvx_status']
    qopt = array(eng.workspace['q'])
    niter = int(eng.workspace['cvx_slvitr'])

    if status != 'Infeasible':
        exitflag = 0
    else:
        exitflag = 1

    return fopt, qopt, exitflag, niter, status, sol_time
