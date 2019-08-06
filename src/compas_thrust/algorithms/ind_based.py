from numpy import abs
from numpy import argmin
from numpy import array
from numpy import float64
from numpy import dot
from numpy import hstack
from numpy import isnan
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import sqrt
from numpy import sum
from numpy import vstack
from numpy import zeros
from numpy import ones
from numpy import append
from numpy.linalg import pinv
from numpy.linalg import matrix_rank
from numpy.random import rand
from numpy.random import randint

from scipy.linalg import svd
from scipy.optimize import fmin_slsqp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import time

from compas_tna.diagrams import FormDiagram

from compas.numerical import connectivity_matrix
from compas.numerical import devo_numpy
from compas.numerical import equilibrium_matrix
from compas.numerical import normrow
from compas.numerical import nonpivots
from compas.utilities import geometric_key

from compas.geometry import intersection_line_line
from compas.geometry import is_point_on_segment

from compas_thrust.plotters.plotters import plot_form

from compas_thrust.algorithms.equilibrium import zlq_from_qid
from compas_thrust.algorithms.equilibrium import reactions

from multiprocessing import Pool
from random import shuffle
from copy import deepcopy


__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'optimise_single',
    'optimise_multi',
    'run_replicate',
    'find_independents',
    'initialize_problem',
]


def optimise_single(form, solver='devo', polish='slsqp', qmin=1e-6, qmax=10, population=300, generations=500,
                    printout=10, tol=0.01, plot=False, frange=[], indset=None, tension=False, planar=False,
                    opt_max=False, t=None, random_ind=False, target=False, use_bounds=None, bounds_width = 5.0,
                    summary = True, objective='loadpath', buttress = False):

    """ Finds the optimised load-path for a FormDiagram.

    Parameters
    ----------
    form : obj
        The FormDiagram.
    solver : str
        Differential Evolution 'devo' or Genetic Algorithm 'ga' evolutionary solver to use.
    polish : str
        'slsqp' polish or None.
    qmin : float
        Minimum qid value.
    qmax : float
        Maximum qid value.
    population : int
        Number of agents for the evolution solver.
    generations : int
        Number of generations for the evolution solver.
    printout : int
        Frequency of print output to the terminal.
    tol : float
        Tolerance on horizontal force balance.
    plot : bool
        Plot progress of the evolution.
    frange : list
        Minimum and maximum function value to plot.
    indset : list
        Independent set to use.
    tension : bool
        Allow tension edge force densities (experimental).
    planar : bool
        Only consider the x-y plane.

    Returns
    -------
    float
        Optimum load-path value.
    list
        Optimum qids

    """

    if printout:
        print('\n' + '-' * 50)
        print('Load-path optimisation started')
        print('-' * 50)

    # Mapping

    k_i  = form.key_index()
    i_k  = form.index_key()
    i_uv = form.index_uv()
    uv_i = form.uv_index()

    # Vertices and edges

    n     = form.number_of_vertices()
    m     = form.number_of_edges()
    fixed = [k_i[key] for key in form.fixed()]
    rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
    sym   = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
    free  = list(set(range(n)) - set(fixed) - set(rol))

    # Constraints

    lb_ind = []
    ub_ind = []
    lb = []
    ub = []
    for key, vertex in form.vertex.items():
        if vertex.get('lb', None):
            lb_ind.append(k_i[key])
            lb.append(vertex['lb'])
        if vertex.get('ub', None):
            ub_ind.append(k_i[key])
            ub.append(vertex['ub'])
    lb = array(lb)
    ub = array(ub)
    lb.shape = (len(lb),1)
    ub.shape = (len(ub),1)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))
    s   = zeros((n, 1))
    w   = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i]  = vertex.get('x')
        y[i]  = vertex.get('y')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i]  = vertex.get('target', 0)
        w[i] = vertex.get('weight', 1.0)

    Wfree = diags(w[free].flatten())
    xy = xyz[:, :2]

    # Temporary anchors location

    anchors = []

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cf  = C[:, fixed]
    Ct = C.transpose()
    Cit = Ci.transpose()
    E   = equilibrium_matrix(C, xy, free, 'csr').toarray()
    uvw = C.dot(xyz)
    U   = uvw[:, 0]
    V   = uvw[:, 1]

    # Independent and dependent branches

    if indset:
        ind = []
        for u, v in form.edges():
            if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
                ind.append(uv_i[(u, v)])
    else:
        _, s, _ = svd(E)
        ind = []
        if random_ind is False:
            # ind = nonpivots(sympy.Matrix(E).rref()[0].tolist())
            ind = find_independents(E)
        else:
            ind = randint(0,m,m - len(s)).tolist()

    k   = len(ind)
    dep = list(set(range(m)) - set(ind))

    for u, v in form.edges():
        form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

    if printout:
        print('Form diagram has {0} (RREF) or {1} (SVD) independent branches '.format(len(ind), m - len(s)))
        print('Shape Equilibrium Matrix: ', E.shape)
        print('Found {0} independents'.format(k))

    # Set-up

    lh     = normrow(C.dot(xy))**2
    Edinv  = -csr_matrix(pinv(E[:, dep]))
    Ei     = E[:, ind]
    p      = vstack([px[free], py[free]])
    q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

    if objective == 'loadpath':
        t = None

    if buttress:
        b = {}
        for key in form.vertices_where({'is_fixed': True}):
            try:
                b[k_i[key]] = form.get_vertex_attributes(key, 'b')
            except:
                pass
                print(b)
    else:
        b = None
    
    print(b)
    
    try:
        joints = form.attributes['joints']
    except:
        joints = None

    print(joints)

    args   = (q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i)

    if use_bounds:
        bounds = []
        for i in ind:
            u,v = i_uv[i]
            qi = form.get_edge_attribute((u,v), 'q')
            if form.get_edge_attribute((u,v), 'is_symmetry') is True:
                bounds.append([qmin, qmax])
            else:
                if qi > bounds_width + qmin:
                    bounds.append([qi-bounds_width,qi+bounds_width])
                else:
                    bounds.append([qmin,qi+bounds_width])
    else:
        bounds = [[qmin, qmax]] * k
    if t:
        bounds.append([-t,t])
    

    # Horizontal checks

    checked = True

    if tol > 0:
        for i in range(10**3):
            q[ind, 0] = rand(k) * qmax
            q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
            Rx = Cit.dot(U * q.ravel()) - px[free].ravel()
            Ry = Cit.dot(V * q.ravel()) - py[free].ravel()
            R  = max(sqrt(Rx**2 + Ry**2))
            if R > tol:
                checked = False
                break

    if checked:

        # Define Objective

        if objective=='loadpath':
            fdevo, fslsqp, fieq = _fint, _fint_, _fieq
        if objective=='target':
            fdevo, fslsqp, fieq = fbf, fbf, _fieq
        if objective == 'min':
            fdevo, fslsqp, fieq = fmin, fmin, _fieq_bounds
        if objective=='max':
            fdevo, fslsqp, fieq = fmax, fmax, _fieq_bounds
        if objective=='bounds':
            fdevo, fslsqp, fieq = fbounds, fbounds, _fieq_bounds

        # Solve

        if solver == 'devo':
            fopt, qopt = _diff_evo(fdevo, bounds, population, generations, printout, plot, frange, args)

        else:
            qi = q[ind]
            if t:
                qopt = append(qi,t)
            else:
                qopt = qi.reshape(k,)
            fopt = fslsqp(qopt,*args)

        if polish == 'slsqp':
            fopt_, qopt_ = _slsqp(fslsqp, qopt, bounds, printout, fieq, args)

            if t:
                qopt_, t_slsqp = qopt_[:k], qopt_[-1]
                qopt, t = qopt[:k], qopt[-1]
            q1 = zlq_from_qid(qopt_, args)[2]
            if fopt_ < fopt:
                if (min(q1) > -0.001 and not tension) or tension:
                    fopt, qopt = fopt_, qopt_
                    if t:
                        t = t_slsqp
        if t:
            form.attributes['offset'] = t

        z, _, q, q_ = zlq_from_qid(qopt, args)

        # Unique key

        gkeys = []
        for i in ind:
            u, v = i_uv[i]
            gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
        form.attributes['indset'] = gkeys

        # Update FormDiagram

        for i in range(n):
            key = i_k[i]
            form.set_vertex_attribute(key=key, name='z', value=float(z[i]))

        for c, qi in enumerate(list(q_.ravel())):
            u, v = i_uv[c]
            form.set_edge_attribute((u, v), 'q', float(qi))

        # Relax

        # q    = array([attr['q'] for u, v, attr in form.edges(True)])
        # Q    = diags(q)
        # CitQ = Cit.dot(Q)
        # Di   = CitQ.dot(Ci)
        # Df   = CitQ.dot(Cf)
        # bx   = px - Df.dot(x[fixed])
        # by   = py - Df.dot(y[fixed])
        # # bz   = pz - Df.dot(z[fixed])
        # x[free, 0] = spsolve(Di, bx)
        # y[free, 0] = spsolve(Di, by)
        # # z[free, 0] = spsolve(Di, bz)

        # for i in range(n):
        #     form.set_vertex_attributes(key=i_k[i], names='xyz', values=[float(j) for j in [x[i], y[i], z[i]]])

        fopt = 0
        for u, v in form.edges():
            if form.get_edge_attribute((u, v), 'is_symmetry') is False:
                qi = form.get_edge_attribute((u, v), 'q')
                li = form.edge_length(u, v)
                fopt += abs(qi) * li**2

        form.attributes['loadpath'] = fopt

        CfQ = Cf.transpose().dot(diags(q.flatten()))

        if objective == 'min':
            fopt = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed]).flatten()[0]
            [[fopt]] = fopt

        if objective == 'max':
            fopt = 1 / (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])[0][0]
            [[fopt]] = fopt

        reactions(form, plot = plot)

        if summary:
            print('\n' + '-' * 50)
            print('qid range : {0:.3f} : {1:.3f}'.format(min(qopt), max(qopt)))
            print('q range   : {0:.3f} : {1:.3f}'.format(min(q), max(q)))
            print('fopt      : {0:.3f}'.format(fopt))
            if t:
                print('t      : {0:.3f}'.format(t))
            print('-' * 50 + '\n')

        return fopt, qopt

    else:

        if printout:
            print('Horizontal equillibrium checks failed')

        return None, None

def _fint(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints = args

    qid, t = qid[:k], qid[-1]
    z, l2, q, q_ = zlq_from_qid(qid, args)
    f = dot(abs(q.transpose()), l2)

    if isnan(f):
        return 10**10

    else:

        if not tension:
            f += sum((q[q < 0] - 5)**4)

        Rx = Cit.dot(U * q_.ravel()) - px[free].ravel()
        Ry = Cit.dot(V * q_.ravel()) - py[free].ravel()
        Rh = Rx**2 + Ry**2
        Rm = max(sqrt(Rh))
        if Rm > tol:
            f += sum(Rh - tol + 5)**4

        return f

def fmin(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

    qid, t = qid[:k], qid[-1]
    z, l2, q, q_ = zlq_from_qid(qid, args)
    CfQ = Cf.transpose().dot(diags(q.flatten()))
    f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])
    # f +=  pz.transpose().dot(z[free])/1000

    if isnan(f) == True or any(qid) == False:
        return 10**10
    else:

        if not tension:
            f += sum((q[q < 0] - 10)**4)
        Rx = Ct.dot(U * q_.ravel()) - px.ravel()
        Ry = Ct.dot(V * q_.ravel()) - py.ravel()
        Rh = Rx**2 + Ry**2
        Rm = max(sqrt(Rh[free]))
        
        if Rm > tol:
            f += sum(Rh - tol + 5)**4

        if lb_ind:
            if t is not 0.0:
                lb = lb - t
            z_lb    = z[lb_ind]
            log_lb  = z_lb < lb
            diff_lb = z_lb[log_lb] - lb[log_lb]
            pen_lb  = sum(abs(diff_lb) + 5)**4
            f += pen_lb

        if ub_ind:
            if t is not 0.0:
                ub = ub - t
            z_ub    = z[ub_ind]
            log_ub  = z_ub > ub
            diff_ub = z_ub[log_ub] - ub[log_ub]
            pen_ub  = sum(abs(diff_ub) + 5)**4
            f += pen_ub

        if b:
            W = C.dot(z)[:,0]
            Rz = Ct.dot(W * q_.ravel()) - pz.ravel()
            for key in b:
                scl = t/Rz[key]
                x_comp = abs(scl * Rx[key])
                y_comp = abs(scl * Ry[key])
                if x_comp > abs(b[key][0][0]):
                    f += (x_comp - b[key][0][0] + 5) **2
                if y_comp > abs(b[key][0][1]):
                    f += (y_comp - b[key][0][1] + 5) **2

        if joints:
            for jt in joints:
                limit = [jt[0],jt[1]]
                u, v = k_i[jt[2][0]] , k_i[jt[2][1]]
                x_ = list(x)
                y_ = list(y)
                z_ = list(z)
                thrust = [[x_[u],y_[u],z_[u]+t],[x_[v],y_[v],z_[v]+t]]
                pt = intersection_line_line(limit,thrust)[0]
                if pt == None or is_point_on_segment(pt,limit,tol=1e-6) == False:
                    f += 10 ** 5
                
        return f

def _fmin(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints = args

    qid, t = qid[:k], qid[-1]
    z, l2, q, q_ = zlq_from_qid(qid, args)
    CfQ = Cf.transpose().dot(diags(q.flatten()))
    f = (CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed])

    if isnan(f) == True or any(qid) == False:
        return 10**10

    return f

def fbounds(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints = args
    qid, t = qid[:k], qid[-1]
    z, l2, q, q_ = zlq_from_qid(qid, args)
    f = 0

    if isnan(f):
        return 10**10

    else:

        if not tension:
            f += sum((q[q < 0] - 5)**4)

        Rx = Cit.dot(U * q_.ravel()) - px[free].ravel()
        Ry = Cit.dot(V * q_.ravel()) - py[free].ravel()
        Rh = Rx**2 + Ry**2
        Rm = max(sqrt(Rh))
        if Rm > tol:
            f += sum(Rh - tol + 5)**4

        if lb_ind:
            if t is not 0.0:
                lb = lb - t
            z_lb    = z[lb_ind]
            log_lb  = z_lb < lb
            diff_lb = z_lb[log_lb] - lb[log_lb]
            pen_lb  = sum(abs(diff_lb) + 5)**4
            f += pen_lb

        if ub_ind:
            if t is not 0.0:
                ub = ub - t
            z_ub    = z[ub_ind]
            log_ub  = z_ub > ub
            diff_ub = z_ub[log_ub] - ub[log_ub]
            pen_ub  = sum(abs(diff_ub) + 5)**4
            f += pen_ub

        if b:
            for key in b:
                scl = t/Rz[key]
                x_comp = abs(scl * Rx[key])
                y_comp = abs(scl * Ry[key])
                # print(b[key],scl,x_comp, y_comp)
                if x_comp > abs(b[key][0][0]):
                    f += (x_comp - b[key][0][0] + 5) **2
                if y_comp > abs(b[key][0][1]):
                    f += (y_comp - b[key][0][1] + 5) **2

        return f

def fmax(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints = args

    qid, t = qid[:k], qid[-1]
    z, l2, q, q_ = zlq_from_qid(qid, args)
    CfQ = Cf.transpose().dot(diags(q_.flatten()))
    f = 1 / ((CfQ.dot(U[:,newaxis])).transpose().dot(x[fixed]) + (CfQ.dot(V[:,newaxis])).transpose().dot(y[fixed]))

    if isnan(f) == True or any(qid) == False:
        return 10**10

    else:

        if not tension:
            f += sum((q[q < 0] - 5)**4)
        W = C.dot(z)[:,0]
        Rx = Ct.dot(U * q_.ravel()) - px.ravel()
        Ry = Ct.dot(V * q_.ravel()) - py.ravel()
        Rz = Ct.dot(W * q_.ravel()) - pz.ravel()
        Rh = Rx**2 + Ry**2
        Rm = max(sqrt(Rh[free]))
        if Rm > tol:
            f += sum(Rh - tol + 5)**4

        if lb_ind:
            if t is not 0.0:
                lb = lb - t
            z_lb    = z[lb_ind]
            log_lb  = z_lb < lb
            diff_lb = z_lb[log_lb] - lb[log_lb]
            pen_lb  = sum(abs(diff_lb) + 2)**4
            f += pen_lb

        if ub_ind:
            if t is not 0.0:
                ub = ub - t
            z_ub    = z[ub_ind]
            log_ub  = z_ub > ub
            diff_ub = z_ub[log_ub] - ub[log_ub]
            pen_ub  = sum(abs(diff_ub) + 2)**4
            f += pen_ub

        if b:
            for key in b:
                scl = t/Rz[key]
                x_comp = abs(scl * Rx[key])
                y_comp = abs(scl * Ry[key])
                # print(b[key],scl,x_comp, y_comp)
                if x_comp > abs(b[key][0][0]):
                    f += (x_comp - b[key][0][0] + 5) **2
                if y_comp > abs(b[key][0][1]):
                    f += (y_comp - b[key][0][1] + 5) **2

        return f

def fbf(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b = args

    qid, t = qid[:k], qid[-1]
    z, l2, q, q_ = zlq_from_qid(qid, args)
    z = z + t
    z_s = abs(z[free] - s[free])
    f = Wfree.dot(z_s.transpose().dot(z_s))
    f = f[0,0]

    if not tension:
            f += sum((q[q < 0] - 5)**4)

    if isnan(f):
        return 10**10

    else:
        return f

def _fint_(qid, *args):

    z, l2, q, q_ = zlq_from_qid(qid, args)
    f = dot(abs(q.transpose()), l2)

    if isnan(f):
        return 10**10

    return f

def _fieq(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints = args

    q[ind, 0], t = qid[:k], qid[-1]
    q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
    q_ = 1 * q
    q[sym] *= 0

    Rx = Cit.dot(U * q_.ravel()) - px[free].ravel()
    Ry = Cit.dot(V * q_.ravel()) - py[free].ravel()
    Rh = Rx**2 + Ry**2
    Rm = max(sqrt(Rh))

    pen_lb = 0.0
    pen_ub = 0.0

    if not tension:
        return hstack([q.ravel() + 10**(-5), tol - Rm])
        # return hstack([Rm + pen_lb + pen_ub])
    return [tol - Rm - pen_lb - pen_ub]

def _fieq_bounds(qid, *args):

    q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b, joints, i_uv, k_i = args

    q[ind, 0], t = qid[:k], qid[-1]
    q[dep] = -Edinv.dot(p - Ei.dot(q[ind]))
    q_ = 1 * q
    q[sym] *= 0
    
    rbut = 0.0

    W = C.dot(z)[:,0]
    Rx = Ct.dot(U * q_.ravel()) - px.ravel()
    Ry = Ct.dot(V * q_.ravel()) - py.ravel()
    Rz = Ct.dot(W * q_.ravel()) - pz.ravel()
    Rh = Rx**2 + Ry**2
    Rm = max(sqrt(Rh[free]))

    pen_lb = 0.0
    pen_ub = 0.0
    joint = 0.0

    if lb_ind:
        if t is not 0.0:
            lb = lb - t
        z_lb    = z[lb_ind]
        log_lb  = z_lb < lb
        diff_lb = z_lb[log_lb] - lb[log_lb]
        pen_lb  = sum(abs(diff_lb) + 2)**4

    if ub_ind:
        if t is not 0.0:
            ub = ub - t
        z_ub    = z[ub_ind]
        log_ub  = z_ub > ub
        diff_ub = z_ub[log_ub] - ub[log_ub]
        pen_ub  = sum(abs(diff_ub) + 2)**4

    if b:
            for key in b:
                scl = t/Rz[key]
                x_comp = abs(scl * Rx[key])
                y_comp = abs(scl * Ry[key])
                if x_comp > abs(b[key][0][0]):
                    rbut += (x_comp - b[key][0][0] + 5) **2
                if y_comp > abs(b[key][0][1]):
                    rbut += (y_comp - b[key][0][1] + 5) **2

    if joints:
        for jt in joints:
            limit = [jt[0],jt[1]]
            u, v = k_i[jt[2][0]] , k_i[jt[2][1]]
            x_ = list(x)
            y_ = list(y)
            z_ = list(z)
            thrust = [[x_[u],y_[u],z_[u]+t],[x_[v],y_[v],z_[v]+t]]
            pt = intersection_line_line(limit,thrust)[0]
            if pt == None or is_point_on_segment(pt,limit,tol=1e-6) == False:
                joint += 10 ** 5

    if not tension:
        return hstack([Rm + pen_lb + pen_ub])
        # return hstack([q.ravel() + 10**(-5), tol - Rm - pen_lb - pen_ub])
        # return hstack([Rm + pen_lb + pen_ub])
    return [tol - Rm - pen_lb - pen_ub - rbut - joint]

def _slsqp(fn, qid0, bounds, printout, fieq, args):

    pout = 2 if printout else 0
    opt  = fmin_slsqp(fn, qid0, args=args, disp=pout, bounds=bounds, full_output=1, iter=1000, f_ieqcons=fieq)

    return opt[1], opt[0]

def _diff_evo(fn, bounds, population, generations, printout, plot, frange, args):

    return devo_numpy(fn=fn, bounds=bounds, population=population, generations=generations, printout=printout,
                      plot=plot, frange=frange, args=args)

def _worker(data):

    try:

        i, form, save_figs, qmin, qmax, population, generations, simple, tension, tol = data
        fopt, qopt = optimise_single(form, qmin=qmin, qmax=qmax, population=population, generations=generations,
                                     printout=0, tension=tension, tol=tol)

        if isnan(fopt):
            fopt = 10**10

        print('Trial: {0} - Optimum: {1:.1f}'.format(i, fopt))

        if save_figs:
            plotter = plot_form(form, radius=0.1, fix_width=False, simple=simple)
            plotter.save('{0}trial_{1}-fopt_{2:.6f}.png'.format(save_figs, i, fopt))
            del plotter

        return fopt, form

    except:

        print('Trial: {0} - FAILED'.format(i))

        return 10**10, None

def optimise_multi(form, trials=10, save_figs='', qmin=0.001, qmax=5, population=300, generations=500, simple=False,
                   tension=False, tol=0.001, opt_max=False):

    """ Finds the optimised load-path for multiple d FormDiagrams.

    Parameters
    ----------
    form : obj
        FormDiagram to analyse.
    trials : int
        Number of trials to perform.
    save_figs : str
        Directory to save plots.
    qmin : float
        Minimum qid value.
    qmax : float
        Maximum qid value.
    population : int
        Number of agents for the evolution solver.
    generations : int
        Number of generations for the evolution solver.
    simple : bool
        Simple red and blue colour plotting.
    tension : bool
        Allow tension edge force densities (experimental).
    tol : float
        Tolerance on horizontal force balance.

    Returns
    -------
    list
        Optimum load-path for each trial.
    list
        Each resulting trial FormDiagram.
    int
        Index of the optimum.

    """

    data = [(i, _form(form), save_figs, qmin, qmax, population, generations, simple, tension, tol)
            for i in range(trials)]

    fopts, forms = zip(*Pool().map(_worker, data))
    best = argmin(fopts)
    if opt_max:
        best = argmax(fopts)

    print('Best: {0} - fopt {1:.1f}'.format(best, fopts[best]))

    return fopts, forms, best

def run_replicate(form, file_complete, delete_face=False, plots=False, save=True):

    print('\nOptimization and Replicate of Form ')

    form = FormDiagram.from_json(file)
    form = _form(form)

    plot_form(form, radius=0.05).show()

    for i in range(1):
            print('Optimisation trial {0}.'.format(i))
            form = _form(form)
            fopt, qopt = optimise_single(form, qmax=5, population=200, generations=200, printout=10, tol=0.01)
            if fopt is not None and fopt < 2000:
                print('Optimisation found a result after {0} trials.'.format(i))
                break

    form.to_json(file)
    form_ = replicate(form,file_complete)

    if plots:
        plot_form(form, radius=0.05).show()
        plot_form(form_, radius=0.05).show()

    z_from_form(form_)
    oveview_forces(form_)

    if save:
        form.to_json(file)
        form_.to_json(file_complete)

    return form, form_

def find_independents(E):

    _, m = E.shape
    Etemp = E[:,[0]]
    ind = []

    for i in range(1,m):
        Etest = hstack([Etemp,E[:,[i]]])
        _ , ncol = Etest.shape
        if matrix_rank(Etest) < ncol:
            ind.append(i)
        else:
            Etemp = Etest

    return ind

def initialize_problem(form, indset = None, printout = None):

    # Mapping

    k_i  = form.key_index()
    i_k  = form.index_key()
    i_uv = form.index_uv()
    uv_i = form.uv_index()

    # Vertices and edges

    n     = form.number_of_vertices()
    m     = form.number_of_edges()
    fixed = [k_i[key] for key in form.fixed()]
    rol   = [k_i[key] for key in form.vertices_where({'is_roller': True})]
    edges = [(k_i[u], k_i[v]) for u, v in form.edges()]
    sym   = [uv_i[uv] for uv in form.edges_where({'is_symmetry': True})]
    free  = list(set(range(n)) - set(fixed) - set(rol))

    # Constraints

    lb_ind = []
    ub_ind = []
    lb = []
    ub = []
    for key, vertex in form.vertex.items():
        if vertex.get('lb', None):
            lb_ind.append(k_i[key])
            lb.append(vertex['lb'])
        if vertex.get('ub', None):
            ub_ind.append(k_i[key])
            ub.append(vertex['ub'])
    lb = array(lb)
    ub = array(ub)
    lb.shape = (len(lb),1)
    ub.shape = (len(ub),1)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))
    s   = zeros((n, 1))
    w   = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i]  = vertex.get('x')
        y[i]  = vertex.get('y')
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        s[i]  = vertex.get('target', 0)
        w[i] = vertex.get('weight', 1.0) # weight used in case of fiting...

    Wfree = diags(w[free].flatten())
    xy = xyz[:, :2]
    z = xyz[:, 2]

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cf  = C[:, fixed]
    Ct = C.transpose()
    Cit = Ci.transpose()
    E   = equilibrium_matrix(C, xy, free, 'csr').toarray()
    uvw = C.dot(xyz)
    U   = uvw[:, 0]
    V   = uvw[:, 1]

    print('Equilibrium Matrix Shape: ', E.shape)
    start_time = time.time()

    # Independent and dependent branches

    if indset:
        ind = []
        for u, v in form.edges():
            if geometric_key(form.edge_midpoint(u, v)[:2] + [0]) in indset:
                ind.append(uv_i[(u, v)])
    else:
        _, s, _ = svd(E)
        ind = find_independents(E)



    k   = len(ind)
    dep = list(set(range(m)) - set(ind))
    elapsed_time = time.time() - start_time
    print('Found {0} independents'.format(k))
    print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

    for u, v in form.edges():
        form.set_edge_attribute((u, v), 'is_ind', True if uv_i[(u, v)] in ind else False)

    if printout:
        print('Form diagram has {0} (RREF) or {1} (SVD) independent branches '.format(len(ind), m - len(s)))
        print('Shape Equilibrium Matrix: ', E.shape)
        print('Found {0} independents'.format(k))

    # Set-up

    try:
        t = form.attributes['offset']
    except:
        t = None

    lh     = normrow(C.dot(xy))**2
    Edinv  = -csr_matrix(pinv(E[:, dep]))
    Ei     = E[:, ind]
    p      = vstack([px[free], py[free]])
    q      = array([attr['q'] for u, v, attr in form.edges(True)])[:, newaxis]

    tol = 0.01
    opt_max = False
    planar = False
    tension = False
    target = False
    anchors = []
    b = None


    args   = (q, ind, dep, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, tol, z, free, fixed, planar, lh, sym, tension, k, lb, ub, lb_ind, ub_ind, opt_max, target, s, Wfree, anchors, x, y, b)
    
    return args










