
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

from scipy.sparse.linalg import spsolve
from scipy.optimize import fmin_slsqp
from scipy.sparse import csr_matrix
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import factorized
from compas.numerical import normrow
from compas.numerical import normalizerow

from compas.numerical import connectivity_matrix
from compas.numerical import devo_numpy
from compas.numerical import equilibrium_matrix
from compas.numerical import normrow
from compas.numerical import nonpivots
from compas.numerical.linalg import _chofactor
from compas_plotters import MeshPlotter
from compas.utilities import geometric_key

from compas_thrust.algorithms.equilibrium import z_from_form
from compas_thrust.algorithms.equilibrium import update_form
from compas_thrust.algorithms.equilibrium import paralelise_form

from compas_thrust.diagrams.form import evaluate_a

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram

from compas_thrust.plotters.plotters import plot_form
from compas_thrust.plotters.plotters import plot_force
from compas_thrust.plotters.plotters import plot_grad

from copy import deepcopy

__author__    = ['Ricardo Maia Avelino <mricardo@ethz.ch>']
__copyright__ = 'Copyright 2019, BLOCK Research Group - ETH Zurich'
__license__   = 'MIT License'
__email__     = 'mricardo@ethz.ch'


__all__ = [
    'lagrangian_scale',
    'evaluate_scale',
    'optimise_fdm',
    'optimise_tna',
]

def lagrangian_scale(form):

    # Mapping

    k_i  = form.key_index()
    uv_i = form.uv_index()

    # Vertices and edges

    vertices = [k_i[key] for key in form.vertices_where({'is_external': False})]
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True, 'is_external' : False})]
    fixed = [k_i[key] for key in form.vertices_where({'is_external': False, 'is_fixed': True})]
    free  = list(set(vertices) - set(fixed))

    m = len(edges)
    n = len(vertices)
    ni = len(free)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))
    q   = zeros((m, 1))

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            i = k_i[key]
            xyz[i, :] = form.vertex_coordinates(key)
            x[i]  = vertex.get('x')
            y[i]  = vertex.get('y')
            z[i]  = vertex.get('z')
            s[i]  = vertex.get('target')
            px[i] = vertex.get('px', 0)
            py[i] = vertex.get('py', 0)
            pz[i] = vertex.get('pz', 0)

    q   = zeros((m, 1))

    for u,v in form.edges():
        if form.get_edge_attribute((u,v), 'is_edge') == True:
            i = uv_i[(u,v)]
            q[i] = form.get_edge_attribute((u,v), 'q')

    p = pz[free]

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cf  = C[:, fixed]
    Cit = Ci.transpose()
    D = Cit.dot(diags(q.flatten())).dot(C)
    Di = D[:, free]
    Df = D[:, fixed]
    Dt = D.transpose()

    args = q, free, fixed, C, Ci, Cit, Cf, pz, s, z

    A11 = diags([1]*n+[0]).toarray()
    A12 = vstack([Dt.toarray(),-1*p.transpose()])
    A21 = hstack([D.toarray(),-1*p])
    A22 = zeros((ni,ni))
    A = vstack([hstack([A11,A12]),hstack([A21,A22])])
    B = vstack([1*s,zeros((1+ni, 1))])

    sol = spsolve(A, B)
    r = sol[n]
    # r = 0.59
    print('The lagragian scale is: {0:.3f}'.format(r))

    # r = 2.0

    # Update zs

    q = q * r
    Q = diags([q.flatten()], [0])
    D = Cit.dot(Q).dot(C)
    Di = D[:, free]
    Df = D[:, fixed]
    A_= Di
    b_ = pz[free, 0] - Df.dot(z[fixed,0])
    z[free, 0] = spsolve(A_, b_)

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            i = k_i[key]
            [zi] = z[i]
            form.set_vertex_attribute(key, 'z', value = zi)

    for u, v in form.edges_where({'is_external': False}):
        if form.get_edge_attribute((u,v),'is_edge') is True:
            i = uv_i[(u,v)]
            [qi] = q[i]
            form.set_edge_attribute((u,v),'q',value=qi)

    print('Final energy {0}'.format(energy(form)))

    return form

def energy(form):

    f = 0
    for key, vertex in form.vertex.items():
            if vertex.get('is_external') == False:
                z  = vertex.get('z')
                s  = vertex.get('target')
                w = vertex.get('weight', 1.0)
                f += w * (z - s)**2

    return f

def loadpath(form):

    lp = 0
    for u, v in form.edges_where({'is_external': False}):
            if form.get_edge_attribute((u,v),'is_edge') is True and form.get_edge_attribute((u,v),'is_symmetry') is False:
                qi = form.get_edge_attribute((u,v),'q')
                li = form.edge_length(u,v)
                lp += qi*li**2

    return lp

def evaluate_scale(form, function, bounds, n = 100, plot = True):

    r0 = bounds[0]
    stp = (bounds[1]-bounds[0])/n
    x = []
    y = []
    q0 = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]

    form_ = deepcopy(form)

    k_i  = form_.key_index()
    uv_i = form_.uv_index()

    for k in range(n):
        r = r0 + stp * k
        q = q0 * r
        x.append(r)

        for u, v in form_.edges_where({'is_edge': True}):
            i = uv_i[(u,v)]
            [qi] = q[i]
            form_.set_edge_attribute((u,v),'q',value=qi)

        form_ = z_from_form(form_)
        y.append(function(form_))

    import matplotlib
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(x, y)

    ax.set(xlabel='Scale (r)', ylabel='Energy (f)',
        title='Evaluate Energy by Scaling')
    ax.grid()

    pos = argmin(y)
    xmin = x[pos]
    ymin = y[pos]
    text= "x={:.3f}, y={:.3f}".format(xmin, ymin)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
            arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmin, ymin), xytext=(0.94,0.96), **kw)

    if plot:
        plt.show()

    return xmin

def optimise_fdm(form, plot= True, surf=False):

    # Mapping

    k_i  = form.key_index()
    i_k  = form.index_key()
    i_uv = form.index_uv()
    uv_i = form.uv_index()

    # Vertices and edges

    vertices = [k_i[key] for key in form.vertices_where({'is_external': False})]
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True, 'is_external' : False})]
    fixed = [k_i[key] for key in form.vertices_where({'is_external': False, 'is_fixed': True})]
    free  = list(set(vertices) - set(fixed))

    m = len(edges)
    n = len(vertices)
    ni = len(free)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))

    for key, vertex in form.vertex.items():
        if vertex.get('is_external') == False:
            i = k_i[key]
            xyz[i, :] = form.vertex_coordinates(key)
            x[i]  = vertex.get('x')
            y[i]  = vertex.get('y')
            z[i]  = vertex.get('z')
            s[i]  = vertex.get('target')
            px[i] = vertex.get('px', 0)
            py[i] = vertex.get('py', 0)
            pz[i] = vertex.get('pz', 0)

    q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True, 'is_external' : False})])[:, newaxis]
    p = pz[free]

    # Interpolate Target Surface

    # if interpolate:
    #     from scipy import interpolate
    #     surf = interpolate.interp2d(x, y, s, kind = 'linear')

    # C and E matrices

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cit = Ci.transpose()

    # From Here calculate Iteractivelly q...

    it = 0
    it_max = 100

    while it < it_max:

        print('-'*10,' Iteration {0}'.format(it))

        D = Cit.dot(diags(q.flatten())).dot(C)
        Di = D[:, free]

        # Grad Calculation

        # (L, values) = _chofactor(Di.toarray())

        W = diags(C.dot(z).flatten())  # dimension m x m
        J = spsolve(Di, Cit.dot(W)) # dimension ni x m

        # Y_ = spsolve(L,Cit.dot(W))
        # X_ = spsolve(L.transpose(),Y_)

        # Gradient of F

        # Grad = 2 * (diags(z[free,0].flatten()) - diags(s[free,0].flatten())).dot(X_)

        grad = (2 * (z[free] - s[free]).transpose() * J).transpose()

        try:
            step = stepsave
            incr = incrsave
        except:
            step = 10.0
            incr = 0.1

        k = 0
        kmax = 300
        q0 = q
        f0 = energy(form)
        form_ = deepcopy(form)
        cont = True

        # Define steplength

        while cont:

            q = q0 + step * grad
            k += 1

            form_ = update_form(form_,q)
            f = energy(form_)

            if f < f0:
                break
            else:
                if step/10 <incr:
                    incr = incr/2

                step = step - incr

                if k%100 == 0:
                    print('Iteration steplength: {0} / Stepsize: {1} / Energy Var: {2:.1f}%'.format(k, step, (f-f0)/f0*100))

                if k>kmax:
                    print('Iterations on steplength reached the maximum {0}'.format(kmax))
                    break

        # Update Function

        if f < f0:
            form = form_
            stepsave = step
            incrsave = incr
            q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]
            for key, vertex in form_.vertex.items():
                i = k_i[key]
                z[i]  = vertex.get('z')
            print('End of Iteration {0} / Stepsize: {1:.5f} / Evergy Var: {2:.1f}% / Energy: {3:.1f}'.format(it, step, (f-f0)/f0*100, f))

            form = z_from_form(form)

            if surf:
                for key, vertex in form.vertex.items():
                    xi, yi, _ = form.vertex_coordinates(key)
                    [si] = surf(xi,yi).tolist()
                    form.set_vertex_attribute(key,'target',value=si)
                print('Updating constraints. Energy: {0:.1f}'.format(energy(form)))

        else:
            print('Optimisation not found - Iteration {0} / Stepsize: {1:.5f} / Energy: {2:.1f}'.format(it, step, f))
            break

        it += 1

        if plot:
            plot_form(form).show()


    return form

def optimise_tna(form, objective ='target', plot= None, it_max=20, alpha=1.0, a_max = 2.5, steplength=None, null=None, save_steps=False):

    # Mapping

    k_i  = form.key_index()
    uv_i = form.uv_index()

    # Vertices and edges

    vertices = [k_i[key] for key in form.vertices()]
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True})]
    fixed = [k_i[key] for key in form.vertices_where({'is_anchor': True})]
    free  = list(set(vertices) - set(fixed))
    if null:
        null_i = [uv_i[(u,v)] for u, v in null]

    n = len(vertices)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))
    w = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i]  = vertex.get('x')
        y[i]  = vertex.get('y')
        z[i]  = vertex.get('z')
        s[i]  = vertex.get('target', 0)
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)
        w[i] = vertex.get('weight', 1.0)

    q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]
    Wfree = diags(w[free].flatten())

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cit = Ci.transpose()

    f = energy(form)
    print('Energy is f: {0}'.format(f))

    it = 0
    step = 100
    print_ = 0
    form.update_default_edge_attributes({'dq': 0.0})

    while it < it_max:

        print('\n','-'*20,'\n',' Iteration {0}'.format(it))

        D = Cit.dot(diags(q.flatten())).dot(C)
        Di = D[:, free]

        # Grad Calculation

        W = diags(C.dot(z).flatten())  # dimension m x m
        J = spsolve(Di, Cit.dot(W)) # dimension ni x m
        grad = ( ( Wfree * - 2 * (z[free] - s[free]) ).transpose() * J).transpose()
        print(grad.shape)

        print('grad range   : {0:.3f} : {1:.3f}'.format(min(grad), max(grad)))
        print('q start iteration range   : {0:.3f} : {1:.3f}'.format(min(q), max(q)))

        if steplength is not None:
            step = steplength
        else:
            for i in range(len(q)):
                if q[i] < 0:
                    print('FD something in tension')
                [t] = q[i]/grad[i]
                if t > 0:
                    if t < step:
                        step = t
        step = max([step, 0.1])
        incr = step/10

        print('Iterations start with steplength: {0}'.format(step))

        k = 0
        kmax = 200
        q0 = deepcopy(q)
        f0 = energy(form)
        print('Value of Energy - start: {0:.1f}'.format(f0))
        form_ = deepcopy(form)
        force_ = ForceDiagram.from_formdiagram(form_)
        cont = True
        a_limit = 2.5

        print('-'*20)

        # Define steplength

        while cont:

            q = q0 - step * grad
            k += 1

            if null:
                q[null_i] *= 0

            if plot and k == 1:
                for u, v in form_.edges_where({'is_edge': True}):
                    [dq] = (q[uv_i[(u,v)]] - q0[uv_i[(u,v)]])
                    form_.set_edge_attribute((u,v),'dq',value = dq)

                plot_grad(form_).show()

            form_, force = paralelise_form(form_,force_, q.flatten(), alpha=alpha) #, plot=True)
            a = evaluate_a(form_, plot=False)

            f = energy(form_)
            # print('f = {0:.2f}'.format(f))
            if a > a_limit:
                f += 10*a
                # print('f penalized = {0:.2f}'.format(f))

            # qout = array([form_.get_edge_attribute((u,v),'q') for u, v in form_.edges_where({'is_edge': True})])[:, newaxis]
            # print('q range after grad  : {0:.3f} : {1:.3f}'.format(min(qout), max(qout)))

            # plot_form(form_).show()

            if f < f0:
                break
            else:
                if step/10 <incr:
                    incr = incr/2

                step = step - incr

                if k%10 == 0:
                    print('Iteration steplength: {0} / Stepsize: {1} / Angle Deviation: {2} / Energy Var: {3:.1f}%'.format(k, step, a, (f-f0)/f0*100))

                if k>kmax:
                    print('Iterations on steplength reached the maximum {0}'.format(kmax))
                    break

        # Verify Deviations

        # a = evaluate_a(form, print=False)

        for u, v in form_.edges_where({'is_edge': True}):
                [dq] = (q[uv_i[(u,v)]] - q0[uv_i[(u,v)]])
                form_.set_edge_attribute((u,v),'dq',value = dq)

        for key, vertex in form_.vertex.items():
            i = k_i[key]
            z[i]  = vertex.get('z')


        # Update Function

        if f < f0:
            form = form_
            q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]
            # stepsave = step
            # incrsave = incr
            print('End of Iteration {0} / Step number: {1} / Stepsize: {2:.5f} / Evergy Var: {3:.1f}% / Energy: {4:.1f} / Angle Dev: {5:.3f}'.format(it, k, step, (f-f0)/f0*100, f, a))
            # form = paralelise_form(form,q)

        else:
            print('Optimisation not found - Iteration {0} / Stepsize: {1:.5f} / Energy: {2:.1f} / Angle Dev: {3}'.format(it, step, f, a))
            plot_form(form_).show()
            plot_force(force).show()
            break

        if save_steps and it%10 == 0:
            print_ +=1
            form.to_json('/Users/mricardo/compas_dev/me/bestfit/iterations/pillow_' + str(print_) + '.json')

        it += 1

        if plot:
            if it%5 == 0:
                plot_form(form).show()
                plot_force(force).show()

            # print('grad is')
            # print(grad)
            # print('q is')
            # print(q[free])
            # form_ = paralelise_form(form_, q0, plot=True)

    return form

def optimise_loadpath(form, plot= None, it_max=20, alpha=1.0, a_max = 2.5, steplength=None):

    # Mapping

    k_i  = form.key_index()
    uv_i = form.uv_index()

    # Vertices and edges

    vertices = [k_i[key] for key in form.vertices()]
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'is_edge': True})]
    fixed = [k_i[key] for key in form.vertices_where({'is_anchor': True})]
    free  = list(set(vertices) - set(fixed))

    n = len(vertices)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x   = zeros((n, 1))
    y   = zeros((n, 1))
    z   = zeros((n, 1))
    s   = zeros((n, 1))
    px  = zeros((n, 1))
    py  = zeros((n, 1))
    pz  = zeros((n, 1))

    for key, vertex in form.vertex.items():
        i = k_i[key]
        xyz[i, :] = form.vertex_coordinates(key)
        x[i]  = vertex.get('x')
        y[i]  = vertex.get('y')
        z[i]  = vertex.get('z')
        s[i]  = vertex.get('target', 0)
        px[i] = vertex.get('px', 0)
        py[i] = vertex.get('py', 0)
        pz[i] = vertex.get('pz', 0)

    q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]

    C   = connectivity_matrix(edges, 'csr')
    Ci  = C[:, free]
    Cit = Ci.transpose()

    f = loadpath(form)
    print('Loadpath is f: {0}'.format(f))

    it = 0
    step = 100
    form.update_default_edge_attributes({'dq': 0.0})

    while it < it_max:

        print('\n','-'*20,'\n',' Iteration {0}'.format(it))

        D = Cit.dot(diags(q.flatten())).dot(C)
        Di = D[:, free]

        # Grad Calculation

        grad1 = (C.dot(x))**2 + (C.dot(y))**2 + (C.dot(z))**2

        W = diags(C.dot(z).flatten())  # dimension m x m
        J = spsolve(Di, Cit.dot(W)) # dimension ni x m
        M = Ci.dot(J) # dimension m x m verify
        K = W.dot(M) # Shape m x 1
        grad2 = ( 2 * q.transpose() * K ).transpose()

        print(grad2.shape)
        print(grad1.shape)

        grad = grad1 + grad2

        print(grad.shape)

        print('grad range   : {0:.3f} : {1:.3f}'.format(min(grad), max(grad)))
        print('q start iteration range   : {0:.3f} : {1:.3f}'.format(min(q), max(q)))

        if steplength is not None:
            step = steplength
        else:
            for i in range(len(q)):
                if q[i] < 0:
                    print('FD something in tension')
                [t] = q[i]/grad[i]
                if t > 0:
                    if t < step:
                        step = t
        incr = step/10

        print('Iterations start with steplength: {0}'.format(step))

        k = 0
        kmax = 200
        q0 = deepcopy(q)
        f0 = loadpath(form)
        print('Value of Loadpath - start: {0:.1f}'.format(f0))
        form_ = deepcopy(form)
        force_ = ForceDiagram.from_formdiagram(form_)
        cont = True
        a_limit = 2.5

        print('-'*20)

        # Define steplength

        while cont:

            q = q0 + step * grad
            k += 1

            form_ = z_from_form(form_)
            # print('q range step_it  : {0:.3f} : {1:.3f}'.format(min(q), max(q)))

            if plot:
                if k == 1:
                    for u, v in form_.edges_where({'is_edge': True}):
                        [dq] = (q[uv_i[(u,v)]] - q0[uv_i[(u,v)]])
                        form_.set_edge_attribute((u,v),'dq',value = dq)

                    plot_grad(form_).show()

            form_, force = paralelise_form(form_,force_, q.flatten(), alpha=alpha) #, plot=True)
            a = evaluate_a(form_, plot=False)

            f = loadpath(form_)
            # print('f = {0:.2f}'.format(f))
            if a > a_limit:
                f += 10*a
                # print('f penalized = {0:.2f}'.format(f))

            # qout = array([form_.get_edge_attribute((u,v),'q') for u, v in form_.edges_where({'is_edge': True})])[:, newaxis]
            # print('q range after grad  : {0:.3f} : {1:.3f}'.format(min(qout), max(qout)))

            # plot_form(form_).show()

            if f < f0:
                break
            else:
                if step/10 <incr:
                    incr = incr/2

                step = step - incr

                if k%10 == 0:
                    print('Iteration steplength: {0} / Stepsize: {1} / Angle Deviation: {2} / Energy Var: {3:.1f}%'.format(k, step, a, (f-f0)/f0*100))

                if k>kmax:
                    print('Iterations on steplength reached the maximum {0}'.format(kmax))
                    break

        # Verify Deviations

        # a = evaluate_a(form, print=False)

        for u, v in form_.edges_where({'is_edge': True}):
                [dq] = (q[uv_i[(u,v)]] - q0[uv_i[(u,v)]])
                form_.set_edge_attribute((u,v),'dq',value = dq)

        for key, vertex in form_.vertex.items():
            i = k_i[key]
            x[i]  = vertex.get('x')
            y[i]  = vertex.get('y')
            z[i]  = vertex.get('z')


        # Update Function

        if f < f0:
            form = form_
            q = array([form.get_edge_attribute((u,v),'q') for u, v in form.edges_where({'is_edge': True})])[:, newaxis]
            # stepsave = step
            # incrsave = incr
            print('End of Iteration {0} / Step number: {1} / Stepsize: {2:.5f} / Evergy Var: {3:.1f}% / Energy: {4:.1f} / Angle Dev: {5:.3f}'.format(it, k, step, (f-f0)/f0*100, f, a))
            # form = paralelise_form(form,q)

        else:
            print('Optimisation not found - Iteration {0} / Stepsize: {1:.5f} / Energy: {2:.1f} / Angle Dev: {3}'.format(it, step, f, a))
            plot_form(form_).show()
            plot_force(force).show()
            break

        it += 1

        if plot:
            if it%5 == 0:
                plot_form(form).show()
                plot_force(force).show()

            # print('grad is')
            # print(grad)
            # print('q is')
            # print(q[free])
            # form_ = paralelise_form(form_, q0, plot=True)

    return form


