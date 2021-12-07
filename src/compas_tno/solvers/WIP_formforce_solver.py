
from numpy import array
# from numpy import max
# from numpy import min
from numpy import newaxis
from numpy import zeros

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags

from compas.numerical import connectivity_matrix

from compas_tno.algorithms import equilibrium_fdm
from compas_tno.algorithms import vertical_equilibrium_fdm

# from compas_tno.diagrams.form import evaluate_a
# from compas_tno.diagrams.form import energy
# from compas_tno.diagrams.form import loadpath

# from compas_tna.diagrams import ForceDiagram

from compas_tno.plotters import plot_form
# from compas_tno.plotters import plot_force
# from compas_tno.plotters import plot_grad

from copy import deepcopy


__all__ = [
    'optimise_fdm',
    'optimise_tna',
]


def optimise_fdm(form, plot=True, surf=False):

    # Mapping

    k_i = form.key_index()

    # Vertices and edges

    vertices = [k_i[key] for key in form.vertices_where({'_is_external': False})]
    edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True, '_is_external': False})]
    fixed = [k_i[key] for key in form.vertices_where({'_is_external': False, 'is_fixed': True})]
    free = list(set(vertices) - set(fixed))

    n = len(vertices)

    # Co-ordinates and loads

    xyz = zeros((n, 3))
    x = zeros((n, 1))
    y = zeros((n, 1))
    z = zeros((n, 1))
    s = zeros((n, 1))
    px = zeros((n, 1))
    py = zeros((n, 1))
    pz = zeros((n, 1))

    for key, vertex in form.vertex.items():
        if vertex.get('_is_external') is False:
            i = k_i[key]
            xyz[i, :] = form.vertex_coordinates(key)
            x[i] = vertex.get('x')
            y[i] = vertex.get('y')
            z[i] = vertex.get('z')
            s[i] = vertex.get('target')
            px[i] = vertex.get('px', 0)
            py[i] = vertex.get('py', 0)
            pz[i] = vertex.get('pz', 0)

    q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True, '_is_external': False})])[:, newaxis]

    # Interpolate Target Surface

    # if interpolate:
    #     from scipy import interpolate
    #     surf = interpolate.interp2d(x, y, s, kind = 'linear')

    # C and E matrices

    C = connectivity_matrix(edges, 'csr')
    Ci = C[:, free]
    Cit = Ci.transpose()

    # From Here calculate Iteractivelly q...

    it = 0
    it_max = 100
    stepsave = 10.0
    incrsave = 0.1

    while it < it_max:

        print('-'*10, ' Iteration {0}'.format(it))

        D = Cit.dot(diags(q.flatten())).dot(C)
        Di = D[:, free]

        # Grad Calculation

        # (L, values) = _chofactor(Di.toarray())

        W = diags(C.dot(z).flatten())  # dimension m x m
        J = spsolve(Di, Cit.dot(W))  # dimension ni x m

        # Y_ = spsolve(L,Cit.dot(W))
        # X_ = spsolve(L.transpose(),Y_)

        # Gradient of F

        # Grad = 2 * (diags(z[free,0].flatten()) - diags(s[free,0].flatten())).dot(X_)

        grad = (2 * (z[free] - s[free]).transpose() * J).transpose()

        step = stepsave
        incr = incrsave

        k = 0
        kmax = 300
        q0 = q
        f0 = form.loadpath()
        form_ = deepcopy(form)
        cont = True

        # Define steplength

        while cont:

            q = q0 + step * grad
            k += 1

            # form_ = update_form(form_, q)
            vertical_equilibrium_fdm(form_)
            f = form.loadpath()

            if f < f0:
                break
            else:
                if step/10 < incr:
                    incr = incr/2

                step = step - incr

                if k % 100 == 0:
                    print('Iteration steplength: {0} / Stepsize: {1} / Energy Var: {2:.1f}%'.format(k, step, (f-f0)/f0*100))

                if k > kmax:
                    print('Iterations on steplength reached the maximum {0}'.format(kmax))
                    break

        # Update Function

        if f < f0:
            form = form_
            stepsave = step
            incrsave = incr
            q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]
            for key, vertex in form_.vertex.items():
                i = k_i[key]
                z[i] = vertex.get('z')
            print('End of Iteration {0} / Stepsize: {1:.5f} / Evergy Var: {2:.1f}% / Energy: {3:.1f}'.format(it, step, (f-f0)/f0*100, f))

            form = equilibrium_fdm(form)

            if surf:
                for key, vertex in form.vertex.items():
                    xi, yi, _ = form.vertex_coordinates(key)
                    [si] = surf(xi, yi).tolist()
                    form.vertex_attribute(key, 'target', value=si)
                print('Updating constraints. Energy: {0:.1f}'.format(form.loadpath()))

        else:
            print('Optimisation not found - Iteration {0} / Stepsize: {1:.5f} / Energy: {2:.1f}'.format(it, step, f))
            break

        it += 1

        if plot:
            plot_form(form).show()

    return form


def optimise_tna(form, objective='target', plot=None, it_max=20, alpha=1.0, a_max=2.5, steplength=None, null=None, save_steps=False):

    # # Mapping

    # k_i = form.key_index()
    # uv_i = form.uv_index()

    # # Vertices and edges

    # vertices = [k_i[key] for key in form.vertices()]
    # edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True})]
    # fixed = [k_i[key] for key in form.vertices_where({'is_anchor': True})]
    # free = list(set(vertices) - set(fixed))
    # if null:
    #     null_i = [uv_i[(u, v)] for u, v in null]

    # n = len(vertices)

    # # Co-ordinates and loads

    # xyz = zeros((n, 3))
    # x = zeros((n, 1))
    # y = zeros((n, 1))
    # z = zeros((n, 1))
    # s = zeros((n, 1))
    # px = zeros((n, 1))
    # py = zeros((n, 1))
    # pz = zeros((n, 1))
    # w = zeros((n, 1))

    # for key, vertex in form.vertex.items():
    #     i = k_i[key]
    #     xyz[i, :] = form.vertex_coordinates(key)
    #     x[i] = vertex.get('x')
    #     y[i] = vertex.get('y')
    #     z[i] = vertex.get('z')
    #     s[i] = vertex.get('target', 0)
    #     px[i] = vertex.get('px', 0)
    #     py[i] = vertex.get('py', 0)
    #     pz[i] = vertex.get('pz', 0)
    #     w[i] = vertex.get('weight', 1.0)

    # q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]
    # Wfree = diags(w[free].flatten())

    # C = connectivity_matrix(edges, 'csr')
    # Ci = C[:, free]
    # Cit = Ci.transpose()

    # if objective == 'target':
    #     grad_exp = grad_target
    #     f_exp = form.distance_target()
    #     print('Optimisation takes as objective minimise distance to target')
    # if objective == 'loadpath':
    #     grad_exp = grad_lp
    #     f_exp = form.loadpath()
    #     print('Optimisation takes as objective minimise loadpath')

    # args = (C, Ci, Cit, x, y, Wfree, s, free)

    # f = f_exp(form)
    # print('Objective Initial is f: {0}'.format(f))

    # it = 0
    # step = 100
    # print_ = 0
    # form.update_default_edge_attributes({'dq': 0.0})

    # while it < it_max:

    #     print('\n', '-'*20, '\n', ' Iteration {0}'.format(it))

    #     grad = grad_exp(q, z, args)

    #     print('grad range   : {0:.3f} : {1:.3f}'.format(min(grad), max(grad)))
    #     print('q start iteration range   : {0:.3f} : {1:.3f}'.format(min(q), max(q)))

    #     if steplength is not None:
    #         step = steplength
    #     else:
    #         for i in range(len(q)):
    #             if q[i] < 0:
    #                 print('FD something in tension')
    #             [t] = q[i]/grad[i]
    #             if t > 0:
    #                 if t < step:
    #                     step = t
    #     step = max([step, 0.1])
    #     incr = step/10

    #     print('Iterations start with steplength: {0}'.format(step))

    #     k = 0
    #     kmax = 200
    #     q0 = deepcopy(q)
    #     f0 = f_exp(form)
    #     print('Value of Energy - start: {0:.1f}'.format(f0))
    #     form_ = deepcopy(form)
    #     force_ = ForceDiagram.from_formdiagram(form_)
    #     cont = True
    #     a_limit = a_max

    #     print('-'*20)

    #     # Define steplength

    #     while cont:

    #         q = q0 - step * grad
    #         k += 1

    #         if null:
    #             q[null_i] *= 0

    #         if plot and k == 1:
    #             for u, v in form_.edges_where({'_is_edge': True}):
    #                 [dq] = (q[uv_i[(u, v)]] - q0[uv_i[(u, v)]])
    #                 form_.edge_attribute((u, v), 'dq', value=dq)

    #             plot_grad(form_).show()

    #         form_, force = paralelise_form(form_, force_, q.flatten(), alpha=alpha)  # , plot=True)
    #         a = max(form.edges_attribute('a'))

    #         f = f_exp(form_)
    #         # print('f = {0:.2f}'.format(f))
    #         if a > a_limit:
    #             f += 10*a
    #             # print('f penalized = {0:.2f}'.format(f))

    #         # qout = array([form_.edge_attribute((u,v),'q') for u, v in form_.edges_where({'_is_edge': True})])[:, newaxis]
    #         # print('q range after grad  : {0:.3f} : {1:.3f}'.format(min(qout), max(qout)))

    #         # plot_form(form_).show()

    #         if f < f0:
    #             break
    #         else:
    #             if step/10 < incr:
    #                 incr = incr/2

    #             step = step - incr

    #             if k % 10 == 0:
    #                 print('Iteration steplength: {0} / Stepsize: {1} / Angle Deviation: {2} / Energy Var: {3:.1f}%'.format(k, step, a, (f-f0)/f0*100))

    #             if k > kmax:
    #                 print('Iterations on steplength reached the maximum {0}'.format(kmax))
    #                 break

    #     # Verify Deviations

    #     # a = evaluate_a(form, print=False)

    #     for u, v in form_.edges_where({'_is_edge': True}):
    #         [dq] = (q[uv_i[(u, v)]] - q0[uv_i[(u, v)]])
    #         form_.edge_attribute((u, v), 'dq', value=dq)

    #     for key, vertex in form_.vertex.items():
    #         i = k_i[key]
    #         z[i] = vertex.get('z')

    #     # Update Function

    #     if f < f0:
    #         form = form_
    #         q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True})])[:, newaxis]
    #         # stepsave = step
    #         # incrsave = incr
    #         print('End of Iteration {0} / Step number: {1} / Stepsize: {2:.5f} / Evergy Var: {3:.1f}% / Energy: {4:.1f} / Angle Dev: {5:.3f}'.format(it, k, step,
    #                                                                                                                                                  (f-f0)/f0*100, f, a))
    #         # form = paralelise_form(form,q)

    #     else:
    #         print('Optimisation not found - Iteration {0} / Stepsize: {1:.5f} / Energy: {2:.1f} / Angle Dev: {3}'.format(it, step, f, a))
    #         plot_form(form_).show()
    #         plot_force(force, form).show()
    #         break

    #     if save_steps and it % 10 == 0:
    #         print_ += 1
    #         form.to_json('/Users/mricardo/compas_dev/me/bestfit/iterations/pillow_' + str(print_) + '.json')

    #     it += 1

    #     if plot:
    #         if it % 5 == 0:
    #             plot_form(form).show()
    #             plot_force(force, form).show()

    print('WIP')

    return form


def grad_target(q, z, args):

    C, Ci, Cit, x, y, Wfree, s, free = args

    D = Cit.dot(diags(q.flatten())).dot(C)
    Di = D[:, free]

    # Grad Calculation

    W = diags(C.dot(z).flatten())  # dimension m x m
    J = spsolve(Di, Cit.dot(W))  # dimension ni x m
    grad = ((Wfree * - 2 * (z[free] - s[free])).transpose() * J).transpose()

    return grad


def grad_lp(q, z, args):

    C, Ci, Cit, x, y, Wfree, s, free = args

    D = Cit.dot(diags(q.flatten())).dot(C)
    Di = D[:, free]

    # Grad Calculation

    grad1 = (C.dot(x))**2 + (C.dot(y))**2 + (C.dot(z))**2
    W = diags(C.dot(z).flatten())  # dimension m x m
    J = spsolve(Di, Cit.dot(W))  # dimension ni x m
    M = Ci.dot(J)  # dimension m x m verify
    K = W.dot(M)  # Shape m x 1
    grad2 = (2 * q.transpose() * K).transpose()
    grad = grad1 + grad2

    return grad
