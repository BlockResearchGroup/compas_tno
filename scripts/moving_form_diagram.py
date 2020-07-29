
import compas_tno

from compas_tno.solvers.WIP_formforce_solver import optimise_tna
from compas_tno.solvers.WIP_formforce_solver import optimise_fdm

from compas_tno.shapes import Shape
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.algorithms import paralelise_form

from numpy import array

from numpy import array
from numpy import max
from numpy import min
from numpy import newaxis
from numpy import zeros

from copy import deepcopy

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
from compas.numerical import normrow
from compas.numerical import connectivity_matrix

from compas_tno.plotters import plot_grad
from scipy.linalg import cho_solve, cho_factor
from numpy.linalg import matrix_rank


def tributary_and_target(XY, form, shape):
    form_ = form.copy()
    zt = array(shape.get_middle_pattern(XY)).reshape(-1, 1)
    pz = zt.copy()
    i = 0
    for key in form_.vertices():
        form_.vertex_attribute(key, 'z', value=zt[i])
        form_.vertex_attribute(key, 'x', value=XY[i][0])
        form_.vertex_attribute(key, 'y', value=XY[i][1])
        i += 1
    i = 0
    for key in form_.vertices():
        pz[i] = form_.vertex_area(key)
        i += 1

    return zt, pz


def create_load_and_target_updater(XY, form, shape):
    total_selfweight = shape.compute_selfweight()
    form_ = form.copy()

    def updater(XY):
        zt, pz = tributary_and_target(XY, form_, shape)
        pz *= total_selfweight / sum(pz)
        return zt, pz

    return updater


def grad_target(q, z, args):

    C, Ci, Cit, x, y, Wfree, s, free = args

    D = Cit.dot(diags(q.flatten())).dot(C)
    Di = D[:, free]
    D_cho = cho_factor(Di.todense())
    print(Di.shape)
    print(matrix_rank(Di.todense()))

    # Grad Calculation
    W = diags(C.dot(z).flatten())  # dimension m x m
    J = spsolve(Di, Cit.dot(W))  # dimension ni x m
    grad = ((Wfree * - 2 * (z[free] - s[free])).transpose() * J).transpose()

    return grad


def f_target(z, s):
    return sum((z - s)**2)

def plot_form_index(form):

    from compas_plotters import MeshPlotter
    plotter = MeshPlotter(form, figsize=(8, 8))
    uv_ix = form.uv_index()
    plotter.draw_edges(text={(u, v): str(uv_ix[(u,v)]) for u, v in form.edges_where({'_is_edge': True})})
    plotter.draw_vertices(radius=0.05)
    plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})], radius=0.10, facecolor='000000')
    plotter.show()

    return


# ----------------------------------
#  Main
# ----------------------------------
thk = 0.50  # thickness on the start in meters

# Basic parameters

# type_structure = 'crossvault'
# type_formdiagram = 'cross_fd'

type_structure = 'dome_polar'
type_formdiagram = 'radial_fd'

if type_structure == 'crossvault':
    discretisation = 10
    discretisation_shape = discretisation * 4
    steplength = 100
else:
    discretisation = [8, 20]
    discretisation_shape = [discretisation[0] * 4, discretisation[1] * 4]
    steplength = None

R = 5.0
exc = 0.00

# ----------------------- Form Diagram -----------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'center': [5.0, 5.0],
    'radius': R,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)

# ----------------------- Shape -----------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'center': [5.0, 5.0],
    'radius': R,
    'xy_span': [[-exc + 0.0, 10.0 + exc], [-exc + 0.0, 10.0 + exc]],
    't': 1.0
}

shape = Shape.from_library(data_shape)
swt = shape.compute_selfweight()
print('Selfweight computed:', swt)

form.selfweight_from_shape(shape)  # Maybe add this to the create_load_and_target

XY0 = form.vertices_attributes('xy')
updater = create_load_and_target_updater(XY0, form, shape)
zt, pz = updater(XY0)

objective = 'target'
form = form.initialise_tna(remove_feet=False)
internal_nodes = [not form.vertex_attribute(key, '_is_external') for key in form.vertices()]
print(internal_nodes)
print(len(internal_nodes))

# INITIALISATION
# --------------------------------------------------------------------
k_i = form.key_index()
uv_i = form.uv_index()
vertices = [k_i[key] for key in form.vertices_where({'_is_external': False})]
edges = [[k_i[u], k_i[v]] for u, v in form.edges_where({'_is_edge': True, '_is_external': False})]
fixed = [k_i[key] for key in form.vertices_where({'is_fixed': True, '_is_external': False})]
free = list(set(vertices) - set(fixed))
n = len(vertices)
print(n, 'vertices')
xyz = zeros((n, 3))
x = zeros((n, 1))
y = zeros((n, 1))
z = zeros((n, 1))
s = zeros((n, 1))
px = zeros((n, 1))
py = zeros((n, 1))
pz = zeros((n, 1))
w = zeros((n, 1))
form.update_default_vertex_attributes({'weight': 1.0})
for i, key in enumerate(form.vertices_where({'_is_external': False})):
    xyz[i, :] = form.vertex_coordinates(key)
    x[i] = form.vertex_attribute(key, 'x')
    y[i] = form.vertex_attribute(key, 'y')
    z[i] = form.vertex_attribute(key, 'z')
    s[i] = form.vertex_attribute(key, 'target')
    px[i] = form.vertex_attribute(key, 'px')
    py[i] = form.vertex_attribute(key, 'py')
    pz[i] = form.vertex_attribute(key, 'pz')
    w[i] = form.vertex_attribute(key, 'weight')
q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True, '_is_external': False})])[:, newaxis]
Wfree = diags(w[free].flatten())
C = connectivity_matrix(edges, 'csr')
print(edges)
print(len(edges))
print(C.shape)
Ci = C[:, free]
Cf = C[:, fixed]
Cit = Ci.transpose()
print(Ci.shape)
print(Cit.shape)
print(q.shape)

# Vertical Update - This has to go into Initialise TNA
zfixed = z[fixed]
Q = diags([q.ravel()], [0])
A = Cit.dot(Q).dot(Ci)
B = Cit.dot(Q).dot(Cf)
z[free] = spsolve(A, pz[free] - B.dot(zfixed)).reshape(-1, 1)

if objective == 'target':
    grad_exp = grad_target
    f_exp = f_target
    print('Optimisation takes as objective minimise distance to target')
# if objective == 'loadpath':
#     grad_exp = grad_lp
#     f_exp = loadpath
#     print('Optimisation takes as objective minimise loadpath')

args = (C, Ci, Cit, x, y, Wfree, s, free)

f = f_exp(z, s)
print('Objective Initial is f: {0}'.format(f))

# print(z)
# print(s)

it = 0
step = 100
print_ = 0
form.update_default_edge_attributes({'dq': 0.0})
# --------------------------------------------------------------------

plot = True
it_max = 1000
alpha = 99.0
a_max = 2.5
null_q = 0.20
save_steps = True

print(form.number_of_edges())

while it < it_max:

    print('\n\n')
    print('-'*20)
    print('Iteration {0}'.format(it))

    grad = grad_exp(q, z, args)

    print('grad range   : {0:.3f} : {1:.3f}'.format(min(grad), max(grad)))
    print('q start iteration range   : {0:.3f} : {1:.3f}'.format(min(q), max(q)))

    if steplength is not None:
        step = steplength
    else:
        for i in range(len(q)):
            if q[i] < 0:
                print('FD something in tension')
            t = q[i]/grad[i]
            if t > 0:
                if t < step:
                    step = t
    step = max([step, 0.1])
    incr = step/10

    print('Iterations start with steplength: {0}'.format(step))

    k = 0
    kmax = 200
    step_min = 10e-4
    q0 = deepcopy(q)
    f0 = f_exp(z, s).item()

    print('Value of Energy - start: {0:.1f}'.format(f0))
    form, force_ = form.reciprocal_from_form(remove_feet=False)
    form_ = form.copy()
    cont = True
    a_limit = a_max

    print('-'*20)

    # Define steplength

    while cont:

        q = q0 - step * grad
        k += 1

        q[q < null_q] = 0.0

        # wrong indexes comparing with the vector ...

        for index, uv in enumerate(form_.edges_where({'_is_edge': True, '_is_external': False})):
            form_.edge_attribute(uv, 'q', q[index].item())

        # if plot and k == 1:
        #     for u, v in form_.edges_where({'_is_edge': True}):
        #         [dq] = (q[uv_i[(u, v)]] - q0[uv_i[(u, v)]])
        #         form_.edge_attribute((u, v), 'dq', value=dq)

        #     plot_grad(form_).show()

        # Update Horizontal
        form_, force_ = form_.reciprocal_from_form(remove_feet=False, alpha=alpha)  # Old paralelise function
        XY = array(form_.vertices_attributes('xy'))[internal_nodes]
        diff_diag = sum(sum((array(XY) - array(XY0))**2))
        if diff_diag > 0.0:
            s, pz = updater(XY)
            XY0 = XY
            # print('UPDATE DIAG XY:', diff_diag)

        qout = array([form_.edge_attribute((u, v), 'q') for u, v in form_.edges_where({'_is_edge': True, '_is_external': False})])[:, newaxis]
        # print('q range after grad  : {0:.3f} : {1:.3f}'.format(min(qout), max(qout)))

        # Update Vertical
        Q = diags([qout.ravel()], [0])
        A = Cit.dot(Q).dot(Ci)
        B = Cit.dot(Q).dot(Cf)
        z[free] = spsolve(A, pz[free] - B.dot(zfixed)).reshape(-1, 1)

        # a = form_.evaluate_a(plot=False)

        # print(s)
        f = f_exp(z, s).item()
        # print('f = {0:.2f}'.format(f))

        # if a > a_limit:
        #     f += 10*a
        #     print('f penalized = {0:.2f}'.format(f))

        if f < f0:
            break
        else:
            if step/10 < incr:
                incr = incr/2

            step = step - incr

            if k % 10 == 0:
                print('Iteration steplength: {0} / Stepsize: {1} / Energy Var: {2:.1f}%'.format(k, step, (f-f0)/f0*100))  # / Angle Deviation: {2}
                print(f0, f)

            if k > kmax or step < step_min:
                print('Iterations reached the maximum or steplength min: It {0} / Steplength: {1}'.format(kmax, step))
                break

    # Verify Deviations

    # a = form_.evaluate_a(plot=False)

    # for u, v in form_.edges_where({'_is_edge': True}):
    #     [dq] = (q[uv_i[(u, v)]] - q0[uv_i[(u, v)]])
    #     form_.edge_attribute((u, v), 'dq', value=dq)

    # for key, vertex in form_.vertex.items():
    #     i = k_i[key]
    #     z[i] = vertex.get('z')

    # Update Function

    if f < f0:
        for i, key in enumerate(form.vertices_where({'_is_external': False})):
            form_.vertex_attribute(key, 'z', z[i].item())
        form = form_
        q = array([form.edge_attribute((u, v), 'q') for u, v in form.edges_where({'_is_edge': True, '_is_external': False})])[:, newaxis]
        # q = qout  # simpler?
        steplength = step
        # incrsave = incr
        print('End of Iteration {0} / Step number: {1} / Stepsize: {2:.5f} / Evergy Var: {3:.1f}% / Energy: {4:.1f}'.format(it, k, step, (f-f0)/f0*100, f))  # / Angle Dev: {5:.3f}
        # form = paralelise_form(form,q)

    else:
        print('Optimisation not found - Iteration {0} / Stepsize: {1:.5f} / Energy: {2:.1f}'.format(it, step, f))  # / Angle Dev: {3}
        plot_form(form).show()
        # plot_force(force, form).show()
        break

    if save_steps and it % 10 == 0:
        print_ += 1
        form.to_json('/Users/mricardo/compas_dev/me/bestfit/dome/dome_' + str(it) + '_alpha=' + str(alpha) + '.json')

    it += 1

    # if plot:
    #     if it % 5 == 0:
    #         plot_form(form).show()
    #         # plot_force(force, form).show()


form.to_json('/Users/mricardo/compas_dev/me/bestfit/dome/dome_'+type_formdiagram+'_final_alpha='+str(alpha)+'.json')
