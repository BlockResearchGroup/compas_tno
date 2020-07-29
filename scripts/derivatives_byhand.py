from scipy.sparse import diags
from numpy import hstack

from compas_tno.algorithms import compute_jacobian
from compas_tno.algorithms import compute_grad
from compas_tno.algorithms import f_objective_pytorch
from compas_tno.algorithms import f_constraints_pytorch
from torch import tensor
from compas_tno.problems.setup import set_symmetry_constraint
from compas_tno.problems.setup import set_joints_constraint
from compas_tno.problems.setup import set_rollers_constraint
from compas_tno.problems.setup import set_cracks_constraint
from compas_tno.problems.setup import set_b_constraint
import matplotlib.pyplot as plt
from numpy import zeros
from numpy import array
from compas_tno.problems.problems import initialise_problem

from compas_tno.problems.objectives import f_min_thrust
from compas_tno.problems.objectives import f_max_thrust

from compas_tno.problems.constraints import constr_wrapper

from numpy import append

from compas_plotters import MeshPlotter
from compas.utilities import i_to_red
from compas.utilities import i_to_rgb

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents

from numpy.linalg import norm
import matplotlib.pyplot as plt
from compas_tno.algorithms import zlq_from_qid
from compas_tno.algorithms import apply_sag
from compas_tno.problems import sensitivities_wrapper
from compas_tno.problems import constr_wrapper_ipopt
from compas_tno.problems import gradient_fmin
from compas.utilities import geometric_key


# import sys
# import numpy
# numpy.set_printoptions(threshold=sys.maxsize)

exitflag = 0  # means that optimisation found a solution
t0 = thk = 0.50  # thickness on the start in meters
thk_reduction = 0.05  # in meters
thk_refined = 0.01
solutions_min = []  # empty lists to keep track of the solutions for min thrust
solutions_max = []  # empty lists to keep track of the solutions for max thrust
size_parameters = []  # empty lists to keep track of  the parameters
thicknesses = []
span = 10.0  # square span for analysis
k = 1
hc = 5.0
# [5.00, 5.92, 6.71, 7.42, 8.06, 8.66, 9.22, 9.75]

type_structure = 'pointed_crossvault'
type_formdiagram = 'fan_fd'  # Try also 'fan_fd'
discretisation = 8
discretisation_shape = discretisation

# type_structure = 'dome'
# type_formdiagram = 'radial_fd'  # Try also 'fan_fd'
# discretisation = [8, 20]
# discretisation_shape = [discretisation[0]*2, discretisation[1]*2]

# ----------------------- Create Form Diagram for analysis ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
    'center': [5.0, 5.0],
    'radius': span/2,
    'r_oculus': 0.0,
    'diagonal': True,
    'partial_diagonal': 'right',
}

import compas_tno
form = FormDiagram.from_library(data_diagram)
apply_sag(form)
# path = '/Users/mricardo/compas_dev/compas_tno/data/dome/Dome_Px=0.24_discr_[8, 20]_min.json'
# path = compas_tno.get('test-min.json')
# path = compas_tno.get('test.json')
# form = FormDiagram.shuffle_diagram(form)
# plot_form(form, show_q=False).show()
# path = '/Users/mricardo/compas_dev/me/shape_comparison/pointed_crossvault/cross_fd/h=5.0/pointed_crossvault_cross_fd_discr_10_min_thk_50.0.json'
# form = FormDiagram.from_json(path)

plot_form(form).show()

# --------------------- Create Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'IPOPT'
optimiser.data['solver'] = 'IPOPT'
optimiser.data['constraints'] = ['funicular', 'envelope']  # Note addition of constraint on rollers
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['printout'] = False
optimiser.data['objective'] = 'min'
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 3000.0
optimiser.data['gradient'] = True
optimiser.data['jacobian'] = True

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'xy_span': [[0, span], [0, k*span]],
    'hc': hc,
    'hm': None,
    'he': None,
    'center': [5.0, 5.0],
    'radius': span/2,
    't': 1.0,
}

vault = Shape.from_library(data_shape)
form.envelope_from_shape(vault)
form.selfweight_from_shape(vault)
swt0 = vault.compute_selfweight()

# INITIALISE OR NOT
form = form.initialise_tna(plot=False)

# for key in form.vertices_where({'is_fixed': True}):
#     form.vertex_attribute(key, 'z', 0.5)
#     print(form.vertex_attribute(key, 'pz'))
#     form.vertex_attribute(key, 'b', [1.0, 1.0])

# rollers_ratio = [0.05, 0.05]
# form.set_boundary_rollers(total_rx=[rollers_ratio[0]*swt0]*2, total_ry=[rollers_ratio[1]*swt0]*2)

print('Initial Selfweight:', swt0)
plot_form(form).show()

### ---------- RUN

optimiser.data['calculation_type'] = 'tensorial'
analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
# analysis.run()

path = compas_tno.get('test.json')
form.to_json(path)

# print(out)

i_k = form.index_key()
i_uv = form.index_uv()
k_i = form.key_index()
uv_i = form.uv_index()

# Setting inds by hand
# ind = [21, 23, 32, 35, 41, 46, 47]
# gkeys = []
# for i in ind:
#     u, v = i_uv[i]
#     gkeys.append(geometric_key(form.edge_midpoint(u, v)[:2] + [0]))
#     form.edge_attribute((u, v), 'is_ind', value=True)
# form.attributes['indset'] = gkeys

find_inds = optimiser.data['find_inds']
printout = optimiser.data['printout']
objective = optimiser.data['objective']
variables = optimiser.data['variables']
qmax = optimiser.data['qmax']
qmin = optimiser.data['qmin']
constraints = optimiser.data['constraints']
indset = form.attributes['indset']

args = initialise_problem(form, indset=indset, printout=printout, find_inds=find_inds)
q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty = args

if 'reac_bounds' in constraints:
    b = set_b_constraint(form, True, True)
else:
    b = None

if 'cracks' in constraints:
    cracks_lb, cracks_ub = set_cracks_constraint(form, True, True)
else:
    cracks_lb, cracks_ub = None, None

if 'rollers' in constraints:
    max_rol_rx, max_rol_ry = set_rollers_constraint(form, True)
else:
    max_rol_rx, max_rol_ry = None, None

if 'joints' in constraints:
    joints = set_joints_constraint(form, True)
else:
    joints = None

if any(el in ['symmetry', 'symmetry-horizontal', 'symmetry-vertical'] for el in constraints):
    Asym = set_symmetry_constraint(form, constraints, True)
else:
    Asym = None

if objective == 'min':
    fobj = f_min_thrust
elif objective == 'max':
    fobj = f_max_thrust

args = (q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, cracks_lb, cracks_ub, free_x, free_y, rol_x, rol_y, Citx, City, Cftx, Cfty, qmin, constraints, max_rol_rx, max_rol_ry, Asym)

fconstr = constr_wrapper

x0 = q[ind]
zb_bounds = [[form.vertex_attribute(i_k[i], 'lb'), form.vertex_attribute(i_k[i], 'ub')] for i in fixed]
bounds = [[-10e-6, qmax]] * k + zb_bounds
x0 = append(x0, z[fixed]).reshape(-1, 1)

# plotter = MeshPlotter(form, figsize=(8,8))
# plotter.draw_edges(text={(u,v): uv_i[(u,v)] for u,v in form.edges()})
# plotter.show()

# plot_independents(form).show()

# CHECK ALSO CASE WITH ZB NON VARIABLE

print('Total of Independents:', len(ind))
print('Number of Variables:', len(x0))
f0 = fobj(x0, *args)
g0 = fconstr(x0, *args)
print('Number of Constraints:', len(g0))

print('Non Linear Optimisation - Initial Objective Value: {0}'.format(f0))
print('Non Linear Optimisation - Initial Constraints Extremes: {0:.3f} to {1:.3f}'.format(max(g0), min(g0)))

#
#  ------------------- Get the torch derivatives
#

EdinvEi = Edinv*Ei
Edinv_p = Edinv.dot(p)

EdinvEi_th = tensor(EdinvEi)
Edinv_p_th = tensor(Edinv_p)
C_th = tensor(C.toarray())
Ci_th = tensor(Ci.toarray())
Cit_th = Ci_th.t()
Cf_th = tensor(Cf.toarray())
pzfree = tensor(pz[free])
xyz = tensor(hstack([x, y, z]))
xy = tensor(hstack([x, y]))
pfixed = tensor(hstack([px, py, pz])[fixed])
U_th = tensor(U.toarray())
V_th = tensor(V.toarray())

args_obj = (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, objective)
args_constr = (Edinv_p_th, EdinvEi_th, ind, dep, C_th, Ci_th, Cit_th, Cf_th, pzfree, xyz, xy, pfixed, k, free, fixed, ub, lb, ub_ind, lb_ind, b, constraints, max_rol_rx, max_rol_ry, rol_x, rol_y, px, py, Asym, U_th, V_th)

print('shape variables:', x0.shape)
print('shape inds:', len(ind))
print('shape fixed:', len(pfixed))

# Restart Variables for Jacobian
variables = tensor(x0, requires_grad=True)
print('variables tensor shape', variables.shape)
g0 = f_constraints_pytorch(variables, *args_constr)
print('g0 shape', g0.shape)
jac = compute_jacobian(variables, g0)
print('jacobian shape', jac.shape)

# Restart Variables for Gradient
variables = tensor(x0, requires_grad=True)
print('len free, ub and lb', len(free), len(ub_ind), len(lb_ind))
f = f_objective_pytorch(variables, *args_obj)
print('f0: ', f)
grad = compute_grad(variables, f)
print('shape gradient', grad.shape)

# print('Jac Rx/Ry dq and dz')
# start = len(dep)+len(ub_ind)+len(ub_ind)
# end = start + 2*len(fixed)
# print(start, end)
# A = jac[start:end]
# for i in range(A.shape[0]):
#     max_ = max(A[i])
#     min_ = min(A[i])
#     print(max_, min_)
# plt.matshow(A)
# plt.colorbar()
# plt.show()

# print('Jac rollers')
# startB = len(dep)+len(ub_ind)+len(ub_ind)
# endB = len(g0)
# print(startB, endB)
# B = jac[startB:endB]
# for i in range(B.shape[0]):
#     max_ = max(B[i])
#     min_ = min(B[i])
#     print(max_, min_)
# plt.matshow(B)
# plt.colorbar()
# plt.show()

#
# --------- BY HAND GRADIENT
#

qid = q[ind]
args_z = q, ind, dep, E, Edinv, Ei, C, Ct, Ci, Cit, Cf, U, V, p, px, py, pz, z, free, fixed, lh, sym, k, lb, ub, lb_ind, ub_ind, s, Wfree, x, y, b, joints, i_uv, k_i
z, l2, q, _ = zlq_from_qid(qid, args_z)

dict_constr = constraints
derivatives = zeros([0, 1])
Q = diags(q.ravel())
CtQC = Ct.dot(Q).dot(C)

xopt = x0

jac_hand = sensitivities_wrapper(x0, *args)
print('Jacobian Hand', jac_hand.shape)

# A_ = jac_hand[start:end]
# B_ = jac_hand[startB:endB]

constr = constr_wrapper_ipopt(x0, *args)
print('Constraints Hand', constr.shape)

grad_comp = set(set(range(len(jac))) - set(ind))
jac_equivalent = jac  # [list(grad_comp)]

# plt.matshow(jac_hand)
# plt.colorbar()
# plt.show()

# plt.matshow(jac)
# plt.colorbar()
# plt.show()

# COMPARISON

# print('Jac REACTIONS')
# error = 0.0
# B = B[:, :len(ind)]
# B_ = B_[:, :len(ind)]
# for i in range(B.shape[0]):
#     max_ = max(B[i])
#     min_ = min(B[i])
#     max__ = max(B_[i])
#     min__ = min(B_[i])
#     norm_i = norm(B[i] - B_[i])
#     error += norm_i
#     # print('diff here:', norm_i)
#     print('ratio', float(max__/max_), float(min__/min_))
#     # print('max/min jac:', float(max_), float(min_), 'max/min jac_hand:', max__, min__)
# print('Error Jac REACTIONS is:', error)
# print(B.shape, B_.shape)
# plt.matshow(hstack([B, B_]))
# plt.colorbar()
# plt.show()

# Comparison of REAC_BOUNDS

# print('Jac Rx/Ry dz')
# A2 = A[:, len(ind):]
# A2_ = A_[:, len(ind):]
# error = 0.0
# for i in range(A2.shape[0]):
#     max_ = max(A2[i])
#     min_ = min(A2[i])
#     max__ = max(A2_[i])
#     min__ = min(A2_[i])
#     norm_i = norm(A2[i] - A2_[i])
#     error += norm_i
#     # print('diff here:', norm_i)
#     # print('max/min jac:', float(max_), float(min_), 'max/min jac_hand:', max__, min__)
# print('Error Jac Rx/Ry dz is:', error)
# # plt.matshow(hstack([A2, A2_]))
# # plt.colorbar()
# # plt.show()
# print(z[fixed])

# print('Jac Rx/Ry dq')
# A1 = A[:, :len(ind)]
# A1_ = A_[:, :len(ind)]
# for i in range(A1.shape[0]):
#     max_ = max(A1[i])
#     min_ = min(A1[i])
#     max__ = max(A1_[i])
#     min__ = min(A1_[i])
#     norm_i = norm(A1[i] - A1_[i])
#     error += norm_i
#     # print('diff here:', norm_i)
#     print('ratio', float(max__/max_), float(min__/min_))
#     # print('max/min jac:', float(max_), float(min_), 'max/min jac_hand:', max__, min__)
# print('Error Jac Rx/Ry dz is:', error)
# plt.matshow(hstack([A1, A1_]))
# plt.colorbar()
# plt.show()

# plt.matshow(A_)
# plt.colorbar()
# plt.show()

# Comparison of All

# print(len(grad_comp))
# plt.matshow(hstack([jac_equivalent, jac_hand]))
# plt.colorbar()
# plt.show()

A = jac_equivalent
B = jac_hand
error = 0.0
for i in range(A.shape[0]):
    norm_i = norm(A[i]-B[i])
    error += norm_i
    # print('Error in column {0} equals: {1}'.format(i, norm_i))
print('---> Error TOTAL in JACOBIAN is:', error)

# Comparison of Qdep - OK!

# A = jac_equivalent[:len(dep)]
# B = jac_hand[:len(dep)]
# error = 0.0
# for i in range(A.shape[0]):
#     norm_i = norm(A[i]-B[i])
#     error += norm_i
# print('Error Qd/qind:', error)
# plt.matshow(hstack([A, B]))
# plt.colorbar()
# plt.show()

# for matrix in [A]:
#     for j in range(matrix.shape[1]):
#         jac_dict = {}
#         for i in range(matrix.shape[0]):
#             aij = matrix[i, j]
#             uv = i_uv[dep[i]]
#             jac_dict[uv] = round(float(aij), 2)
#         max_, min_ = max(jac_dict), min(jac_dict)
#         plotter = MeshPlotter(form, figsize=(8, 8))
#         plotter.draw_edges(text={(u, v): str(jac_dict[(u, v)]) for u, v in form.edges_where({'is_ind': False})}) # color={(u, v): i_to_red(abs(jac_dict[(u, v)])) for u, v in form.edges()}
#         plotter.show()

# Comparison of dz/dqind -

A = jac_equivalent[len(dep):(len(dep)+len(lb_ind)), :len(ind)]
B = jac_hand[len(dep):(len(dep)+len(lb_ind)), :len(ind)]
print('Shape of A, B', A.shape, B.shape, form.number_of_edges())
for matrix in [A]:
    for j in range(matrix.shape[1]):
        jac_dict = {}
        for i in range(matrix.shape[0]):
            aij = matrix[i, j]
            jac_dict[k_i[free[i]]] = float(aij)
        max_, min_ = max(jac_dict.values()), min(jac_dict.values())
        print(min_, max_)
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges()
        plotter.draw_edges(keys=[i_uv[ind[j]]], color='FF0000', width=4)
        plotter.draw_vertices(text={key: '{0:.1f}'.format((jac_dict[key]-min_)/(max_ - min_)) for key in form.vertices_where({'is_fixed': False})}, facecolor={key: i_to_red((jac_dict[key]-min_)/(max_ - min_)) for key in form.vertices_where({'is_fixed': False})})
        plotter.show()
error = 0.0
for i in range(A.shape[0]):
    norm_i = norm(A[i]-B[i])
    error += norm_i
print('Error is dz/dqind:', error)
plt.matshow(hstack([A, B]))
plt.colorbar()
plt.show()

# Comparison of dz/dzf -

# print('q inside')
# print(q)

# print(lb_ind)
# print(free)

A = jac_equivalent[len(dep):(len(dep)+len(lb_ind)), len(ind):]
B = jac_hand[len(dep):(len(dep)+len(lb_ind)), len(ind):]
print('Shape of A, B', A.shape, B.shape, form.number_of_vertices())
for matrix in [A, B]:
    for j in range(matrix.shape[1]):
        jac_dict = {}
        for i in range(matrix.shape[0]):
            aij = matrix[i, j]
            jac_dict[k_i[free[i]]] = round(float(aij), 2)
        plotter = MeshPlotter(form, figsize=(10, 10))
        plotter.draw_edges()
        plotter.draw_vertices(text={key: str(jac_dict[key]) for key in form.vertices_where({'is_fixed': False})}, facecolor={key: i_to_red(abs(jac_dict[key])) for key in form.vertices_where({'is_fixed': False})})
        plotter.show()
error = 0.0
for i in range(A.shape[0]):
    norm_i = norm(A[i]-B[i])
    error += norm_i
print('Error is dz/dqind:', error)
plt.matshow(hstack([A, B]))
plt.colorbar()
plt.show()

# Comparison of SYM -

# A = jac_equivalent[-1*Asym.shape[0]:]
# B = jac_hand[-1*Asym.shape[0]:]
# error = 0.0
# for i in range(A.shape[0]):
#     norm_i = norm(A[i]-B[i])
#     error += norm_i
# print('Error is:', error)
# plt.matshow(hstack([A, B]))
# plt.colorbar()
# plt.show()


grad_hand = gradient_fmin(x0, *args)
highlights = []
print('Gradient comparison (Hand/AD)')
for i in range(len(grad_hand)):
    ghi = float(grad_hand[i])
    gi = float(grad[i])
    if abs(gi) > 10e-20:
        div = ghi/gi
    else:
        div = 0.0
    dif = abs(ghi-gi)
    if dif < 1e-6:
        dif = 0.0
    else:
        highlights.append(i)
    print('{0} : {1:.3f} | {2:.3f} | {3:.3f} | {4}'.format(i, ghi, gi, div, dif))
print('1/n_supports:', 1/len(fixed))

plt.matshow(hstack([grad_hand, grad]))
plt.colorbar()
plt.show()

plot_independents(form, highlights=highlights).show()

# print('Gradient Hand')
# print(grad_hand.flatten())
# print('Gradient AD')
# print(array(grad).flatten())
# print('divide')
# print(grad_hand[:len(ind)]/grad[:len(ind)])
error_grad = norm(grad_hand-array(grad))
error_grad_ind = norm(grad_hand[:len(ind)]-array(grad)[:len(ind)])
print('Error Grad: ', error_grad)
print('Error Grad (IND): ', error_grad_ind)

# normtotal = 0
# for i in range(len(dep), len(dep) + 2):
#     norm_i = norm(jac_equivalent[i]-jac_hand[i])
#     print(i)
#     print(norm_i)
#     if norm_i > 0.0:
#         print(jac_equivalent[i])
#         print(jac_hand[i])
#     normtotal += norm_i
# print('Norm diff', normtotal)

# TODO
# Fiz constraint dq/dz
# Code the constraint Rx/zf
