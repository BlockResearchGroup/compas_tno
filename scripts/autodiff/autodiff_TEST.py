
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser

from compas_tno.problems import initialize_loadpath
from compas_tno.problems import constr_wrapper
from compas_tno.problems import sensitivities_wrapper

from compas_tno.autodiff.jax_swt import weights_from_xyz_jax
from compas_tno.autodiff.jax_objectives import f_min_thrust_jax
from compas_tno.autodiff.jax_constraints import f_constraints_jax
from compas_tno.autodiff.jax_equilibrium import xyz_from_q_jax
from compas_tno.autodiff.jax_equilibrium import z_from_q_jax
from compas_tno.autodiff.jax_constraints import f_constraints_jax
from compas_tno.autodiff.jax_equilibrium import q_from_variables_jax

from compas_tno.algorithms import weights_from_xyz

from jax.numpy import diag
from jax.numpy import vstack

from jax.numpy import dot
from numpy import zeros
from numpy import identity

from jax.numpy import array

from jax import jacfwd, jacrev, jit, grad

import matplotlib.pyplot as plt
import numpy as np


def make_form_shape():
    # Parameters Geometries

    type_formdiagram = 'cross_fd'
    type_structure = 'crossvault'
    thk = 0.5
    discretisation = 10
    discretisation_shape = discretisation * 2

    # Create form diagram

    data_diagram = {
        'type': type_formdiagram,
        'xy_span': [[0.0, 10.0], [0.0, 10.0]],
        'fix': 'corners',
        'discretisation': discretisation
    }

    form = FormDiagram.from_library(data_diagram)

    # Create shape

    data_shape = {
        'type': type_structure,
        'thk': thk,
        'discretisation': discretisation_shape,
        'xy_span': [[0.0, 10.0], [0.0, 10.0]],
        't': 0.0,
    }

    dome = Shape.from_library(data_shape)

    return form, dome


def z_from_variables(variables, M):

    k = M.k
    nb = M.nb
    check = k
    qid = variables[:k].reshape(-1, 1)
    pz = M.P[array(M.free),  2].reshape(-1, 1)
    zb = M.X[array(M.fixed), 2].reshape(-1, 1)

    q = dot(M.B, qid) + M.d

    if 'zb' in M.variables:
        zb = variables[check: check + nb]
    # Xfree = xyz_from_q_jax(q, P[array(M.free)], X[array(M.fixed)], M.Ci, M.Cit, M.Cb)
    zfree = z_from_q_jax(q, pz, zb, M.Ci, M.Cit, M.Cb)
    z = dot(M.Pmatrix, vstack([zfree, zb]))

    return z


obj = 'min'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']
starting_point = 'current'

form, shape = make_form_shape()

F, V0, V1, V2 = form.tributary_matrices()

xyz = array(form.vertices_attributes('xyz'))
pz = weights_from_xyz(xyz, F, V0, V1, V2)

print('pz plan:', sum(pz))
print(pz.shape)

i = 0
for key in form.vertices():
    form.vertex_attribute(key, 'pz', -1 * pz[i])
    i += 1

M = initialize_loadpath(form)

xyz = M.X

pz = weights_from_xyz(xyz, F, V0, V1, V2)

print('pz 3d:', sum(pz))
print(pz.shape)

optimiser = Optimiser()
optimiser.settings['objective'] = obj
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['starting_point'] = starting_point
optimiser.settings['printout'] = True
optimiser.settings['derivative_test'] = True
optimiser.settings['autodiff'] = True

analysis = Analysis.from_elements(shape, form, optimiser)
analysis.apply_envelope()
analysis.set_up_optimiser()

M = optimiser.M
x0 = optimiser.x0

xyz = M.X

M.V0 = V0
M.V1 = V1
M.V2 = V2
M.F = F

# pz_jax = weights_from_xyz_jax(xyz, F, V0, V1, V2)
# dx = jacfwd(weights_from_xyz_jax)
# dpzdx = dx(xyz, F, V0, V1, V2)

# print('jax weights', sum(pz_jax), max(pz_jax), min(pz_jax))

# print(dpzdx.shape, xyz.shape)

# pz_anal = weights_from_xyz(xyz, F, V0, V1, V2)

# print('anal weights', sum(pz_anal), max(pz_anal), min(pz_anal))

# print(M.q)

z_anal = M.X[:, 2]
jac_z = jacfwd(z_from_variables)
z_jac_jax = jac_z(x0, M)[:, 0, :, 0]

print(z_jac_jax.shape)

n = M.n
m = M.m

jac_anal = sensitivities_wrapper(x0, M)
z_jac_anal = jac_anal[2*m:2*m+n, :]

print(z_jac_jax.shape, z_jac_anal.shape)

print('max z jac anal:', max(z_jac_anal.flatten()))
print('min z jac anal:', min(z_jac_anal.flatten()))
print('max z jac jax:', max(z_jac_jax.flatten()))
print('min z jac jax:', min(z_jac_jax.flatten()))

diff_jac_z = abs(z_jac_anal - z_jac_jax) ** 2

print('Diff jacobian:', sum(sum(diff_jac_z)))
print('-'*20)

plt.matshow(diff_jac_z)
plt.show()

constraints_jax = f_constraints_jax(x0, M)
f_jacobian_jax_ = jacfwd(f_constraints_jax)
constraint_anal = constr_wrapper(x0, M)

print('constr shape', constraint_anal.shape, constraints_jax.shape)

print('Diff constr', sum(constraints_jax - constraint_anal)**2)

jacobian_jax = f_jacobian_jax_(x0, M)[:, :, 0]
jacobian_anal = sensitivities_wrapper(x0, M)

print('max constr jac jax:', max(jacobian_jax.flatten()))
print('min constr jac jax:', min(jacobian_jax.flatten()))
print('max constr jac anal:', max(jacobian_anal.flatten()))
print('min constr jac anal:', min(jacobian_anal.flatten()))

print(jacobian_jax.shape, jacobian_anal.shape)

diff_jax = abs(jacobian_jax - jacobian_anal)

print(sum(sum((jacobian_jax - jacobian_anal)**2)))  # jacobian has differences!

import matplotlib.pyplot as plt
import numpy as np

plt.matshow(diff_jax)
plt.show()

# xyz_noupdate = xyz_from_q_jax(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)
# dx_noupdate = jacfwd(xyz_from_q_jax)
# dxdq_noupdate = dx_noupdate(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)[:, :, :, 0]
# print(dxdq_noupdate.shape, M.q.shape)

# # xyz_noupdate = xyz_from_q_jax(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)
# # dx_noupdate = jacfwd(xyz_from_q_jax)
# # dxdq_noupdate = dx_noupdate(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)[:, :, :, 0]
# # print(dxdq_noupdate.shape, M.q.shape)

# for i, key in enumerate(form.vertices_where({'is_fixed': False})):
#     x, y, z = xyz_noupdate[i]
#     form.vertex_attribute(key, 'x', x)
#     form.vertex_attribute(key, 'y', y)
#     form.vertex_attribute(key, 'z', z)

# view = Viewer(form)
# view.draw_thrust()
# view.draw_shape()
# view.show()

# xyz_withupdate = xyz_from_q_jax(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb, update_loads=True, M=M)
# dx_withupdate = jacfwd(xyz_from_q_jax)
# dxdq_withupdate = dx_withupdate(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb, update_loads=True, M=M)[:, :, :, 0]
# print(dxdq_withupdate.shape, M.q.shape)

# for i, key in enumerate(form.vertices_where({'is_fixed': False})):
#     x, y, z = xyz_withupdate[i]
#     form.vertex_attribute(key, 'x', x)
#     form.vertex_attribute(key, 'y', y)
#     form.vertex_attribute(key, 'z', z)

# view = Viewer(form)
# view.draw_thrust()
# view.draw_shape()
# view.show()

# from numpy import savez
# fn = '/Users/mricardo/compas_dev/compas_tno/data/test.npz'
# savez(fn, xyz_noupdate=xyz_noupdate, xyz_withupdate=xyz_withupdate, dxdq_noupdate=dxdq_noupdate, dxdq_withupdate=dxdq_withupdate)
# print('Saved matrices to:', fn)

# # constr = f_constraints_jax(x0, M)

# # jac = jacfwd(f_constraints_jax)
# # j = jac(x0, M)[:, :, 0]

# # print(j.shape)

# # print(j)

# # # ANALYTIC

# # fobj = optimiser.fobj
# # fgrad = optimiser.fgrad

# # fmin = fobj(x0, M)
# # print('analytic', fmin)

# # fmingrad_anal = fgrad(x0, M)
# # print(fmingrad_anal)

# # analysis.run()

# # fmin = f_min_thrust_jax(x0, M)
# # print('jax', fmin)

# # fgrad = grad(f_min_thrust_jax)
# # fmingrad_jax = fgrad(x0, M)

# # print(fmingrad_jax)

# # diff = abs(fmingrad_jax.flatten() - fmingrad_anal.flatten())
# # print('Diff', max(diff))

# # view = Viewer(form)
# # view.draw_thrust()
# # view.draw_shape()
# # view.show()

# # jac_hand = deriv_weights_from_matrices(xyz, F, V0, V1, V2, features=[])

# # jac_X = jacfwd(weights_from_matrices_jax)

# # jacsol = jac_X(xyz, F, V0, V1, V2)

# # print('jac_hand.shape', jac_hand.shape)
# # print('jacsol.shape', jacsol.shape)

# # import compas_tno
# # outfile = compas_tno.get('test.npz')

# # from numpy import savez
# # savez(outfile, A=jac_hand, B=jacsol)

# # import compas_tno
# # outfile = compas_tno.get('test.npz')

# # from numpy import load
# # matrices = load(outfile)
# # A = matrices['A']
# # B = matrices['B']

# # jac_xyz = jacfwd(weights_from_matrices_jax_xyz)

# # jacsol = jac_xyz(x, y, z, F, V0, V1, V2)

# # print(jacsol.shape) ##  check this dimension

# # print(jacsol[1])

# # print(xyz)

# # @partial(jit, static_argnums=1)
# # def weights_from_xyz(xyz, tributary_dict):

# #     pz = zeros((len(xyz), 1))

# #     for i in tributary_dict:
# #         p0 = xyz[i]
# #         area = 0
# #         for j in tributary_dict[i]:
# #             pj = xyz[j]
# #             for face in tributary_dict[i][j]:
# #                 face_coord = [xyz[pt] for pt in face]
# #                 pf = centroid_points(face_coord)
# #                 v1 = subtract_vectors(pj, p0)
# #                 v2 = subtract_vectors(pf, p0)
# #                 area += length_vector(cross_vectors(v1, v2))
# #         pz[i] = 0.25 * area

# #     return pz

# # weights_from_xyz(xyz, tributary_dict)
# # a = [0, 0]
# # a = {0:0, 1:0}
# # pz = test(xyz, a)
# # print(pz)
# # grad_pz = grad(test)

# # gradval= grad_pz(xyz, a)

# # print(gradval)

# # pz = weights_from_x_y_z(x, y, z, tributary_dict)

# # print(sum(pz))

# # jac = jacfwd(weights_from_x_y_z)(x, y, z, tributary_dict)

# # print(jac)
# # print(jac.shape)

# # pz = weights_from_xyz(xyz, tributary_dict)

# # print(sum(pz))

# # grad_pz = grad(weights_from_xyz, argnums=1)

# # dpzdX = grad_pz(xyz, tributary_dict)

# # print(dpzdX)

# # Parameters Optimisations

# # obj = 'min'
# # solver = 'SLSQP'
# # constraints = ['funicular', 'envelope', 'reac_bounds']
# # variables = ['q', 'zb']
# # features = ['fixed']
# # starting_point = 'loadpath'

# # Optimisation

# # view = Viewer(form, dome)
# # view.draw_thrust()
# # view.draw_shape()
# # view.show()

# # opt = Optimiser()
# # opt.settings['objective'] = obj
# # opt.settings['solver'] = solver
# # opt.settings['constraints'] = constraints
# # opt.settings['variables'] = variables
# # opt.settings['features'] = features
# # opt.settings['starting_point'] = starting_point
# # opt.settings['printout'] = True

# # print(dome.datashape)

# # analysis = Analysis.from_elements(dome, form, opt)
# # analysis.apply_selfweight()
# # analysis.apply_envelope()
# # analysis.apply_reaction_bounds()
# # analysis.set_up_optimiser()

# # ------------

# # M = opt.M
# # x0 = opt.x0


# # from compas_tno.problems import f_min_thrust
# # from compas_tno.problems import gradient_fmin

# # print('-'* 10, 'FROM ANALYTIC', '-'* 10)
# # pt1 = f_min_thrust(x0, M)
# # a1 = gradient_fmin(x0, M)
# # print(a1)
# # print(pt1)

# # # make non sparse

# # M.C = M.C.toarray()
# # M.Ci = M.Ci.toarray()
# # M.Cit = M.Cit.toarray()
# # M.Cb = M.Cb.toarray()

# # from compas_tno.autodiff.jax import xyz_from_q_jax
# # from compas_tno.autodiff.jax import f_min_thrust_jax
# # from jax import grad

# # # grad_xyz = grad(xyz_from_q_jax)
# # grad_min = grad(f_min_thrust_jax)

# # print('-'* 10, 'FROM JAX', '-'* 10)
# # pt2 = f_min_thrust_jax(x0, M)
# # a2 = grad_min(x0, M)
# # print(a2)
# # print(pt2)
# # # b = grad_xyz(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)


# # print(a1.shape)
# # print(a2.shape)

# # diff = abs(a1.reshape(-1, 1) - a2)
# # print(max(diff))

# # # print(a)

# # # analysis.run()

# # # view = Viewer(form, dome)
# # # view.draw_thrust()
# # # view.draw_shape()
# # # view.draw_cracks()
# # # view.show()
