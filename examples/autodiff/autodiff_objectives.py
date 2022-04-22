from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer

from compas_tno.problems import initialize_loadpath

from compas.geometry import subtract_vectors
from compas.geometry import length_vector
from compas.geometry import cross_vectors
from compas.geometry import centroid_points

from numpy import tri, zeros

from compas_tno.utilities.loads import weights_from_x_y_z, deriv_weights_from_matrices
from compas_tno.utilities.loads import weights_from_matrices
from compas_tno.utilities.loads import test

from compas_tno.autodiff.jax_swt import weights_from_matrices_jax, weights_from_matrices_jax_xyz

from jax import grad
from numpy import array
from jax import jit

from functools import partial

from jax import jacfwd, jacrev


def make_form_shape():
    # Parameters Geometries

    type_formdiagram = 'cross_fd'
    type_structure = 'crossvault'
    thk = 0.5
    discretisation = 14
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


form, shape = make_form_shape()

# tributary_dict = form.tributary_dict()

# print('Number of vertices:', form.number_of_vertices())

F, V0, V1, V2 = form.tributary_matrices()

xyz = array(form.vertices_attributes('xyz'))
pz = weights_from_matrices(xyz, F, V0, V1, V2)

print('pz plan:', sum(pz))
print(pz.shape)

i = 0
for key in form.vertices():
    form.vertex_attribute(key, 'pz', -1 * pz[i])
    i += 1

M = initialize_loadpath(form)

xyz = M.X
pz = weights_from_matrices(xyz, F, V0, V1, V2)

print('pz 3d:', sum(pz))
print(pz.shape)

view = Viewer(form)
view.draw_thrust()
view.show()

jac_hand = deriv_weights_from_matrices(xyz, F, V0, V1, V2, features=[])

jac_X = jacfwd(weights_from_matrices_jax)

jacsol = jac_X(xyz, F, V0, V1, V2)

print('jac_hand.shape', jac_hand.shape)
print('jacsol.shape', jacsol.shape)

import compas_tno
outfile = compas_tno.get('test.npz')

from numpy import savez
savez(outfile, A=jac_hand, B=jacsol)

import compas_tno
outfile = compas_tno.get('test.npz')

from numpy import load
matrices = load(outfile)
A = matrices['A']
B = matrices['B']

# jac_xyz = jacfwd(weights_from_matrices_jax_xyz)

# jacsol = jac_xyz(x, y, z, F, V0, V1, V2)

# print(jacsol.shape) ##  check this dimension

# print(jacsol[1])

# print(xyz)

# @partial(jit, static_argnums=1)
# def weights_from_xyz(xyz, tributary_dict):

#     pz = zeros((len(xyz), 1))

#     for i in tributary_dict:
#         p0 = xyz[i]
#         area = 0
#         for j in tributary_dict[i]:
#             pj = xyz[j]
#             for face in tributary_dict[i][j]:
#                 face_coord = [xyz[pt] for pt in face]
#                 pf = centroid_points(face_coord)
#                 v1 = subtract_vectors(pj, p0)
#                 v2 = subtract_vectors(pf, p0)
#                 area += length_vector(cross_vectors(v1, v2))
#         pz[i] = 0.25 * area

#     return pz

# weights_from_xyz(xyz, tributary_dict)
# a = [0, 0]
# a = {0:0, 1:0}
# pz = test(xyz, a)
# print(pz)
# grad_pz = grad(test)

# gradval= grad_pz(xyz, a)

# print(gradval)

# pz = weights_from_x_y_z(x, y, z, tributary_dict)

# print(sum(pz))

# jac = jacfwd(weights_from_x_y_z)(x, y, z, tributary_dict)

# print(jac)
# print(jac.shape)

# pz = weights_from_xyz(xyz, tributary_dict)

# print(sum(pz))

# grad_pz = grad(weights_from_xyz, argnums=1)

# dpzdX = grad_pz(xyz, tributary_dict)

# print(dpzdX)

# Parameters Optimisations

# obj = 'min'
# solver = 'SLSQP'
# constraints = ['funicular', 'envelope', 'reac_bounds']
# variables = ['q', 'zb']
# features = ['fixed']
# starting_point = 'loadpath'

# Optimisation

# view = Viewer(form, dome)
# view.draw_thrust()
# view.draw_shape()
# view.show()

# opt = Optimiser()
# opt.settings['objective'] = obj
# opt.settings['solver'] = solver
# opt.settings['constraints'] = constraints
# opt.settings['variables'] = variables
# opt.settings['features'] = features
# opt.settings['starting_point'] = starting_point
# opt.settings['printout'] = True

# print(dome.datashape)

# analysis = Analysis.from_elements(dome, form, opt)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
# analysis.set_up_optimiser()

# ------------

# M = opt.M
# x0 = opt.x0


# from compas_tno.problems import f_min_thrust
# from compas_tno.problems import gradient_fmin

# print('-'* 10, 'FROM ANALYTIC', '-'* 10)
# pt1 = f_min_thrust(x0, M)
# a1 = gradient_fmin(x0, M)
# print(a1)
# print(pt1)

# # make non sparse

# M.C = M.C.toarray()
# M.Ci = M.Ci.toarray()
# M.Cit = M.Cit.toarray()
# M.Cb = M.Cb.toarray()

# from compas_tno.autodiff.jax import xyz_from_q_jax
# from compas_tno.autodiff.jax import f_min_thrust_jax
# from jax import grad

# # grad_xyz = grad(xyz_from_q_jax)
# grad_min = grad(f_min_thrust_jax)

# print('-'* 10, 'FROM JAX', '-'* 10)
# pt2 = f_min_thrust_jax(x0, M)
# a2 = grad_min(x0, M)
# print(a2)
# print(pt2)
# # b = grad_xyz(M.q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)


# print(a1.shape)
# print(a2.shape)

# diff = abs(a1.reshape(-1, 1) - a2)
# print(max(diff))

# # print(a)

# # analysis.run()

# # view = Viewer(form, dome)
# # view.draw_thrust()
# # view.draw_shape()
# # view.draw_cracks()
# # view.show()
