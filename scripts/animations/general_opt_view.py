
formadress = '/Users/mricardo/compas_dev/me/freeform/tom/form_solution_iter200.json'
forceadress = '/Users/mricardo/compas_dev/me/freeform/tom/force_solution_iter200.json'

formadress_1200 = '/Users/mricardo/compas_dev/me/freeform/tom/form_solution_iter1200.json'
forceadress_1200 = '/Users/mricardo/compas_dev/me/freeform/tom/force_solution_iter1200.json'

file_qs = '/Users/mricardo/compas_dev/compas_tno/data/output_200it.json'
file_qs_add = '/Users/mricardo/compas_dev/compas_tno/data/output_1000it.json'
file_Xform = '/Users/mricardo/compas_dev/compas_tno/data/output_200it_form.json'
file_Xforce = '/Users/mricardo/compas_dev/compas_tno/data/output_200it_force.json'


# from compas_tno.diagrams import FormDiagram
# from compas_tno.diagrams import ForceDiagram
# from compas.datastructures import Mesh
# from compas_tno.optimisers import Optimiser
# from compas_tno.analysis import Analysis
# from compas_tno.shapes import Shape

# from compas_plotters import MeshPlotter

# from compas_tno.viewers import view_shapes
# from compas_tno.viewers import view_thrust
# from compas_tno.viewers import view_bestfit_solution
# from compas_tno.viewers import view_thrust_as_lines
# from compas_tno.plotters import plot_form
# from compas_tno.plotters import plot_superimposed_diagrams
# from compas_tno.algorithms import equilibrium_fdm
# import compas
# import json


# form = FormDiagram.from_json(formadress)
# force = ForceDiagram.from_formdiagram(form)
# key_index = form.key_index()
# _key_index = force.key_index()

# form, force = form.reciprocal_from_form(plot=False)

# with open(file_qs, mode='r', encoding='utf-8') as f:
#     data = json.load(f)

# Xform = {}
# Xforce = {}

# iterations = len(data['iterations'])
# print('initial for # it:', iterations)

# from compas_tno.problems import initialise_problem_general
# from numpy import array
# M = initialise_problem_general(form)
# from compas_tno.algorithms import xyz_from_q

# for it in range(iterations):
#     print(it)
#     q = array(data['iterations'][str(it)]).reshape(-1, 1)
#     M.X[M.free] = xyz_from_q(q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)
#     for key in form.vertices():
#         index = key_index[key]
#         x, y, z = M.X[index]
#         form.vertex_attribute(key, 'x', x)
#         form.vertex_attribute(key, 'y', y)
#         form.vertex_attribute(key, 'z', z)
#     for index, (u, v) in enumerate(form.edges_where({'_is_edge': True})):
#         qi = q[index].item()
#         form.edge_attribute((u, v), 'q', qi)
#     form, force = form.reciprocal_from_form(plot=False)
#     Xforce[it] = force.vertices_attributes('xyz')
#     Xform[it] = M.X.tolist()

# print('end')

# with open(file_qs_add, mode='r', encoding='utf-8') as f:
#     data = json.load(f)

# iterations = len(data['iterations'])
# print('initial for # it:', iterations)

# for it in range(iterations):
#     it_ = it + 200
#     if it_ % 5 == 0 or it == iterations - 1:
#         pass
#     else:
#         continue
#     print(it_)
#     q = array(data['iterations'][str(it)]).reshape(-1, 1)
#     M.X[M.free] = xyz_from_q(q, M.P[M.free], M.X[M.fixed], M.Ci, M.Cit, M.Cb)
#     for key in form.vertices():
#         index = key_index[key]
#         x, y, z = M.X[index]
#         form.vertex_attribute(key, 'x', x)
#         form.vertex_attribute(key, 'y', y)
#         form.vertex_attribute(key, 'z', z)
#     for index, (u, v) in enumerate(form.edges_where({'_is_edge': True})):
#         qi = q[index].item()
#         form.edge_attribute((u, v), 'q', qi)
#     form, force = form.reciprocal_from_form(plot=False)
#     length = len(Xforce)
#     Xforce[length] = force.vertices_attributes('xyz')
#     Xform[length] = M.X.tolist()

# print(len(Xform))
# print(len(Xforce))
# print(length)

# with open(file_Xform, mode='w', encoding='utf-8') as f:
#     data = json.dump(Xform, f)

# with open(file_Xforce, mode='w', encoding='utf-8') as f:
#     data = json.dump(Xforce, f)

# print('Saved jsons to:')
# print(file_Xform)
# print(file_Xforce)

# form.to_json(formadress_1200)
# force.to_json(forceadress_1200)

# ----------------------------------------

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas.datastructures import Mesh
from compas_view2 import app
import time
import compas
import json

from compas_tno.viewers import animation_from_optimisation

form = FormDiagram.from_json(formadress)
force = ForceDiagram.from_json(forceadress)

animation_from_optimisation(form, file_Xform=file_Xform, force=force, file_Xforce=file_Xforce)

# form = Mesh.from_json(formadress)
# force = Mesh.from_json(forceadress)
# key_index = form.key_index()
# _key_index = force.key_index()

# # Simply scalling on (0,0,0)
# scale = 8.0
# y_trans = 6.0

# faces_delete = []

# for fkey in form.faces():
#     if form.face_attribute(fkey, '_is_loaded') is False:
#         faces_delete.append(fkey)

# for fkey in faces_delete:
#     form.delete_face(fkey)

# _X0 = force.vertices_attributes('xyz')
# for key in force.vertices():
#     index = _key_index[key]
#     force.vertex_attribute(key, 'x', _X0[index][0]/scale)
#     force.vertex_attribute(key, 'y', _X0[index][1]/scale+y_trans)
#     force.vertex_attribute(key, 'z', _X0[index][2]/scale)

# viewer = app.App()

# obj = viewer.add(form)
# _obj = viewer.add(force)

# with open(file_Xform, mode='r', encoding='utf-8') as f:
#     Xform = json.load(f)

# with open(file_Xforce, mode='r', encoding='utf-8') as f:
#     Xforce = json.load(f)

# iterations = len(Xform)
# print('Total iterations:', iterations)


# @viewer.on(interval=100, frames=iterations - 1)
# def update(f):
#     print(f)
#     if f == 0:
#         time.sleep(5)

#     Xf = Xform[str(f)]
#     _Xf = Xforce[str(f)]
#     for vertex in form.vertices():
#         index = key_index[vertex]
#         form.vertex_attribute(vertex, 'x', Xf[index][0])
#         form.vertex_attribute(vertex, 'y', Xf[index][1])
#         form.vertex_attribute(vertex, 'z', Xf[index][2])
#     for vertex in force.vertices():
#         index = _key_index[vertex]
#         force.vertex_attribute(vertex, 'x', _Xf[index][0]/scale)
#         force.vertex_attribute(vertex, 'y', _Xf[index][1]/scale+y_trans)
#         force.vertex_attribute(vertex, 'z', _Xf[index][2]/scale)
#     obj.update()
#     _obj.update()


# viewer.run()
