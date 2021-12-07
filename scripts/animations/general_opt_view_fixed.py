
formadress = '/Users/mricardo/compas_dev/me/freeform/tom/form_solution_iter200.json'
forceadress = '/Users/mricardo/compas_dev/me/freeform/tom/force_solution_iter200.json'

file_qs = '/Users/mricardo/compas_dev/compas_tno/data/output_fixed.json'

file_Xform = '/Users/mricardo/compas_dev/compas_tno/data/output_fixed_it_form.json'
file_Xforce = '/Users/mricardo/compas_dev/compas_tno/data/output_fixed_it_force.json'

form0fixed = '/Users/mricardo/compas_dev/me/freeform/tom/form0_fixed.json'
force0fixed = '/Users/mricardo/compas_dev/me/freeform/tom/force0_fixed.json'

# ----------------------------------------

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

# print(len(Xform))
# print(len(Xforce))

# with open(file_Xform, mode='w', encoding='utf-8') as f:
#     data = json.dump(Xform, f)

# with open(file_Xforce, mode='w', encoding='utf-8') as f:
#     data = json.dump(Xforce, f)

# print('Saved jsons to:')
# print(file_Xform)
# print(file_Xforce)

# ----------------------------------------

from compas_tna.diagrams import FormDiagram
from compas_tna.diagrams import ForceDiagram
from compas.datastructures import Mesh
from compas_view2 import app
import time
import compas
import json

form = Mesh.from_json(form0fixed)
force = Mesh.from_json(force0fixed)
key_index = form.key_index()
_key_index = force.key_index()

# Simply scalling on (0,0,0) - Already applied to force @ initial step
scale = 8.0
y_trans = 6.0

viewer = app.App()

obj = viewer.add(form)
_obj = viewer.add(force)

with open(file_Xform, mode='r', encoding='utf-8') as f:
    Xform = json.load(f)

with open(file_Xforce, mode='r', encoding='utf-8') as f:
    Xforce = json.load(f)

iterations = len(Xform)
iterations = 230  # Cap at iteration 230

@viewer.on(interval=100, frames=iterations - 1)
def update(f):
    print(f)
    if f == 0:
        time.sleep(5)

    Xf = Xform[str(f)]
    _Xf = Xforce[str(f)]
    for vertex in form.vertices():
        index = key_index[vertex]
        form.vertex_attribute(vertex, 'x', Xf[index][0])
        form.vertex_attribute(vertex, 'y', Xf[index][1])
        form.vertex_attribute(vertex, 'z', Xf[index][2])
    for vertex in force.vertices():
        index = _key_index[vertex]
        force.vertex_attribute(vertex, 'x', _Xf[index][0]/scale)
        force.vertex_attribute(vertex, 'y', _Xf[index][1]/scale+y_trans)
        force.vertex_attribute(vertex, 'z', _Xf[index][2]/scale)
    obj.update()
    _obj.update()


viewer.run()
