from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.problems import initialise_form
from compas_tno.datastructures import MeshDos
from compas_tno.viewers import view_mesh
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_thrust
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_solution
from compas_tno.algorithms import z_from_form
from compas_tno.algorithms import zq_from_qid

import compas_rv2

import json
import compas
# from compas_rv2.datastructures import FormDiagram  # this can also just be a mesh
FILE = '/Users/mricardo/compas_dev/me/freeform/Gene/formfinding.rv2'
with open(FILE, 'r') as f:
    session = json.load(f)
data = session['data']

for key in data:
    print(key)

form = FormDiagram.from_data(data['form'])

form.delete_face(363)
form.delete_face(362)

MESH = '/Users/mricardo/compas_dev/me/freeform/Gene/mesh.json'
# middle = MeshDos.from_mesh(form)
middle = MeshDos.from_json(MESH)
middle.delete_face(8)
middle.store_normals()

count_z = 0
count_keys = []
for key in middle.vertices():
    if middle.vertex_attribute(key, 'target') is None:
        z = middle.vertex_attribute(key, 'z')
        middle.vertex_attribute(key, 'target', z)
        count_z += 1
        count_keys.append(key)

print(middle.vertex_attribute(126, 'target'))
print('count z', count_z)

# plotter = MeshPlotter(middle)
# plotter.draw_faces(text={key: key for key in middle.faces()})
# plotter.draw_edges()
# plotter.draw_vertices(text={key: key for key in middle.vertices()})
# plotter.show()

for key in form.vertices():
    pz = form.vertex_area(key)
    form.vertex_attribute(key, 'pz', pz)
    if form.vertex_attribute(key, 'is_anchor'):
        form.vertex_attribute(key, 'is_fixed', True)

form = z_from_form(form)
# plot_form(form, show_q=False, simple=True, radius=0.4).show()

# # ------ This intiate and save
# form.initialise_loadpath(printout=False)
# form.overview_forces()
# path = '/Users/mricardo/compas_dev/me/freeform/Gene/form_lp.json'
# form.to_json(path)
# # ------

path = '/Users/mricardo/compas_dev/me/freeform/Gene/form_lp.json'
form = FormDiagram.from_json(path)

qs = []
for u, v in form.edges():
    q = form.edge_attribute((u, v), 'q')
    qs.append(q)

print('q max/min:', max(qs), min(qs))

ind = []
i = 0
for u, v in form.edges():
    if form.edge_attribute((u, v), 'is_ind'):
        ind.append(i)
    i+=1

print(ind)

# plot_form(form, show_q=False, simple=True, radius=0.4).show()
# view_mesh(form).show()

thk = 2.0
gradients = True

# intrados = middle.offset_mesh(n=thk/2, direction='down')
# extrados = middle.offset_mesh(n=thk/2, direction='up')
# data = {'thk': thk, 'type': 'general', 't': 0.0}
# shape = Shape.from_meshes(intrados, extrados, middle=middle, data=data)
# shape.ro = 0.5

pts = []# = middle.vertices_attributes('xyz')
height = middle.vertices_attributes('target')
# print(height)
print(len(pts), len(height))
i = 0
for key in middle.vertices():
    x, y, _ = middle.vertex_coordinates(key)
    z = middle.vertex_attribute(key, 'target')
    pts.append([x, y, z])

# for i in range(len(pts)):
#     pts[i][2] = height[i]

# print(pts)

shape = Shape.from_middle_pointcloud(pts, topology=form, thk=thk, printout=True)

# view_shapes(shape).show()

# plot_form(form, show_q=False, simple=True).show()
# view_thrust(form).show()

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'SLSQP'
optimiser.data['constraints'] = ['funicular', 'envelope']
# optimiser.data['constraints'] = ['funicular', 'envelope', 'symmetry']
optimiser.data['variables'] = ['ind', 't']
optimiser.data['objective'] = 't'
optimiser.data['thickness_type'] = 'constant'
optimiser.data['min_thk'] = 0.01
optimiser.data['max_thk'] = thk*1.0
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients

form.overview_forces()

analysis = Analysis.from_elements(shape, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

SAVE = '/Users/mricardo/compas_dev/me/freeform/Gene/solution.json'
form.to_json(SAVE)

plot_form(form, show_q=False, simple=True, radius=0.4).show()

view_solution(form, shape).show()
