from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import plot_form
from compas_tno.problems import initialise_form
from compas_tno.shapes import MeshDos
from compas_tno.viewers import view_mesh
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_thrust
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_solution
from compas_tno.algorithms import z_from_form

# path = '/Users/mricardo/compas_dev/me/freeform/Gene/mesh.json'
# mesh = Mesh.from_json(path)
# # form = FormDiagram.from_mesh(mesh)

# path_form = '/Users/mricardo/compas_dev/me/freeform/Gene/form.json'
# form = FormDiagram.from_json(path_form)

# # form.edges_attributes(names='q', values=1.0)
# for key in form.vertices():
#     pz = form.vertex_area(key)
#     form.vertex_attribute(key, 'pz', pz)
#     _, _, z_mesh = mesh.vertex_coordinates(key)
#     form.vertex_attribute(key, 'target', z_mesh)
#     if mesh.vertex_attribute(key, 'is_fixed'):
#         form.vertex_attribute(key, 'is_fixed', True)
# for u, v in form.edges():
#     form.edge_attribute((u, v), 'q', 1.0)
#     if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
#         form.edge_attribute((u, v), '_is_edge', False)

# plot_form(form, show_q=False, simple=True, radius=0.3, max_width=5).show()

# form = initialise_form(form, printout=True)

# middle = MeshDos.from_formdiagram_attribute(form, attribute='target')
# middle.to_json('/Users/mricardo/compas_dev/me/freeform/Gene/target_mesh.json')

path_form2 = '/Users/mricardo/compas_dev/me/freeform/Gene/form_init.json'
path_lp = '/Users/mricardo/compas_dev/me/freeform/Gene/form_lp.json'
# form.to_json(path_form2)

# form = FormDiagram.from_json(path_form2)

# pzt = 0
# pz = []
# for key in form.vertices():
#     pz = form.vertex_attribute(key, 'pz')
#     pzt += pz
# print('total weight:', pzt)

# form.initialise_loadpath()
# plot_form(form, show_q=False, simple=True).show()

# q = []
# for u, v in form.edges():
#     q = form.edge_attribute((u, v), 'q')
#     form.edge_attribute((u, v), 'q', q/1000)
#     q.append(q/1000)

# form = z_from_form(form)

# print('max/min q:', max(q), min(q))
# print('max/min pz:', max(pz), min(pz))
# plot_form(form, show_q=False, simple=True).show()
# form.to_json(path_lp)

# view_thrust(form).show()

path_middle = '/Users/mricardo/compas_dev/me/freeform/Gene/target_mesh.json'
path_lp = '/Users/mricardo/compas_dev/me/freeform/Gene/form_lp.json'
form = FormDiagram.from_json(path_lp)
form.overview_forces()

middle = MeshDos.from_json(path_middle)
middle.store_normals()

thk = 2.0
gradients = True

intrados = middle.offset_mesh(n=thk/2, direction='down')
extrados = middle.offset_mesh(n=thk/2, direction='up')
data = {'thk': thk, 'type': 'general', 't': 0.0}
shape = Shape.from_meshes(intrados, extrados, middle=middle, data=data)
shape.ro = 1.0

# view_shapes(shape).show()

# plot_form(form, show_q=False, simple=True).show()
# view_thrust(form).show()

optimiser = Optimiser()
optimiser.data['library'] = 'IPOPT'
optimiser.data['solver'] = 'IPOPT'
optimiser.data['constraints'] = ['funicular', 'envelope']
# optimiser.data['constraints'] = ['funicular', 'envelope', 'symmetry']
optimiser.data['variables'] = ['ind', 'zb', 't']
optimiser.data['objective'] = 't'
optimiser.data['thickness_type'] = 'constant'
optimiser.data['min_thk'] = 0.0
optimiser.data['max_thk'] = thk*2.0
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients

analysis = Analysis.from_elements(shape, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

view_solution(form, shape).show()
