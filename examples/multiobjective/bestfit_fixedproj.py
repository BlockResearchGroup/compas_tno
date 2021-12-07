from compas_tno.diagrams import FormDiagram
from compas.datastructures import Mesh
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.shapes import Shape

from compas_plotters import MeshPlotter

from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_bestfit_solution
from compas_tno.viewers import view_thrust_as_lines
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.algorithms import equilibrium_fdm
import compas

# Load and print useful INFO

formtno = '/Users/mricardo/compas_dev/me/freeform/tom/form_solution_entropy.json'
solution = '/Users/mricardo/compas_dev/me/freeform/tom/form_solution_fixedproj.json'
form = FormDiagram.from_json(formtno)

bbox = form.bounding_box()
bbox_xy = form.bounding_box_xy()

print('Number of edges (form):', form.number_of_edges())
print('Number of vertices (form):', form.number_of_vertices())
print('Bounding Box:', bbox)
print('Bounding Box xy:', bbox_xy)

# # Add bounds as intra- extra-dos and create a SHAPE (WIP):

offset = 0.1
# for key in form.vertices():
#     zt = form.vertex_attribute(key, 'target')
#     form.vertex_attribute(key, 'ub', zt + offset)
#     form.vertex_attribute(key, 'lb', zt - offset)

# for u, v in form.edges_where({'_is_edge': True}):
#     q = form.edge_attribute((u, v), 'q')
#     form.edge_attribute((u, v), 'q', q * 15)
#     faces = form.edge_faces(u, v)
#     for face in faces:
#         if face:
#             if form.face_attribute(face, '_is_loaded') is False:
#                 form.edge_attribute((u, v), 'q', q * 10 * 15)

plot_form(form, show_q=False, max_width=3.0).show()

# vault = Shape.from_formdiagram_and_attributes(form)  TODO: fix this
vault = Shape()

vault.data['thk'] = 2*offset
vault.data['t'] = 0.0

# Setting optimisation data

obj = 'bestfit'
solver = 'SLSQP'
constraints = ['envelopexy']
variables = ['q']
features = []  # ['sym']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'current'
gradients = True

# Apply Selfweight and Envelope

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q

apply_envelope_from_shape(form, vault)
apply_selfweight_from_shape(form, vault)

apply_bounds_on_q(form, qmax=0.0)

form_base = form.copy()

form_base = form.copy()

# equilibrium_fdm(form)
# plot_superimposed_diagrams(form, form_base, max_width=3.0).show()
# view_thrust(form).show()

# form = FormDiagram.from_json('/Users/mricardo/compas_dev/me/freeform/tom/form_solution_iter200.json')

# Creating Analysis and Running

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['axis_symmetry'] = axis_sym
optimiser.settings['objective'] = obj
optimiser.settings['plot'] = True
optimiser.settings['printout'] = True
optimiser.settings['max_iter'] = 500
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['derivative_test'] = False
optimiser.settings['save_iterations'] = True

optimiser.settings['starting_point'] = starting_point

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.set_up_optimiser()
analysis.run()

form.to_json(solution)

# view_bestfit_solution(form).show()
view_thrust(form).show()
