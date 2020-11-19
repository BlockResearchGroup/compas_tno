import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_normals
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_mesh
from compas_tno.datastructures import MeshDos
import json
from scipy import rand
from numpy import array

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.50
error = 0.0
span = 10.0
k = 1.0
n = 2
type_structure = 'crossvault'
type_formdiagram = 'fan_fd'
discretisation = 10
gradients = True  # False

# ----------------------- Point Cloud -----------------------

file_name = 'jeronimos2'
pointcloud = '/Users/mricardo/compas_dev/me/min_thk/pointcloud/' + file_name + '.json'

middle = []
tol = 10e-4

with open(pointcloud) as json_file:
    data = json.load(json_file)
    for key, pt in data['target'].items():
        middle.append(pt)

triangulated_shape = Shape.from_middle_pointcloud(middle, thk=thk)
# view_normals(triangulated_shape).show()
# view_shapes_pointcloud(triangulated_shape).show()


# ----------------------- Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
# plot_form(form, show_q=False, fix_width=False).show()

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
# vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub)
# more improved, considers the real middle
points_lb = triangulated_shape.intrados.vertices_attributes('xyz')
points_ub = triangulated_shape.extrados.vertices_attributes('xyz')

i = 0
for pt in points_lb:
    if abs(pt[0] - 10) < 10e-3 and abs(pt[1] - 10) < 10e-3:
        points_lb[i][:2] = [10.0, 10.0]
    if abs(pt[0] - 8) < 10e-3 and abs(pt[1] - 10) < 10e-3:
        points_lb[i][:2] = [8.0, 10.0]
    i += 1
i = 0
for pt in points_ub:
    if abs(pt[0] - 10) < 10e-3 and abs(pt[1] - 10) < 10e-3:
        points_ub[i][:2] = [10.0, 10.0]
    if abs(pt[0] - 8) < 10e-3 and abs(pt[1] - 10) < 10e-3:
        points_ub[i][:2] = [8.0, 10.0]
    i += 1

vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub, data={'type': 'general', 't': 0.0, 'thk': thk})
vault.store_normals()
# view_shapes_pointcloud(vault).show()
# view_normals(vault).show()

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is: {0:.2f}'.format(swt))
print('Area is: {0:.2f}'.format(area))

# view_shapes(vault).show()

form.selfweight_from_shape(vault)
# form.selfweight_from_shape(analytical_shape)

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.initialise_loadpath()
plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'SLSQP'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb', 'n']
optimiser.data['objective'] = 'n'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

n_reduction = - 1 * analysis.optimiser.fopt
thk_min = thk - 2 * n_reduction
print('Approx. Minimum THK:', thk_min)
# data_shape['thk'] = thk_min

plot_form(form, show_q=False, cracks=True).show()

# analytical_shape = Shape.from_library(data_shape)
# form.envelope_from_shape(analytical_shape)

# plot_form(form, show_q=False, cracks=True).show()

form.to_json(compas_tno.get('test.json'))

view_solution(form, vault).show()
