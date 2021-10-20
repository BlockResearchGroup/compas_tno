import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
import json


# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.44
error = 0.0
span_x = 5.97447  # 5.76 - 7.2 - 5.97448
span_y = 11.2  # 9.85 - 12.31 - 11.2
k = 1.0
n = 2
type_structure = 'crossvault'
type_formdiagram = 'fan_fd'
discretisation = 14

objective = 'n'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['ind', 'zb', 'n']
features = ['fixed']  # , 'symmetry']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'current'
gradients = True
max_iter = 500

# ----------------------- Point Cloud -----------------------

file_name = 'fanvaulting_t=50'
file_name = 'amiens_internet'
pointcloud = '/Users/mricardo/compas_dev/me/min_thk/pointcloud/' + file_name + '.json'

points_ub = []
points_lb = []

tol = 10e-4

with open(pointcloud) as json_file:
    data = json.load(json_file)
    for key, pt in data['UB'].items():
        points_ub.append(pt)
    for key, pt in data['LB'].items():
        points_lb.append(pt)

# triangulated_shape = Shape.from_pointcloud(points_lb, points_ub)
# view_shapes_pointcloud(triangulated_shape).show()

data_diagram_base = {
    'type': 'cross_fd',
    'xy_span': [[0, span_x], [0, span_y]],
    'discretisation': discretisation*2,
    'fix': 'corners',
}
form_base = FormDiagram.from_library(data_diagram_base)

# structured_shape = Shape.from_pointcloud_and_topology(form_base, points_lb, points_ub)
# view_shapes_pointcloud(structured_shape).show()
# view_meshes([structured_shape.intrados, structured_shape.extrados]).show()

# ----------------------- Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span_x], [0, span_y]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
# plot_form(form, show_q=False, fix_width=False).show()

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
vault = Shape.from_pointcloud_and_topology(form, points_lb, points_ub, data={'type': 'general', 't': 0.0, 'thk': thk})
# more improved, considers the real middle
# # vault = Shape.from_meshes_and_formdiagram(form, structured_shape.intrados, structured_shape.extrados, data={'type': 'general', 't': 0.0, 'thk': thk})
vault.store_normals(plot=False)
# view_shapes_pointcloud(vault).show()
# view_normals(vault).show()

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is: {0:.2f}'.format(swt))
print('Area is: {0:.2f}'.format(area))

# view_shapes(vault).show()

from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_from_shape

apply_selfweight_from_shape(form, vault)
apply_envelope_from_shape(form, vault)

# --------------------- 3. Create Starting point with TNA ---------------------

from compas_tno.problems import initialize_loadpath
from compas_tno.algorithms import form_update_with_parallelisation

form_update_with_parallelisation(form, plot=True)
# initialize_loadpath(form)
# form = form.form_update_with_parallelisation(plot=False)
# form.initialise_loadpath()
print('back_here')
plot_form(form, show_q=True).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['objective'] = objective
optimiser.settings['features'] = features
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['max_iter'] = max_iter

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
