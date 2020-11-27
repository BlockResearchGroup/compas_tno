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
import json
from scipy import rand
from compas_tno.viewers import view_mesh
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv_row

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
type_formdiagram = 'cross_fd'
discretisation = 14
gradients = True  # False
n_step = 0.01

# ----------------------- Point Cloud -----------------------

# file_name = 'fanvaulting_t=50'
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

# THIS IS USEFUL ONLY TO VISUALISE THE INITIAL POINT CLOUD AND DETECT PROBLEMS
# triangulated_shape = Shape.from_pointcloud(points_lb, points_ub)
# view_shapes_pointcloud(triangulated_shape).show()

# ----------------------- Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span_x], [0, span_y]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
plot_form(form, show_q=False, fix_width=False).show()

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
# vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub)
# more improved, considers the real middle
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
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'SLSQP'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb', 'n']
optimiser.data['objective'] = 'n'
optimiser.data['printout'] = False
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients

# --------------------- 5. Set up and run analysis ---------------------

folder = os.path.join('/Users/mricardo/compas_dev/me', 'max_n', file_name, type_structure, type_formdiagram, 'min_max')
os.makedirs(folder, exist_ok=True)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation) + '_offset-method'

forms_address = os.path.join(folder, title)

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.max_n_minmax_GSF(n_step=n_step, save_forms=forms_address, plot=False)
thicknesses, solutions = results

# ----------------------- Save output data --------------------------

csv_file = os.path.join(folder, title + '_data.csv')
save_csv_row(thicknesses, solutions, path=csv_file, title=title, limit_state=False)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True, limit_state=False).show()
