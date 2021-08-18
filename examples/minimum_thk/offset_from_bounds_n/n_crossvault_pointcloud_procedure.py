import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_solution
from compas_tno.shapes import MeshDos
from compas.datastructures import mesh_delete_duplicate_vertices
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv_row

from scipy import rand

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
error = 0.0
span = 10.0
k = 1.0
n = 2
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
file_name = 'pointcloud_from_shape'
discretisation = 10
n_step = 0.02
gradients = True  # False

# ----------------------- Shape Analytical ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
}

analytical_shape = Shape.from_library(data_shape)

area_analytical = analytical_shape.middle.area()
swt_analytical = analytical_shape.compute_selfweight()

print('Analytical Self-weight is:', swt_analytical)
print('Analytical Area is:', area_analytical)

# ----------------------- Point Cloud -----------------------

xy = []
points_ub = []
points_lb = []

for i in range(n * discretisation + 1):
    for j in range(n * discretisation + 1):
        xy.append([i * span / (n * discretisation), j * span / (n * discretisation)])

z_ub = analytical_shape.get_ub_pattern(xy).reshape(-1, 1) + error * (2 * rand(len(xy), 1) - 1)
z_lb = analytical_shape.get_lb_pattern(xy).reshape(-1, 1) + error * (2 * rand(len(xy), 1) - 1)

for i in range(len(xy)):
    points_lb.append([xy[i][0], xy[i][1], float(z_lb[i])])
    points_ub.append([xy[i][0], xy[i][1], float(z_ub[i])])

# triangulated_shape = Shape.from_pointcloud(points_lb, points_ub)
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
plot_form(form, show_q=False, fix_width=False).show()

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
# vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub)
# more improved, considers the real middle
vault = Shape.from_pointcloud_and_formdiagram(form, points_lb, points_ub, middle=analytical_shape.middle, data={'type': 'general', 't': 0.0, 'thk': thk})
vault.store_normals()

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

# view_shapes(vault).show()

form.selfweight_from_shape(vault)
# form.selfweight_from_shape(analytical_shape)

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb', 'n']
optimiser.settings['objective'] = 'n'
optimiser.settings['printout'] = False
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients

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
