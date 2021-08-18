
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_normals
from compas_tno.viewers import view_mesh
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_solution
from compas_tno.shapes import MeshDos
from compas.geometry import normalize_vector
from compas.geometry import sum_vectors
from compas.geometry import norm_vector
from compas.geometry import scale_vector
from compas.geometry import angle_vectors
from compas.geometry import centroid_points
import math
from copy import deepcopy

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------


thk = 1.0
A = 1.05  # Parameter similar to a ""
span_x = 10.0
span_y = 10.0
span = max(span_x, span_y)
xy_span = [[0, span_x], [0, span_y]]
xy_span_shape = [[-span_x/2*(A - 1), span_x*(1 + (A - 1)/2)], [-span_y/2*(A - 1), span_y*(1 + (A - 1)/2)]]
k = 1.0
n = 4
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
gradients = True

# --------------------
# Shape
# --------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': xy_span_shape,
    't': 0.0,
}

analytical_shape = Shape.from_library(data_shape)
swt = analytical_shape.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')

# --------------------
# Diagrams
# --------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': xy_span,
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')

# --------------------
# PointCloud
# --------------------

xy = []
points_central = []

for i in range(n * discretisation + 1):
    for j in range(n * discretisation + 1):
        xy.append([i * span / (n * discretisation), j * span / (n * discretisation)])

zt = analytical_shape.get_middle_pattern(xy).reshape(-1, 1)

for i in range(len(xy)):
    points_central.append([xy[i][0], xy[i][1], float(zt[i])])

vault = Shape.from_middle_pointcloud(points_central, topology=form, thk=thk, treat_creases=True)

# --------------------
# Initialise
# --------------------

# form = form.initialise_tna(plot=False)
form.selfweight_from_shape(vault)
form.envelope_from_shape(vault)
form.initialise_loadpath()
# plot_form(form).show()

# --------------------
# Optimiser
# --------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb', 't']
optimiser.settings['objective'] = 't'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
print(optimiser.settings)

# --------------------
# Analysis
# --------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()
analysis.run()

plot_form(form, simple=True, cracks=True).show()

view_solution(form, vault).show()

# view_shapes(vault).show()
# view_mesh(vault.intrados, normals=True).show()

# vault.store_normals()
# view_normals(vault).show()
