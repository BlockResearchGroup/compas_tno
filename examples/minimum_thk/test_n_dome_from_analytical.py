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
from compas_tno.datastructures import MeshDos
import math
from scipy import rand

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome_polar'
type_formdiagram = 'radial_fd'
discretisation = [8, 20]
center = [5.0, 5.0]
gradients = False
n = 1
error = 0.00
ro = 20.0

# ----------------------- Shape Analytical ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*n, discretisation[1]*n],
    'center': center,
    'radius': radius,
    't': 0.0,
    'expanded': True
}

analytical_shape = Shape.from_library(data_shape)
analytical_shape.ro = ro
area_analytical = analytical_shape.middle.area()
swt_analytical = analytical_shape.compute_selfweight()

print('Analytical Self-weight is:', swt_analytical)
print('Analytical Area is:', area_analytical)

analytical_shape.store_normals()
# view_shapes(analytical_shape).show()

# ----------------------- Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
# plot_form(form, show_q=False, fix_width=False).show()

# ------- Create shape given a topology and a shape defined by meshes --------

vault = Shape.from_meshes_and_formdiagram(form, analytical_shape.intrados, analytical_shape.extrados, middle=analytical_shape.middle, data={'type': 'general', 't': 0.0, 'thk': thk})

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is: {0:.2f} diff ({1:.2f}%)'.format(swt, 100*(swt - swt_analytical)/(swt_analytical)))
print('Area is: {0:.2f} diff ({1:.2f}%)'.format(area, 100*(area - area_analytical)/(area_analytical)))

form.selfweight_from_shape(vault)
# view_shapes(vault).show()

# Plot the normals:

normals = {}
for key in vault.intrados.vertices():
    normals[key] = vault.intrados.vertex_attribute(key, 'n')
    print(vault.intrados.vertex_attribute(key, 'n'))

from compas_plotters import MeshPlotter
from compas.utilities import i_to_red
plotter = MeshPlotter(vault.intrados, figsize=(10, 10))
plotter.draw_edges()
plotter.draw_vertices(text=normals)
plotter.show()

view_shapes(vault).show()
# view_shapes_pointcloud(vault).show()

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'SLSQP'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
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

n_reduction = -1 * analysis.optimiser.fopt
thk_min = thk - 2*n_reduction*thk
print('Approx. Minimum THK:', thk_min)
data_shape['thk'] = thk_min

plot_form(form, show_q=False, cracks=True).show()

# analytical_shape = Shape.from_library(data_shape)
# form.envelope_from_shape(analytical_shape)

# plot_form(form, show_q=False, cracks=True).show()

form.to_json(compas_tno.get('test.json'))

view_solution(form, vault).show()
