
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from copy import deepcopy

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.50
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [4, 12]
ro = 1.0
type_diagonal = 'straight'
gradients = True
n = 1

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*n, discretisation[1]*n],
    'center': [5.0, 5.0],
    'radius': radius,
    't': 0.0
}

dome = Shape.from_library(data_shape)
dome.ro = ro
swt = dome.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': type_diagonal,
}

form = FormDiagram.from_library(data_diagram)

plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Normal Starting point ---------------------

form.selfweight_from_shape(dome)
form.envelope_from_shape(dome)
form.initialise_loadpath()
# form = form.form_update_with_parallelisation(plot=False)
plot_form(form).show()

# Print of SWT at Nodes
from compas_plotters import MeshPlotter
tol = 6 #10e-4
plotter = MeshPlotter(form, figsize=(10, 10))
plotter.draw_edges()
plotter.draw_vertices(radius=0.10)
plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 2) for key in form.vertices() if radius - tol < form.vertex_attribute(key, 'x') < radius + tol}, radius=0.10)
plotter.show()

# Print of UB at Nodes
plotter = MeshPlotter(form, figsize=(10, 10))
plotter.draw_edges()
plotter.draw_vertices(radius=0.10)
plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'ub'), 2) for key in form.vertices() if radius - tol < form.vertex_attribute(key, 'x') < radius + tol}, radius=0.10)
plotter.show()

# Print of LB at Nodes
plotter = MeshPlotter(form, figsize=(10, 10))
plotter.draw_edges()
plotter.draw_vertices(radius=0.10)
plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'lb'), 2) for key in form.vertices() if radius - tol < form.vertex_attribute(key, 'x') < radius + tol}, radius=0.10)
plotter.show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

solvers = [['Scipy', 'SLSQP'], ['MMA', 'MMA'], ['IPOPT', 'IPOPT']]
i = 0

optimiser = Optimiser()
optimiser.settings['library'] = solvers[i][0]
optimiser.settings['solver'] = solvers[i][1]
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds', 'symmetry']
optimiser.settings['variables'] = ['ind', 'zb', 't']
optimiser.settings['objective'] = 't'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
print(optimiser.settings)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
# analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

thk_min = form.attributes['thk']
print(thk_min)
data_shape['thk'] = thk_min
dome = Shape.from_library(data_shape)
form.envelope_from_shape(dome)

plot_form(form, show_q=False, cracks=True).show()

view_solution(form, dome).show()

