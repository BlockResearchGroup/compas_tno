
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from compas_tno.viewers.thrust import view_solution
from copy import deepcopy
from compas_plotters import MeshPlotter

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8, 20]
ro = 1.0
gradients = True
qmax = 1000.0
t = 0.0
n = 1
plot_init = False

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation[0]*n, discretisation[1]*n],
    'center': [5.0, 5.0],
    'radius': radius,
    't': t
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
    'partial_diagonal': 'right',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
if plot_init:
    plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

form.selfweight_from_shape(dome)
form.initialise_loadpath()
# form = form.initialise_tna(plot=False)
form.envelope_from_shape(dome)
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

solvers = [['Scipy', 'SLSQP'], ['MMA', 'MMA'], ['IPOPT', 'IPOPT']]
i = 0

optimiser = Optimiser()
optimiser.data['library'] = solvers[i][0]
optimiser.data['solver'] = solvers[i][1]
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb', 't']
optimiser.data['objective'] = 't'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = qmax
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients
print(optimiser.data)

# tol = 10e-3
# # Print of SWT at Nodes
# plotter = MeshPlotter(form, figsize=(10, 10))
# plotter.draw_edges()
# plotter.draw_vertices(radius=0.10)
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 2) for key in form.vertices() if radius - tol < form.vertex_attribute(key, 'x') < radius + tol}, radius=0.10)
# plotter.show()

# # Print of UB at Nodes
# plotter = MeshPlotter(form, figsize=(10, 10))
# plotter.draw_edges()
# plotter.draw_vertices(radius=0.10)
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'ub'), 2) for key in form.vertices() if radius - tol < form.vertex_attribute(key, 'x') < radius + tol}, radius=0.10)
# plotter.show()

# # Print of LB at Nodes
# plotter = MeshPlotter(form, figsize=(10, 10))
# plotter.draw_edges()
# plotter.draw_vertices(radius=0.10)
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'lb'), 2) for key in form.vertices() if radius - tol < form.vertex_attribute(key, 'x') < radius + tol}, radius=0.10)
# plotter.show()

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
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

