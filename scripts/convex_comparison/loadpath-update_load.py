from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.utilities import apply_bounds_on_q
from compas_tno.algorithms import weights_from_xyz
from numpy import array

span = 10.0
k = 1.0

discretisation = [20, 16]
type_formdiagram = 'radial_fd'
type_structure = 'dome'
discretisation_shape = [2 * discretisation[0], 2 * discretisation[1]]

# discretisation = 10
# type_formdiagram = 'cross_fd'
# type_structure = 'crossvault'
# discretisation_shape = 10 * discretisation

save = True
solutions = {}

obj = 'loadpath'
solver = 'MATLAB'
constraints = ['funicular']
variables = ['q', 'zb']
features = ['fixed']
starting_point = 'current'

thk = 0.50  # thickness of the problem
qmax = 1000.0

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': span/2,
    'discretisation': discretisation,
    'xy_span': [[0, span], [0, k*span]],
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'xy_span': [[0, span], [0, k*span]],
    'center': [5.0, 5.0],
    'radius': span/2,
    't': 1.0,
}

vault = Shape.from_library(data_shape)
vault.ro = 20.0

# ------------------------------------------------------------
# -----------------------  INITIALISE   ----------------------
# ------------------------------------------------------------

apply_bounds_on_q(form, qmin=-qmax)

form_base = form.copy()

# ------------------------------------------------------------
# ------------------- Proper Implementation ------------------
# ------------------------------------------------------------

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['objective'] = obj
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = False
optimiser.settings['starting_point'] = starting_point
optimiser.settings['save_iterations'] = True
optimiser.settings['normalize_loads'] = True

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()

analysis.run()

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()

pzt0 = form.lumped_swt()
print('Lumped Load:', pzt0)

# --- Iteration 01

F, V0, V1, V2 = form.tributary_matrices()
xyz = array(form.vertices_attributes('xyz'))
pz = -1 * weights_from_xyz(xyz, F, V0, V1, V2, thk=thk, density=vault.ro)
for i, key in enumerate(form.vertices()):
    form.vertex_attribute(key, 'pz', pz[i])

pzt = form.lumped_swt()
print('New Lumped Load: {0} | Difference (%): {1}'.format(pzt, (pzt - pzt0)/pzt0*100))

pzt0 = pzt

analysis.set_up_optimiser()
analysis.run()

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()

# --- Iteration 02

xyz = array(form.vertices_attributes('xyz'))
pz = -1 * weights_from_xyz(xyz, F, V0, V1, V2, thk=thk, density=vault.ro)
for i, key in enumerate(form.vertices()):
    form.vertex_attribute(key, 'pz', pz[i])

pzt = form.lumped_swt()
print('New Lumped Load: {0} | Difference (%): {1}'.format(pzt, (pzt - pzt0)/pzt0*100))

pzt0 = pzt

analysis.set_up_optimiser()
analysis.run()

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()

# --- Iteration 03

xyz = array(form.vertices_attributes('xyz'))
pz = -1 * weights_from_xyz(xyz, F, V0, V1, V2, thk=thk, density=vault.ro)
for i, key in enumerate(form.vertices()):
    form.vertex_attribute(key, 'pz', pz[i])

pzt = form.lumped_swt()
print('New Lumped Load: {0} | Difference (%): {1}'.format(pzt, (pzt - pzt0)/pzt0*100))

pzt0 = pzt

analysis.set_up_optimiser()
analysis.run()

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()
