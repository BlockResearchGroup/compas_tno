from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.utilities import apply_bounds_on_q
import compas_tno

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

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

analysis.run()

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()

q_matlab = optimiser.M.q[:]

# ------------

optimiser.settings['solver'] = 'CVXPY'

analysis.set_up_optimiser()
analysis.run()

q_cvxpy = optimiser.M.q[:]

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()


# ------------

optimiser.settings['find_inds'] = True
optimiser.settings['solver'] = 'MATLAB'

analysis.set_up_optimiser()
analysis.run()

q_matlab_ind = optimiser.M.q[:]

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()

# ------------

optimiser.settings['find_inds'] = True
optimiser.settings['solver'] = 'CVXPY'

analysis.set_up_optimiser()
analysis.run()

q_cvxpy_ind = optimiser.M.q[:]

view = Viewer(form, vault)
view.draw_thrust()
view.draw_force()
view.show()

# # ---- MAKE THE VIDEO ----

# from compas_tno.viewers import animation_from_optimisation
# from compas_tno.algorithms import reciprocal_from_form

# DATA_XFORM = compas_tno.get('Xform.json')
# DATA_XFORCE = compas_tno.get('Xforce.json')

# force = reciprocal_from_form(form)

# animation_from_optimisation(form, DATA_XFORM, force, DATA_XFORCE, interval=150)
