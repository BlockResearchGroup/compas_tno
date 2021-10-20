from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.viewers import Viewer

from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

span = 10.0
k = 1.0
lamd = None
discretisation = 14
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 1 * discretisation

type_ = None
pattern_load = None  # '/Users/mricardo/compas_dev/me/loadpath/corner/topology/' + type_ + '_complete.json'

c = 0.1

lambd = 0.1

objective = 't'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'
gradients = True

if objective == 't':
    variables.append(objective)
if objective == 'lambd':
    variables.append(objective)

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)

# plot_form(form, show_q=False).show()

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'xy_span': [[0, span], [0, k*span]],
    # 'hc': hc,
    # 'hm': None,
    # 'he': None,
    'center': [5.0, 5.0],
    'radius': span/2,
    't': 0.0,
}

vault = Shape.from_library(data_shape)

# ------------------------------------------------------------
# -----------------------  INITIALISE   ----------------------
# ------------------------------------------------------------

if 'lambd' in variables:
    apply_horizontal_multiplier(form, lambd=lambd)

if 'envelopexy' in constraints:
    apply_envelope_on_xy(form, c=c)

form_base = form.copy()

# ------------------------------------------------------------
# ------------------- Proper Implementation ------------------
# ------------------------------------------------------------

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['features'] = features
optimiser.settings['variables'] = variables
optimiser.settings['objective'] = objective
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = False
optimiser.settings['printout'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
optimiser.settings['derivative_test'] = False
optimiser.settings['max_iter'] = 500

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

data = optimiser.to_data()
optimiser = Optimiser.from_data(data)

if objective in ['min', 'max']:
    weight = 0
    for key in form.vertices():
        weight += form.vertex_attribute(key, 'pz')

    thrust = abs(optimiser.fopt)
    print('Ratio Thrust/Weight:', thrust/weight)

import compas_tno
json_form = compas_tno.get('form.json')
json_shape = compas_tno.get('shape.json')

print('Saving problem data to:')
print(json_form)
print(json_shape)

form.to_json(json_form)
vault.to_json(json_shape)

# Viewing
plot_form(form, show_q=False, cracks=True).show()
view = Viewer(form)
view.show_solution()
