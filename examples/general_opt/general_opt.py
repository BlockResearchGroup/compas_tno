from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis


# added a silly change


span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 1 * discretisation

type_ = None
pattern_load = None  # '/Users/mricardo/compas_dev/me/loadpath/corner/topology/' + type_ + '_complete.json'

c = 0.1

objective = 'min'
solver = 'SLSQP'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'
gradients = True

if objective == ['t']:
    variables.append(objective[0])
if objective == ['lambd']:
    variables.append(objective[0])

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)

plot_form(form, show_q=False).show()

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

# Apply Selfweight and Envelope

apply_envelope_from_shape(form, vault)
apply_selfweight_from_shape(form, vault)

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
# optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
optimiser.settings['derivative_test'] = False
optimiser.settings['max_iter'] = 500

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

print(optimiser.settings)

weight = 0
for key in form.vertices():
    weight += form.vertex_attribute(key, 'pz')

thrust = abs(optimiser.fopt)
print('Ratio Thrust/Weight:', thrust/weight)

# Viewing
plot_superimposed_diagrams(form, form_base).show()
plot_form(form, show_q=False, cracks=True).show()
view_solution(form).show()
