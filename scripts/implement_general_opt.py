from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution
from compas_tno.problems import initialise_problem
from compas_tno.algorithms import zlq_from_q
from compas_tno.algorithms import xyz_from_q
from compas_tno.problems import f_min_thrust_general
from compas_tno.problems import constr_wrapper_general
from compas_tno.problems import sensitivities_wrapper_general
from compas_tno.problems import gradient_fmin_general

from compas_plotters import MeshPlotter

from compas_tno.algorithms import apply_sag

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis

from numpy import zeros
from numpy import ones

from numpy import hstack
from numpy import vstack


span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'cross_fd'
type_structure = 'crossvault'
thk = 0.15
discretisation_shape = 4 * discretisation

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

# Apply Selfweight and Envelope

form.envelope_from_shape(vault)
form.selfweight_from_shape(vault)

form.envelope_on_x_y(c=0.50)

form_base = form.copy()

apply_sag(form)

# form.initialise_loadpath()

# ------------------------------------------------------------
# ------------------- Proper Implementation ------------------
# ------------------------------------------------------------

optimiser = Optimiser()
optimiser.data['library'] = 'SLSQP'
optimiser.data['solver'] = 'SLSQP'
# optimiser.data['library'] = 'IPOPT'
# optimiser.data['solver'] = 'IPOPT'
optimiser.data['constraints'] = ['funicular', 'envelope']
# optimiser.data['variables'] = ['q', 'sym']
optimiser.data['variables'] = ['q', 'sym', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = False
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = True
optimiser.data['jacobian'] = True
optimiser.data['derivative_test'] = True

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

weight = 0
for key in form.vertices():
    weight += form.vertex_attribute(key, 'pz')

thrust = abs(optimiser.fopt)
print('Ratio Thrust/Weight:', thrust/weight)

# print(max(zfree), min(zfree))

for key in form.vertices_where({'is_fixed': True}):
    print(form.vertex_coordinates(key))

# Viewing
plot_superimposed_diagrams(form, form_base).show()
plot_form(form, show_q=False, cracks=True).show()
view_solution(form).show()