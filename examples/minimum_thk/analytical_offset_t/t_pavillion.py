
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
# ----------- EXAMPLE OF MIN THRUST FOR Vault WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.50
span = 10.0
k = 1.0
n = 4
type_structure = 'pavillionvault'
type_formdiagram = 'ortho_fd'
discretisation = 10
gradients = True

# ----------------------- 1. Create Vault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()
# vault.ro = 1.0
print('Selfweight computed:', swt)

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'all',
}

form = FormDiagram.from_library(data_diagram)
plot_form(form, show_q=False, fix_width=False).show()
print('Form Diagram Created!')
print(form)

# --------------------- 3. Create Starting point with TNA ---------------------

form.selfweight_from_shape(vault)
# form = form.initialise_tna(plot=True)
form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb', 't']
optimiser.settings['objective'] = 't'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 3000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
print(optimiser.settings)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

thk_min = form.attributes['thk']
data_shape['thk'] = thk_min
vault = Shape.from_library(data_shape)
form.envelope_from_shape(vault)
file_address = compas_tno.get('test.json')
form.to_json(file_address)

plot_form(form, show_q=False, cracks=True).show()

view_solution(form, vault).show()
