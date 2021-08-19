
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers import view_thrusts
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes
from copy import deepcopy
import math

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
span = 10.0
k = 1.0
n = 1
R = 0.70 * span
hc = math.sqrt(R**2 - (R - span/2)**2)
print('Analysis with hc:', hc)
# hc = 6.0  # 5.0*math.sqrt(2)
type_structure = 'pointed_crossvault'
type_formdiagram = 'fan_fd'
discretisation = 14
gradients = True

# ----------------------- 1. Create Vault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 1.0,
    'hc': hc,
    'hm': None,
    'he': None,  # [5.0, 5.0, 5.0, 5.0],
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')
# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
# plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.selfweight_from_shape(vault)

keys_diagonal = [2, 61, 121, 181, 241, 301, 361, 392, 362, 302, 242, 182, 122, 62, 4]
keys_open_edge = [2, 3, 65, 125, 185, 245, 305, 365, 306, 246, 186, 126, 66, 5, 4]

pz_diagonal = []
pz_open_edge = []

for key in keys_diagonal:
    pz = form.vertex_attribute(key, 'pz')
    pz_diagonal.append(pz)
for key in keys_open_edge:
    pz = form.vertex_attribute(key, 'pz')
    pz_open_edge.append(pz)

print('-------- Loads ---------')
print(pz_open_edge)


save_lp = compas_tno.get('lp.json')
try:
    form.from_json(save_lp)
except:
    form = form.initialise_loadpath()
    form.to_json(save_lp)
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

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

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

thk_min = form.attributes['thk']
data_shape['thk'] = thk_min
print('Min THK = ', thk_min)
vault = Shape.from_library(data_shape)
form.envelope_from_shape(vault)
file_address = compas_tno.get('test.json')
form.to_json(file_address)

plot_form(form, show_q=False, cracks=True).show()
view_solution(form, vault).show()
