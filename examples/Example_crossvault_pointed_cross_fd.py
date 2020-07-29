import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust
from compas_tno.viewers.thrust import view_solution
import math

# ----------------------------------------------------------------------------
# ---------- EXAMPLE OF MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 0.5
type_structure = 'pointed_crossvault'
type_formdiagram = 'fan_fd'
discretisation = 10
hc = 5 * math.sqrt(2)

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0,
    'hc': hc,
    'hm': None,
    'he': None,
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Crossvault created!')
# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 3. Create Initial point with TNA ---------------------

form = form.initialise_tna(plot=False)
plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'max'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 900.0
print(optimiser.data)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

plot_form(form, show_q=False, simple=True, cracks=True).show()

import os
file_address = os.path.join(compas_tno.get('/rqe/'), type_structure + '_' + type_formdiagram + '_t=50_'+ optimiser.data['objective'] + '.json')
form.to_json(file_address)

view_solution(form, vault).show()
