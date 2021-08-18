import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust

# ----------------------------------------------------------------------
# ---- ATTENTION: THIS SCRIPT REQUIRES IPOPT OPTIMISATION PACKAGE ------
# ----------------------------------------------------------------------
# ------ EXAMPLE OF MIN/MAX THRUST FOR ARCH WITH REDUCED THK -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10

# ----------------------- 1. Create CrossVault Shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation,
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Crossvault created!')
# view_shapes(vault).show()

# ------------------- 2. Create and initiate Form Diagram with TNA ------------------------

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

form = form.initialise_tna(plot=False)

# --------------------- 3.1 Create Minimisation for minimum thrust ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'IPOPT'
optimiser.settings['solver'] = 'IPOPT'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['objective'] = 'min'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
print(optimiser.settings)

# --------------------------- 3.2 Run optimisation with ipopt ---------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

import time
start_time = time.time()

analysis.run()

elapsed_time = time.time() - start_time
print('Elapsed Time: {0:.1f} sec'.format(elapsed_time))

# --------------------------- 4. Save and print ---------------------------

file_address = compas_tno.get('test.json')
form.to_json(file_address)
plot_form(form, show_q=False)
view_thrust(form).show()
