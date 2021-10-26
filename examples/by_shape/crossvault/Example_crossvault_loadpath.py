from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
import compas_tno

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN THRUST FOR CROSS FORM DIAGRAM ----------------
# ----------------------------------------------------------------------

# Basic parameters

thk = 1.0
radius = 5.0
type_structure = 'crossvault'
type_formdiagram = 'fan_fd'

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': 1.0,
    'discretisation': [10, 10],
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Crossvault created!')

# view_shapes(vault).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0,10],[0,10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)
# plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 3. Create Convex Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'MATLAB'
optimiser.settings['solver'] = 'SDPT3'
optimiser.settings['constraints'] = ['funicular']
optimiser.settings['variables'] = ['ind']
optimiser.settings['objective'] = 'loadpath'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 3000.0
print(optimiser.settings)

# -------------- 4. Create Analysis Model and Run Convex Opt --------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.set_up_optimiser() # Find independent edges
analysis.run()
plot_form(form, show_q=False).show()

import os
file_address = os.path.join(compas_tno.get('/rqe/'), type_structure + '_' + type_formdiagram + '_t=50_'+ optimiser.settings['objective'] + '.json')
form.to_json(file_address)

# file_address = compas_tno.get('test.json')
# form.to_json(file_address)
