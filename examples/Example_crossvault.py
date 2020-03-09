from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis

# ----------------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR CROSSVAULT --------------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 1.0
radius = 5.0
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'

# ----------------------- Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': 0.5,
    'discretisation': [10, 10],
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Crossvault created!')
# view_shapes(vault).show()

# ----------------------- Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, 10], [0, 10]],
    'discretisation': 10,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
form.overview_forces()
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=True).show()

# --------------------- Create Convex Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'MATLAB'
optimiser.data['solver'] = 'SDPT3'
optimiser.data['constraints'] = ['funicular']
optimiser.data['variables'] = ['ind']
optimiser.data['objective'] = 'loadpath'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 900.0
print(optimiser.data)

# -------------- Create Analysis Model and Run Convex Opt --------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.set_up_optimiser()  # Find independent edges
analysis.run()
plot_form(form, show_q=False).show()

# --------------------- Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 900.0
print(optimiser.data)

# --------------------- Create Minimisation Optimiser ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()  # Find independent edges
analysis.run()

plot_form(form, show_q=False).show()

file_address = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form.to_json(file_address)
