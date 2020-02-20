
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters.plotters import plot_form
from compas_tno.analysis.analysis import Analysis

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR DOME --------------------
# ----------------------------------------------------------------------

# Basic parameters

thk = 1.0
radius = 5.0
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'

# ----------------------- Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': 0.5,
    'discretisation': [50, 50],
    'xy_span': [[0.0,10.0],[0.0,10.0]],
    't' : 0.0
}

vault = Shape.from_library(data_shape)
swt = vault.compute_selfweight()

print('Crossvault created!')

# view_shapes(vault).show()

# ----------------------- Create Form Diagram ---------------------------

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
plot_form(form, show_q=False, fix_width=True).show()

# ----------------------- Create Optimiser ----------------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['qmax'] = 100.0
print(optimiser.data)

# ----------------------- Create Analysis Model ------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser() # Find independent edges
# plot_form(form, show_q=False).show()
analysis.run()

plot_form(form, show_q=False).show()
