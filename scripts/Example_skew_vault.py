from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust
from compas_tno.algorithms import z_from_form

# ----------------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR CROSSVAULT --------------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 0.5
type_structure = 'crossvault'
type_formdiagram = 'ortho'
discretisation = [10, 10]

# ----------------------- 1. Create CrossVault shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [10, 10],
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0
}

vault = Shape.from_library(data_shape)
print('Crossvault created!')

# ----------------------- 2. Create Form Diagram ---------------------------

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
# form.plot()

# --------------------- Create Convex Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'MATLAB'
optimiser.settings['solver'] = 'SDPT3'
optimiser.settings['constraints'] = ['funicular']
optimiser.settings['variables'] = ['ind']
optimiser.settings['objective'] = 'loadpath'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 150.0
print(optimiser.settings)

# -------------- Create Analysis Model and Run Convex Opt --------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.set_up_optimiser() # Find independent edges
analysis.run()
# plot_form(form, show_q=False).show()

for key in form.vertices_where({'is_fixed': True}):
    form.vertex_attribute(key, 'z', 2.0)
    break

z_from_form(form)

file_adress = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form.to_json(file_adress)

view_thrust(form).show()
