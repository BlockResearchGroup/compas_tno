from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust

# ----------------------------------------------------------------------
# ------ ATTENTION: THIS SCRIPT REQUIRES MMA OPTIMISATION PACKAGE ------
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

# form = form.form_update_with_parallelisation(plot=False)

# --------------------- 3. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'MMA'
optimiser.settings['solver'] = 'MMA'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['objective'] = 'min'
optimiser.settings['printout'] = True
optimiser.settings['solver_options']['derivatives'] = 'PyTorch'  # 'DF_brute' 'DF_reduced' and 'analytical' in process.
optimiser.settings['plot'] = True
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
print(optimiser.settings)

# --------------------- 4.2 Create Minimisation Optimiser ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()  # Find independent edges
analysis.run()

form = analysis.form
plot_form(form, show_q=False).show()

file_adress = '/Users/mricardo/compas_dev/me/reformulation/test.json'
# form.to_json(file_adress)

view_thrust(form).show()
