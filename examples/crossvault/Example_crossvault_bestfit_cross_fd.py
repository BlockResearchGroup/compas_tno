import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers import view_meshes
import os

# ----------------------------------------------------------------------------
# ---------- EXAMPLE OF MIN THRUST FOR CROSSVAULT WITH CROSS FD --------------
# ----------------------------------------------------------------------------

# Basic parameters

thk = 0.50
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
gradients = False

# ----------------------- 1. Create CrossVault shape ---------------------------

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
# plot_form(form, show_q=False, fix_width=True).show()

# --------------------- 3. Create Initial point with TNA ---------------------

form = form.form_update_with_parallelisation(plot=False)
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular']
optimiser.settings['variables'] = ['ind']
optimiser.settings['objective'] = 'target'
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 900.0
print(optimiser.settings)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# OPT SOL: 13.040

# print('Print target Result - only ind variables')
plot_form(form, show_q=False, cracks=True).show()
# file_address = os.path.join(compas_tno.get(''),'test1.json')
# form.to_json(file_address)

form1 = form.copy()


# --------------------- 6. Run it for MAX ---------------------

optimiser.settings['variables'] = ['ind', 'zb']
analysis.set_up_optimiser()
analysis.run()

# OPT SOL: 6.314

# print('Print target Result - only ind and zb variables')
plot_form(form, show_q=False, cracks=True).show()
# file_address = os.path.join(compas_tno.get(''),'test2.json')
# form.to_json(file_address)

# form, force = form.reciprocal_from_form(plot=True)

view_meshes([form, form1]).show()
