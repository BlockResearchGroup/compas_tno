import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import Viewer
from compas_tno.shapes import MeshDos

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
span = 10.0
k = 1.0
n = 1
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'
discretisation = 10
gradients = True

# ----------------------- Shape Analytical ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation*n,
    'xy_span': [[0, span], [0, k*span]],
    't': 0.0,
}

analytical_shape = Shape.from_library(data_shape)

area_analytical = analytical_shape.middle.area()
swt_analytical = analytical_shape.compute_selfweight()

print('Self-weight is:', swt_analytical)
print('Area is:', area_analytical)

# ----------------------- Form Diagram ----------------------------

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners',
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
# plot_form(form, show_q=False, fix_width=False).show()

form.selfweight_from_shape(analytical_shape)
form.envelope_from_shape(analytical_shape)

# ----------------------- Build Identical Intrados/extrados ---------------------------

intrados = MeshDos.from_formdiagram_attribute(form, attribute='lb')
extrados = MeshDos.from_formdiagram_attribute(form, attribute='ub')

data_general = {
    'type': 'general',
    't': 0.0,
}

vault = Shape.from_meshes(intrados, extrados, middle=None, data=data_general)

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is:', swt)
print('Area is:', area)

# view_shapes(vault).show()

# --------------------- 3. Create Starting point with TNA ---------------------

# form.form_update_with_parallelisation(plot=False)
form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb', 'n']
optimiser.settings['objective'] = 'n'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

n_reduction = -1*analysis.optimiser.fopt

vault.update_dos_from_form(form)
plot_form(form, show_q=False, cracks=True).show()

thk_min = thk - 2*n_reduction*thk
print('Minimum THK:', thk_min)
data_shape['thk'] = thk_min

# analytical_shape = Shape.from_library(data_shape)
# form.envelope_from_shape(analytical_shape)

# plot_form(form, show_q=False, cracks=True).show()

view_solution(form, vault).show()
