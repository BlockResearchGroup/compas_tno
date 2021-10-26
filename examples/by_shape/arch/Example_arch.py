import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis.analysis import Analysis

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR ARCH --------------------
# ----------------------------------------------------------------------

# Basic parameters

H = 1.0
L = 2.0
thk = 0.2
discretisation = 20
b = 0.5  # Out of plane dimension  of arch
t = 1.0
type_structure = 'arch'
type_formdiagram = 'arch'

# ----------------------- 1. Create Arch shape ---------------------------

data_shape = {
    'type': type_structure,
    'H': H,
    'L': L,
    'thk': thk,
    'discretisation': discretisation,
    'b': b,
    't': t,
    'x0': 0.0
}

arch = Shape.from_library(data_shape)
area = arch.middle.area()
swt = arch.compute_selfweight()
print('Arch created!')
print('Self-weight is:', swt)
print('Area is:', area)
# view_shapes(arch).show()

gradients = True

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'H': H,
    'L': L,
    'total_nodes': discretisation,
    'x0': 0.0
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)

# --------------------- 3.1 Create Minimisation for minimum thrust ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'SLSQP'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['objective'] = 'min'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['derivative_test'] = True
print(optimiser.settings)

# --------------------------- 3.2 Run optimisation with scipy ---------------------------

analysis = Analysis.from_elements(arch, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()
form = analysis.form
# file_address = compas_tno.get('test.json')
# form.to_json(file_address)
optimiser = analysis.optimiser
fopt = optimiser.fopt
save_photo = False
blocks_on_plot = False
plot_form_xz(form, arch, show_q=False, plot_reactions='simple', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, save=save_photo, hide_negative=True).show()

# --------------------- 4.1 Modify Minimisation for maximum thrust ---------------------

optimiser.settings['objective'] = 'max'

# ------------------------- 4.2 Run optimisation with scipy ---------------------------

analysis.set_up_optimiser()
analysis.run()
form = analysis.form
file_address = compas_tno.get('test.json')
form.to_json(file_address)
save_photo = False
plot_form_xz(form, arch, show_q=False, plot_reactions='simple', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, save=save_photo, hide_negative=True).show()

