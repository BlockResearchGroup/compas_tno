import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_shapes
from compas_tno.viewers import Viewer

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR ARCH --------------------
# ----------------------------------------------------------------------

# Basic parameters

incl = 1.0
H = 1.0*incl
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
    't': 0.0,
    'x0': 0.0
}

arch = Shape.from_library(data_shape)

data_general = {
    'type': 'general',
    't': 0.0,
}

vault = Shape.from_meshes(arch.intrados, arch.extrados, middle=None, data=data_general)

area = vault.middle.area()
swt = vault.compute_selfweight()

area_arch = arch.middle.area()
swt_arch = arch.compute_selfweight()

print('Arch created!')

print('Self-weight is:', swt)
print('Area is:', area)

print('Self-weight is:', swt_arch)
print('Area is:', area_arch)

# view_shapes(vault).show()

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
optimiser.settings['constraints'] = ['funicular', 'envelope']
optimiser.settings['variables'] = ['ind', 'zb', 's']
optimiser.settings['objective'] = 's'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['thk'] = thk
print(optimiser.settings)

# --------------------------- 3.2 Run optimisation with scipy ---------------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

form = analysis.form
file_address = compas_tno.get('test.json')
form.to_json(file_address)
optimiser = analysis.optimiser
fopt = -1 * optimiser.fopt  # inverse because we are min inv
save_photo = False
blocks_on_plot = False

vault.update_dos_from_form(form)
view_solution(form, vault).show()

plot_form(form).show()

# view_shapes(vault).show()

thk_min = thk - 2*fopt*thk
print('Minimum THK:', thk_min)
data_shape['thk'] = thk_min
arch = Shape.from_library(data_shape)
form.envelope_from_shape(arch)

plot_form_xz(form, arch, show_q=False, plot_reactions='all', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, save=save_photo, hide_negative=True).show()
