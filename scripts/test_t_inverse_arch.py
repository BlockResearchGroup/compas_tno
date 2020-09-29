import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.analysis.analysis import Analysis
from compas_tno.plotters import diagram_of_thrust
from compas_tno.plotters import save_csv
from compas_tno.plotters import save_csv_row

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

vault = Shape.from_library(data_shape)
area = vault.middle.area()
swt = vault.compute_selfweight()
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
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'SLSQP'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb', 't']
optimiser.data['objective'] = 't'
optimiser.data['printout'] = False
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 5000.0
optimiser.data['gradient'] = gradients
optimiser.data['jacobian'] = gradients
optimiser.data['thk'] = thk
print(optimiser.data)

# --------------------------- 3.2 Run optimisation with scipy ---------------------------

# analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
# analysis.set_up_optimiser()
# analysis.run()

folder = os.path.join('/Users/mricardo/compas_dev/me', 'shape_comparison', type_structure, type_formdiagram)
os.makedirs(folder, exist_ok=True)
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)
# if sag:
#     title = title + 'sag_' + str(sag)
# if smooth:
#     title = title + 'smooth_'
forms_address = os.path.join(folder, title)

thk_max = 0.30

analysis = Analysis.from_elements(vault, form, optimiser)
results = analysis.thk_minmax_GSF(thk_max, save_forms=forms_address, plot=False)
thicknesses, solutions = results

# ----------------------- Save output data --------------------------

csv_file = os.path.join(folder, title + '_data.csv')
save_csv_row(thicknesses, solutions, path=csv_file, title=title)

img_graph = os.path.join(folder, title + '_diagram.pdf')
diagram_of_thrust(thicknesses, solutions, save=img_graph, fill=True)




form = analysis.form
file_address = compas_tno.get('test.json')
form.to_json(file_address)
optimiser = analysis.optimiser
fopt = optimiser.fopt
save_photo = False
blocks_on_plot = False

thk_min = form.attributes['thk']
print(thk_min)
data_shape['thk'] = thk_min
arch = Shape.from_library(data_shape)
form.envelope_from_shape(arch)

plot_form_xz(form, arch, show_q=False, plot_reactions='simple', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, save=save_photo, hide_cracks=False, hide_negative=True).show()
