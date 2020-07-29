import os
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form_xz
from compas_tno.plotters import plot_forms_xz
from compas_tno.analysis.analysis import Analysis
from copy import deepcopy

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR ARCH --------------------
# ----------------------------------------------------------------------

# Basic parameters

H = 1.0
L = 2.0
thk = 0.2 # 0.108
discretisation = 20
b = 0.5
t = 10.0
type_structure = 'arch'
type_formdiagram = 'arch'
forms = []

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
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 1000.0
print(optimiser.data)

# --------------------------- 3.2 Run optimisation with scipy ---------------------------

analysis = Analysis.from_elements(arch, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()

analysis.set_up_optimiser()
analysis.run()

forms.append(deepcopy(form))

# save_photo =  os.path.join(compas_tno.get('/imgs/'), 'arch_minthk_noblocks_' + optimiser.data['objective'] + '.pdf')
blocks_on_plot = False
# plot_form_xz(form, arch, show_q=False, plot_reactions='simple', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, save=save_photo, hide_negative=True).show()

# --------------------- 4.1 Modify Minimisation for maximum thrust ---------------------

# optimiser.data['objective'] = 'bestfit'
# analysis.set_up_optimiser()
# analysis.run()
# forms.append(deepcopy(form))


# optimiser.data['objective'] = 'max'
# optimiser.data['qmax'] = 102.6
# analysis.set_up_optimiser()
# analysis.run()
# forms.append(deepcopy(form))

# optimiser.data['qmax'] = 102.6
# Constraint [ 87.6 / 102.6 / 117.6 ]

# optimiser.data['objective'] = 'min'
# for key in form.vertices():
#     form.vertex_attribute(key, 'ub', 0.97)
# analysis.set_up_optimiser()
# analysis.run()
# for key in form.vertices():
#     form.vertex_attribute(key, 'ub', 1.0 + thk/2)
# forms.append(deepcopy(form))

# optimiser.data['objective'] = 'min'
# analysis.set_up_optimiser()
# analysis.run()
# forms.append(deepcopy(form))

optimiser.data['objective'] = 'max'
analysis.set_up_optimiser()
analysis.run()
forms.append(deepcopy(form))

colours = ['0000FF', 'FF0000']

save_photo = os.path.join(compas_tno.get('/imgs/'), 'arch_noblocks_minmax_' + optimiser.data['objective'] + '.pdf')
# plot_form_xz(form, arch, show_q=False, plot_reactions='simple', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, save=save_photo, hide_negative=True).show()

# print(forms)

# save_photo = os.path.join(compas_tno.get('/imgs/'), 'arch_various_.pdf')
plot_forms_xz(forms, arch, colours=colours, plot_reactions='simple', fix_width=True, max_width=5, radius=0.02, stereotomy=blocks_on_plot, hide_negative=True, hide_cracks=True, save=save_photo).show()



# optimiser.data['qmax'] = 102.6
# Constraint [ 87.6 / 102.6 / 117.6 ]
