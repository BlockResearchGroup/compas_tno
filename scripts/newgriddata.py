import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from copy import deepcopy

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8, 20]
n = 1

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [n*discretisation[0], n*discretisation[1]],
    'center': [5.0, 5.0],
    'radius': radius,
    't' : 1.0
}

dome = Shape.from_library(data_shape)
swt = dome.compute_selfweight()
print('Selfweight computed:', swt)
print('Vault geometry created!')

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': radius,
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)


form.selfweight_from_shape(dome)


# plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form = form.initialise_loadpath()
# plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = 'Scipy'
optimiser.settings['solver'] = 'slsqp'
optimiser.settings['constraints'] = ['funicular', 'envelope', 'reac_bounds', ' symmetry']
optimiser.settings['variables'] = ['ind', 'zb']
optimiser.settings['objective'] = 'min'
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['qmax'] = 3000.0
print(optimiser.settings)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

from numpy import array
XY = form.vertices_attributes('xy')
# print(XY[0][1])
# print(XY.shape)

# zt = dome.get_middle_pattern(XY)
# zub = dome.get_ub_pattern(XY)
# zlb = dome.get_lb_pattern(XY)

# print(zt)
# print(zub)
# print(zlb)

# print(max(zub), min(zlb))

# form.vertices_attributes('target', zt)
# form.vertices_attributes('ub', zub)
# form.vertices_attributes('lb', zlb)

analysis.run()
