import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_thrust
from compas_tno.viewers import view_thrusts
from compas_tno.viewers import view_shapes
from copy import deepcopy
import math

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome_spr'
type_formdiagram = 'radial_fd'
discretisation = [8, 16]
center = [5.0, 5.0]
k = 0.66
spr_angle = k * math.pi/2
radius_projected = radius * math.sin(spr_angle)
print('Radius projected:', radius_projected)

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [16, 32],  # Always keep the discretisation higher but also a multiplier. Ex: times 2.0
    'center': center,
    'radius': radius,
    'theta': [0, spr_angle],
    't': 10.0,
}

dome = Shape.from_library(data_shape)
print('Vault geometry created!')
swt = dome.compute_selfweight()
print('Selfweight computed:', swt)
# view_shapes(dome).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': center,
    'radius': radius_projected,  # Attention here, use radius_projected
    'discretisation': discretisation,
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

form = form.initialise_tna(plot=False)
plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 10000.0
print(optimiser.data)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

print('Print min Result')
plot_form(form, show_q=False, cracks=True).show()
form_min = deepcopy(form)

# # --------------------- 6. Run it for MIN ---------------------

# optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
# analysis.set_up_optimiser()
# analysis.run()


# --------------------- 6. Run it for MAX ---------------------

optimiser.data['objective'] = 'max'
analysis.set_up_optimiser()
analysis.run()

print('Print max Result')
plot_form(form, show_q=False, cracks=True).show()

file_address = compas_tno.get('test.json')
form.to_json(file_address)

view_thrusts([form, form_min]).show()