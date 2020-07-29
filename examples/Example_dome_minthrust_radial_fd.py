import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
from copy import deepcopy

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8, 20]
n = 2
plot = False

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [n*discretisation[0], n*discretisation[1]],
    'center': [5.0, 5.0],
    'radius': radius,
    't': 1.0
}

dome = Shape.from_library(data_shape)
swt = dome.compute_selfweight()
print('Selfweight computed:', swt)
print('Dome created!')

# from compas_tno.viewers import view_shapes
# view_shapes(dome).show()

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

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.selfweight_from_shape(dome)
form = form.initialise_tna(plot=False)
# form = form.initialise_loadpath()
if plot:
    plot_form(form, show_q=False).show()

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
optimiser.data['qmax'] = 3000.0
optimiser.data['gradient'] = True
optimiser.data['jacobian'] = True
print(optimiser.data)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

if plot:
    plot_form(form, show_q=False, cracks=True).show()
    print('Print min Result')
form_min = deepcopy(form)

file_address = compas_tno.get('test-min.json')
form.to_json(file_address)

form, force = form.reciprocal_from_form(plot=plot)

# --------------------- 6. Run it for MAX ---------------------

optimiser.data['objective'] = 'max'
analysis.set_up_optimiser()
analysis.run()

if plot:
    print('Print max Result')
    plot_form(form, show_q=False, cracks=True).show()

file_address = compas_tno.get('test-max.json')
form.to_json(file_address)

form, force = form.reciprocal_from_form(plot=plot)

view_thrusts([form, form_min]).show()
