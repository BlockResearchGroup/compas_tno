
import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrusts
import os

# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

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
# plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3. Create Starting point with TNA ---------------------

# form = form.initialise_tna(plot=False)
form.selfweight_from_shape(dome)
orm = form.initialise_tna(plot=False)
# form = form.initialise_loadpath()
plot_form(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', ' symmetry']
optimiser.data['variables'] = ['ind']
optimiser.data['objective'] = 'target'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 3000.0
print(optimiser.data)

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

print('Print target Result - only ind variables')
plot_form(form, show_q=False, cracks=True).show()
form1 = form.copy()
file_address = os.path.join(compas_tno.get(''),'test1.json')
form.to_json(file_address)

# form, force = form.reciprocal_from_form(plot=True)

# --------------------- 6. Run it for MAX ---------------------

optimiser.data['variables'] = ['ind', 'zb']
analysis.set_up_optimiser()
analysis.run()

print('Print target Result - only ind and zb variables')
plot_form(form, show_q=False, cracks=True).show()

file_address = os.path.join(compas_tno.get(''),'test2.json')
form.to_json(file_address)

# form, force = form.reciprocal_from_form(plot=True)

view_thrusts([form, form1]).show()
