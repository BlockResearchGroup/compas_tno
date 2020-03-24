from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis.analysis import Analysis
from compas_tno.viewers.thrust import view_thrust

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR DOME --------------------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.5
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'
discretisation = [8,16]

# ----------------------- 1. Create Dome shape ---------------------------

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [8, 16],
    'center': [5.0, 5.0],
    'radius': radius,
    't' : 10.0
}

dome = Shape.from_library(data_shape)
swt = dome.compute_selfweight()
print('Dome created!')
# view_shapes(dome).show()

# ----------------------- 2. Create Form Diagram ---------------------------

data_diagram = {
    'type': type_formdiagram,
    'center': [5.0, 5.0],
    'radius': radius,
    'discretisation': [8, 16],
    'r_oculus': 0.0,
    'diagonal': False,
    'partial_diagonal': False,
}

form = FormDiagram.from_library(data_diagram)
print('Form Diagram Created!')
print(form)
# plot_form(form, show_q=False, fix_width=False).show()

# --------------------- 3.1. Create Convex Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'MATLAB'
optimiser.data['solver'] = 'SDPT3'
optimiser.data['constraints'] = ['funicular']
optimiser.data['variables'] = ['ind']
optimiser.data['objective'] = 'loadpath'
optimiser.data['printout'] = True
optimiser.data['plot'] = False
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 150.0
print(optimiser.data)

# # # -------------- 3.2. Create Analysis Model and Run Convex Opt --------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.set_up_optimiser() # Find independent edges
analysis.run()
# plot_form(form, show_q=False).show()
# view_thrust(form).show()

file_adress = '/Users/mricardo/compas_dev/me/reformulation/test.json'
form.to_json(file_adress)
# form.from_json(file_adress)

# --------------------- 4.1 Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.data['library'] = 'MMA'
optimiser.data['solver'] = 'MMA'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['printout'] = True
optimiser.data['solver_options']['derivatives'] = 'DF_reduced'  # 'DF_brute' 'DF_reduced' and 'analytical' in process.
optimiser.data['plot'] = True
optimiser.data['find_inds'] = True
optimiser.data['qmax'] = 50.0
print(optimiser.data)

# --------------------- 4.2 Create Minimisation Optimiser ---------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser() # Find independent edges
analysis.run()

form = analysis.form
plot_form(form, show_q=False).show()

file_adress = '/Users/mricardo/compas_dev/me/reformulation/test.json'
# form.to_json(file_adress)

view_thrust(form).show()