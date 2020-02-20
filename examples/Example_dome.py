
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes.shape import Shape
from compas_tno.viewers.shapes import view_shapes
from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.plotters.plotters import plot_form
from compas_tno.analysis.analysis import Analysis

# ----------------------------------------------------------------------
# -----------EXAMPLE OF MIN and MAX THRUST FOR DOME --------------------
# ----------------------------------------------------------------------

# Basic parameters

thk = 1.0
radius = 5.0
type_structure = 'dome'
type_formdiagram = 'radial_fd'

# ----------------------- Create Dome shape ---------------------------

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

# ----------------------- Create Form Diagram ---------------------------

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
# form = FormDiagram.from_json('/Users/mricardo/compas_dev/me/minmax/dome/radial/radial_discr_8_16_lp.json')
form.overview_forces()
print('Form Diagram Created!')
print(form)
plot_form(form, show_q=False, fix_width=True).show()

# ----------------------- Create Optimiser ----------------------------

optimiser = Optimiser()
optimiser.data['library'] = 'Scipy'
optimiser.data['solver'] = 'slsqp'
optimiser.data['constraints'] = ['funicular', 'envelope', 'reac_bounds']
optimiser.data['variables'] = ['ind', 'zb']
optimiser.data['objective'] = 'min'
optimiser.data['qmax'] = 100.0
print(optimiser.data)

# ----------------------- Create Analysis Model ------------------------

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser() # Find independent edges
# plot_form(form, show_q=False).show()
analysis.run()

plot_form(form, show_q=False).show()
