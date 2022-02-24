from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer

from compas_tno.plotters import plot_form_xz

type_formdiagram = 'arch'  # write the type of form diagram you want and is in the file shape
type_structure = 'arch'
discretisation = 20

H = 0.5
L = 1.0
b = 0.30
x0 = 0.0
thk = 0.10

objective = 'max'  # try 'max'
solver = 'SLSQP'  # try SLSQP
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb']  # in the future add 'tlb' as variables
features = ['fixed']
starting_point = 'loadpath'

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'H': H,
    'L': L,
    'x0': x0,
    'total_nodes': discretisation,
}

form = FormDiagram.from_library(data_diagram)

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'H': H,
    'L': L,
    'x0': x0,
    'discretisation': discretisation*10,
    'b': b,
    't': 0.1,
}

arch = Shape.from_library(data_shape)

# Define optimisation settings

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['objective'] = objective
optimiser.settings['plot'] = False
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
optimiser.settings['printout'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['sym_loads'] = False

# Set up and run analysis

analysis = Analysis.from_elements(arch, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# Visualise solution

plot_form_xz(form, arch, radius=0.02, fix_width=True, max_width=5.0).show()

view = Viewer(form, arch)
view.settings['scale.loads'] = 1.0
view.view_loads()
view.view_shape()
view.view_cracks()
view.draw_thrust()
view.view_reactions()
view.show()
