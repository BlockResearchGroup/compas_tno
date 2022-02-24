from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer


# Geometry parameters

radius = 5.0
thk = 0.5
discretisation = [8, 10]
discretisation_shape = [2*discretisation[0], 2*discretisation[1]]

# Parameters Optimisations

obj = 'min'
solver = 'SLSQP'
constraints = ['funicular', 'envelope', 'reac_bounds']
variables = ['q', 'zb']
features = ['fixed']
starting_point = 'tna'

# Create shape/diagram

dome = Shape.create_dome(thk=thk, radius=radius, discretisation=discretisation_shape, t=0.5)

form = FormDiagram.create_circular_radial_form(discretisation=discretisation, radius=radius)

# Create optimiser

optimiser = Optimiser()
optimiser.settings['objective'] = obj
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['starting_point'] = starting_point
optimiser.settings['derivative_test'] = False
optimiser.settings['printout'] = True
optimiser.settings['plot'] = True

# Create analysis

analysis = Analysis.from_elements(dome, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()

analysis.run()

view = Viewer(form, dome)
view.show_solution()

# x0 = optimiser.x0

# print(x0)
