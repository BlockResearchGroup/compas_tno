from compas_tno.diagrams import FormDiagram
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer

from compas_tno.problems import initialise_form
from compas_tno.utilities import rollers_on_openings


type_structure = 'crossvault'
type_formdiagram = 'fan_fd'
thk = 0.50
xyspan = [[0.0, 10.0], [0.0, 10.0]]
discretisation = 4
c = 0.1

data_formdiagram = {
    'type': type_formdiagram,
    'xy_span': xyspan,
    'x0': 0,
    'discretisation': discretisation,
    'fix': 'corners'
}

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': 20,
    'xy_span': xyspan,
    't': 0.0,
}

fopts = []
cs_ok = []

form = FormDiagram.from_library(data_formdiagram)
form_base = FormDiagram.from_library(data_formdiagram)

rollers_on_openings(form, xy_span=xyspan, max_f=5.0, constraint_directions='all')
# M = initialise_form(form)

obj = 'min'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']
features = ['fixed']
starting_point = 'tna'
make_video = False
autodiff = False

# Create shape/diagram

vault = Shape.from_library(data_shape)
# vault.ro = 0.2

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
optimiser.settings['save_iterations'] = make_video
optimiser.settings['autodiff'] = autodiff

# Create analysis

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
# analysis.apply_envelope_on_xy(c=c)

analysis.set_up_optimiser()

# # to view starting point
# view = Viewer(form, vault)
# view.draw_thrust()
# view.draw_shape()
# view.draw_force()
# view.show()

analysis.run()

# to view solution
view = Viewer(form, vault)
view.settings['force.scale'] = 0.05*10
view.draw_thrust()
view.draw_shape()
view.draw_cracks()
view.draw_force()
view.show()

if analysis.optimiser.exitflag == 0:
    fopts.append(analysis.optimiser.fopt)
    cs_ok.append(c)

plotter = TNOPlotter(form)
plotter.settings['scale.reactions'] = 0.05/10
plotter.draw_base_form(form_base=form_base)
plotter.draw_form()
plotter.draw_supports()
plotter.draw_cracks()
plotter.draw_constraints()
plotter.draw_reactions()
plotter.show()

print(fopts)
print(cs_ok)
