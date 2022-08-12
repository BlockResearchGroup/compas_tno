from compas_tno.diagrams import FormDiagram
from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis
from compas_tno.plotters import TNOPlotter
from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities.form import displacement_map_4parabolas

xy_span_pattern = [[0, 10.0], [0, 10.0]]
xy_span_shape = [[-1, 11.0], [-1, 11.0]]

form = FormDiagram.create_cross_form(xy_span=xy_span_pattern)
form_base = FormDiagram.create_cross_form(xy_span=xy_span_pattern)
shape = Shape.create_crossvault(xy_span=xy_span_shape)

objective = 'min'
solver = 'SLSQP'
variables = ['q', 'zb', 'delta']
constraints = ['funicular', 'envelope', 'envelopexy', 'displ_map']
features = ['']
starting_point = 'loadpath'
max_iter = 500
c = 0.25

# new info for rdisplacement map:
dX = displacement_map_4parabolas(form, tol=0.1)

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['objective'] = objective
optimiser.settings['features'] = features
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
# optimiser.settings['find_inds'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['qmax'] = 10000.0
optimiser.settings['max_iter'] = max_iter

optimiser.settings['displ_map'] = dX
optimiser.settings['max_delta'] = 0.5
optimiser.settings['delta0'] = 0.1

analysis = Analysis.from_elements(shape, form, optimiser)
analysis.apply_envelope()
analysis.apply_selfweight(normalize_loads=False)
analysis.apply_envelope_on_xy(c=c)
analysis.set_up_optimiser()
analysis.run()

print(optimiser.xopt)
deltaopt = optimiser.xopt[-1]
print('Delta Optimum:', deltaopt)

view = Viewer(form)
view.show_solution()

plotter = TNOPlotter(form)
plotter.draw_base_form(form_base=form_base)
plotter.draw_form()
plotter.draw_supports()
plotter.draw_cracks()
plotter.draw_constraints()
plotter.show()
