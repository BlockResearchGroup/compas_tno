from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import TNOPlotter
from compas_plotters import Plotter
from compas_tno.viewers import Viewer

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

import compas_tno
import os

thk = 0.50
discretisation_shape = 28
discretisation = 14
type_structure = 'crossvault'
type_formdiagram = 'fan_fd'  # write the type of form diagram you want and is in the file shape
delta = 0.5
xyspan_shape = [[0.0 - delta, 10.0+delta], [0.0 - delta, 10.0 + delta]]
xyspan = [[0.0, 10.0], [0.0, 10.0]]

save = False
solutions = {}

objective = 'min'  # try 'max'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']  # in the future add 'tlb' as variables
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'loadpath'

# Create form diagram

data_formdiagram = {
    'type': type_formdiagram,
    'xy_span': xyspan,
    'x0': 0,
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_formdiagram)
# plot_simple_form(form).show()

# from compas_tno.algorithms import apply_sag
# apply_sag(form)

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation_shape, discretisation_shape],
    'xy_span': xyspan_shape,
    't': 0.0,
}

crossvault = Shape.from_library(data_shape)

# ------------------------------------------------------------
# ------------------- Proper Implementation ------------------
# ------------------------------------------------------------

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['features'] = features
optimiser.settings['objective'] = objective
optimiser.settings['plot'] = True
optimiser.settings['find_inds'] = False
optimiser.settings['max_iter'] = 500
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
optimiser.settings['printout'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['sym_loads'] = False

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(crossvault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
analysis.set_up_optimiser()

# ad = '/Users/mricardo/compas_dev/compas_tno/data/analysis-cross14.json'
# analysis.to_json(ad)

analysis.run()

# ----------------------- 6. additional thickness ----------------------

view = Viewer(form, show_grid=False)
view.settings['scale.reactions'] = 0.005
view.settings['camera.distance'] = 30
view.settings['camera.target'] = [5, 5, 0]
# view.draw_loads()
view.draw_shape()
view.draw_cracks()
view.draw_thrust()
view.draw_reactions()
view.show()

# print (pz_tot)
# print(pz3)
