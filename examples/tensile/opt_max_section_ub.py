from numpy.lib.function_base import blackman
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.viewers import Viewer
from compas_plotters import MeshPlotter

from compas.utilities import i_to_red

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_on_q

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

from compas.geometry import Polygon

import compas_tno
import os

span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'arch'  # write the type of form diagram you want and is in the file shape
type_structure = 'arch'
thk = 0.1
discretisation_shape = 10 * discretisation

H = 2.0
L = 5

save = False
solutions = {}

objective = 'max_section'  # try 'max'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 'tub']  # in the futture add 'tlb' as variables
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'loadpath'
tubmax = 0.5

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'H': H,
    'L': L,
    'x0': 0,
    'total_nodes': discretisation,
}

form = FormDiagram.from_library(data_diagram)
# plot_simple_form(form).show()

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'H': H,
    'L': L,
    'x0': 0,
    'discretisation': discretisation_shape,
    'b': 0.3,
    't': 0.0,
}

vault = Shape.from_library(data_shape)
vault.ro = 20.0

# ------------------------------------------------------------
# -----------------------  INITIALISE   ----------------------
# ------------------------------------------------------------

# Apply Selfweight and Envelope

apply_envelope_from_shape(form, vault)
apply_selfweight_from_shape(form, vault)
apply_bounds_on_q(form, qmax=1e-6)

# add point load in vertex 5
pz3 = form.vertex_attribute(5, 'pz')
form.vertex_attribute(5, 'pz', pz3 - 0.5)

# view = Viewer(form, vault)
# view.view_shape()
# view.show()

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
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = False
optimiser.settings['max_iter'] = 500
optimiser.settings['gradient'] = True
optimiser.settings['jacobian'] = True
optimiser.settings['printout'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['tubmax'] = tubmax
optimiser.settings['sym_loads'] = False

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.set_up_optimiser()
analysis.run()

# ----------------------- 6. additional thickness ----------------------

path = compas_tno.get('')
address = os.path.join(path, 'form.json')
form.to_json(address)
print('Form Saved to:', address)

from compas_tno.plotters import plot_form_xz
plot_form_xz(form, vault).show()

# plotter = MeshPlotter(form)
# plotter.draw_edges()
# plotter.draw_vertices(
#     text={key: round(form.vertex_attribute(key, 'tub'), 2) for key in form.vertices() if form.vertex_attribute(key, 'tub') > 0.001},
#     facecolor={key: i_to_red(form.vertex_attribute(key, 'tub')/tubmax) for key in form.vertices()}
#     )
# plotter.show()

plot_form(form).show()
view = Viewer(form, vault)
view.show_solution()
