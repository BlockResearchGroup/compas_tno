from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.plotters import TNOPlotter
from compas_plotters import Plotter
from compas_tno.viewers import Viewer

from compas.utilities import i_to_red, i_to_blue

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_on_q
from compas_tno.utilities import apply_bounds_tub_tlb

from compas_tno.optimisers import Optimiser
from compas_tno.analysis import Analysis

import compas_tno
import os

thk = 0.50
discretisation_shape = 10
discretisation = 10
type_structure = 'crossvault'
type_formdiagram = 'cross_fd'  # write the type of form diagram you want and is in the file shape

save = False
solutions = {}

objective = 'max_section'  # try 'max'
solver = 'IPOPT'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb', 'tub']  # in the future add 'tlb' as variables
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'tna'

# Set the maximum allowable increase of each thickness
tubmax = 0.5
tlbmax = 0.5
disct_reac_max = 0.5

# Create form diagram

data_formdiagram = {
    'type': type_formdiagram,
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    'x0': 0,
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_formdiagram)
# plot_simple_form(form).show()

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': [discretisation_shape, discretisation_shape],
    'xy_span': [[0.0, 10.0], [0.0, 10.0]],
    't': 0.0,
}

crossvault = Shape.from_library(data_shape)
crossvault.ro = 20.0
vault_swt = crossvault.compute_selfweight()
# print (vault_swt)

plotter = Plotter()
plotter.fontsize = 12
artist = plotter.add(form)
artist.draw_vertexlabels()
plotter.zoom_extents()
plotter.show()

# ------------------------------------------------------------
# -----------------------  INITIALISE   ----------------------
# ------------------------------------------------------------

# Apply Selfweight and Envelope

apply_envelope_from_shape(form, crossvault)
apply_selfweight_from_shape(form, crossvault)
apply_bounds_on_q(form, qmax=1e-6)
apply_bounds_tub_tlb(form, tubmax=tubmax, tlbmax=tlbmax)

# for key in form.vertices_where({'is_fixed': True}):
#     x = form.vertex_coordinates(key)[0]
#     if abs(x - x0) < 10e-3:
#         form.vertex_attribute(key, 'tub_reacmax', [- disct_reac_max, 0])  # left reaction @ (0, 0)
#     if abs(x - L) < 10e-3:
#         form.vertex_attribute(key, 'tub_reacmax', [+ disct_reac_max, 0])  # right reaction @ (L, 0)

# add load

pc = 20
lambda_mult = pc*vault_swt/100
pzv_add = -1

# for key in range (6):
#     pzkey = form.vertex_attribute(key+54, 'pz')
#     form.vertex_attribute(key+54, 'pz', pzkey + lambda_mult*pzv_add)

pzkey = form.vertex_attribute(59, 'pz')
form.vertex_attribute(59, 'pz', pzkey + lambda_mult*pzv_add)

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
optimiser.settings['sym_loads'] = False

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(crossvault, form, optimiser)
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

# ----------------------- 6. additional thickness ----------------------

# path = compas_tno.get('')
# address = os.path.join(path, 'form_tub_tlb.json')
# form.to_json(address)
# print('Form Saved to:', address)

# path = compas_tno.get('')
# address = os.path.join(path, 'arch20.json')
# crossvault.to_json(address)
# print('arch Saved to:', address)

# folder = compas_tno.get('')
# folder = os.path.join(folder, 'test')
# os.makedirs(folder, exist_ok=True)
# title = 'crossvault' + '_' + 'cross_fd' + '_discr_' + str(discretisation) + '_' + optimiser.settings['objective'] + 'bottom_line' + '_thk_' + str(100*thk) + '_pct_stw_' + str(pc) + '.json'
# save_form = os.path.join(folder, title)
# form.to_json(save_form)

# print('Solution Saved at:', save_form)

plotter = Plotter()
plotter.fontsize = 12
artist = plotter.add(form)
print(artist)
artist.draw_vertices(
    text={key: round(form.vertex_attribute(key, 'tub'), 2) for key in form.vertices() if form.vertex_attribute(key, 'tub') > 0.001},
    facecolor={key: i_to_red(form.vertex_attribute(key, 'tub')/tubmax) for key in form.vertices()}
    )
plotter.show()

view = Viewer(form, crossvault)
view.settings['scale.loads'] = 0.5
view.draw_loads()
view.draw_shape()
view.draw_cracks()
view.draw_thrust()
view.show()

# print (pz_tot)
# print(pz3)
