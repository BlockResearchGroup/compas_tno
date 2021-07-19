from compas_tno.diagrams import FormDiagram
from compas.datastructures import Mesh
from compas_tno.shapes import Shape
from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_superimposed_diagrams
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_thrust

from compas_tno.optimisers.optimiser import Optimiser
from compas_tno.analysis.analysis import Analysis
import os
from compas_tno.plotters import save_csv
from compas_tno.plotters import diagram_of_thrust

from compas_tno.algorithms import apply_sag

span = 10.0
k = 1.0
discretisation = 10
type_formdiagram = 'fan_fd'
type_structure = 'crossvault'
thk = 0.50
discretisation_shape = 4 * discretisation

type_ = 'D5'
pattern_load = '/Users/mricardo/compas_dev/me/loadpath/corner/topology/' + type_ + '_complete.json'

c = 0.1

thk = 0.50

# Create form diagram

data_diagram = {
    'type': type_formdiagram,
    'xy_span': [[0, span], [0, k*span]],
    'discretisation': discretisation,
    'fix': 'corners'
}

form = FormDiagram.from_library(data_diagram)
from compas.geometry import distance_point_point_xy

if pattern_load:
    mesh = Mesh.from_json(pattern_load)
    form = FormDiagram.from_mesh(mesh)

    v, f = form.to_vertices_and_faces()
    form = FormDiagram.from_mesh(Mesh.from_vertices_and_faces(v, f))

    corner_pts = [[0.0, 0.0], [0.0, 10.0], [10.0, 10.0], [10.0, 0.0]]
    tol = 0.0001
    corner_keys = []

    for i in range(len(corner_pts)):
        for key in form.vertices():
            pt = form.vertex_coordinates(key)
            if distance_point_point_xy(corner_pts[i], pt) < tol:
                corner_keys.append(key)
                break

    for key in corner_keys:
        form.vertex_attribute(key, 'is_fixed', True)

from compas_plotters import MeshPlotter
plotter = MeshPlotter(form, fig_size=(6, 6))
plotter.draw_edges(width=2.0)
# plotter.draw_vertices(text={key: key for key in form.vertices()})
plotter.draw_vertices(radius=0.03)
plotter.show()

pattern_out = '/Users/mricardo/compas_dev/me/loadpath/corner/topology/' + type_ + '.json'
print('Saved to: ', pattern_out)
form.to_json(pattern_out)
plot_form(form, show_q=False).show()

# Create shape

data_shape = {
    'type': type_structure,
    'thk': thk,
    'discretisation': discretisation_shape,
    'xy_span': [[0, span], [0, k*span]],
    # 'hc': hc,
    # 'hm': None,
    # 'he': None,
    'center': [5.0, 5.0],
    'radius': span/2,
    't': 0.0,
}

vault = Shape.from_library(data_shape)

# ------------------------------------------------------------
# -----------------------  INITIALISE   ----------------------
# ------------------------------------------------------------

# Apply Selfweight and Envelope

from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_on_xy
from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_bounds_on_q

apply_envelope_from_shape(form, vault)
apply_selfweight_from_shape(form, vault)
if 'lambd' in variables:
    apply_horizontal_multiplier(form, lambd=lambd)

if 'envelopexy' in constraints:
    apply_envelope_on_xy(form, c=c)
apply_bounds_on_q(form, qmax=0.0)

form_base = form.copy()

# ------------------------------------------------------------
# ------------------- Proper Implementation ------------------
# ------------------------------------------------------------

optimiser = Optimiser()
# optimiser.data['library'] = 'SLSQP'
# optimiser.data['solver'] = 'SLSQP'
optimiser.data['library'] = 'IPOPT'
optimiser.data['solver'] = 'IPOPT'
optimiser.data['constraints'] = ['funicular', 'envelope']
# optimiser.data['variables'] = ['q', 'sym']
optimiser.data['variables'] = ['q', 'zb', 'sym']
optimiser.data['objective'] = 'min'
optimiser.data['plot'] = False
optimiser.data['find_inds'] = False
optimiser.data['printout'] = True
optimiser.data['qmax'] = 1000.0
optimiser.data['gradient'] = True
optimiser.data['jacobian'] = True
optimiser.data['derivative_test'] = False
optimiser.data['max_iter'] = 500

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

weight = 0
for key in form.vertices():
    weight += form.vertex_attribute(key, 'pz')

thrust = abs(optimiser.fopt)
print('Ratio Thrust/Weight:', thrust/weight)

folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', type_structure, type_formdiagram, 'with_Xb_mov_c_' + str(c))
title = type_structure + '_' + type_formdiagram + '_discr_' + str(discretisation)

if type_:
    folder = os.path.join('/Users/mricardo/compas_dev/me', 'general_opt', type_structure, type_, 'mov_c_' + str(c))
    title = type_structure + '_' + type_ + '_discr_' + str(discretisation)

os.makedirs(folder, exist_ok=True)

save_form = os.path.join(folder, title)
address = save_form + '_' + optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'

if optimiser.exitflag == 0 or optimiser.data['library'] == 'IPOPT':
    form.to_json(address)
    print('Form Saved to:', address)

# Viewing
plot_superimposed_diagrams(form, form_base).show()
plot_form(form, show_q=False, cracks=True).show()
view_solution(form).show()
