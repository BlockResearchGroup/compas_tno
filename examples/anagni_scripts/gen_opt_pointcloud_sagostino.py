import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import view_solution
from compas_tno.viewers import view_shapes_pointcloud
from compas_tno.viewers import view_shapes
from compas_tno.viewers import view_meshes
from compas_tno.viewers import view_mesh
from compas.datastructures import Mesh
from compas_plotters import MeshPlotter
import json
import os


# ----------------------------------------------------------------------
# ----------- EXAMPLE OF MIN THRUST FOR DOME WITH RADIAL  FD -----------
# ----------------------------------------------------------------------

# Basic parameters

thk = 0.55
error = 0.0
x0 = 2.346
xf = 5.261
y0 = 14.97
yf = 18.450
k = 1.0
n = 2
ro_fill = 14.0
masonry_ro = 16.0
help_ = 0.50

objective = 'min'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['ind', 'zb']
features = ['fixed']  # , 'symmetry']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'current'
gradients = True
max_iter = 5000

# ------------------- FormDiagram --------------------

folder = '/Users/mricardo/compas_dev/me/anagni/'

file_formdiagram = os.path.join(folder, 'sagostino_vault_formdiagram.json')

mesh = Mesh.from_json(file_formdiagram)
vertices, faces = mesh.to_vertices_and_faces()
form = FormDiagram.from_vertices_and_faces(vertices, faces)

# Find the supports on the boundary
boundary = [0, 1, 2, 3, 4, 9, 20, 31, 57, 89, 111, 137,
            148, 159, 167, 166, 163, 164, 165, 158, 147, 136, 110, 78, 56, 30, 19, 8]

inner_ring = list(form.vertices_on_boundary())
print(inner_ring)
for key in boundary:  # form.vertices_on_boundary():
    form.vertex_attribute(key, 'is_fixed', True)

for u, v in form.edges():
    if form.vertex_attribute(u, 'is_fixed') and form.vertex_attribute(v, 'is_fixed'):
        form.edge_attribute((u, v), '_is_edge', False)

# plotter = MeshPlotter(form, figsize=(5, 5))
# plotter.draw_edges()
# # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})])
# plotter.draw_vertices(text={key: key for key in form.vertices()})
# plotter.draw_faces()
# plotter.show()

# ----------------------- Point Cloud -----------------------

file_name = 'sagostino_dome'  # This point cloud considers the extrados with 25 cm only
# Note: The fill loads must be added
pointcloud = folder + file_name + '.json'

points_ub = []
points_lb = []
tol = 10e-4

with open(pointcloud) as json_file:
    data = json.load(json_file)
    for key, pt in data['UB'].items():
        points_ub.append(pt)
    for key, pt in data['LB'].items():
        points_lb.append(pt)

print(points_ub)
print(points_lb)

# triangulated_shape = Shape.from_pointcloud(points_lb, points_ub)
# view_shapes_pointcloud(triangulated_shape).show()

form_base = form

# structured_shape = Shape.from_pointcloud_and_topology(form_base, points_lb, points_ub)
# view_shapes_pointcloud(structured_shape).show()

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
vault = Shape.from_pointcloud_and_topology(form, points_lb, points_ub, data={'type': 'general', 't': 0.0, 'thk': thk})
# vault.store_normals(plot=False)

vault.ro = masonry_ro

area = vault.middle.area()
swt = vault.compute_selfweight()

print('Interpolated Volume Data:')
print('Self-weight is: {0:.2f}'.format(swt))
print('Area is: {0:.2f}'.format(area))

# view_shapes(vault).show()

from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_from_shape

apply_selfweight_from_shape(form, vault)
apply_envelope_from_shape(form, vault)

# plotter = MeshPlotter(form, figsize=(5, 5))
# plotter.draw_edges()
# # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})])
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz')) for key in form.vertices()})
# plotter.draw_faces()
# plotter.show()

lantern_weight = 522.10
lantern_i = lantern_weight/len(inner_ring)
for key in inner_ring[:-1]:
    pz = form.vertex_attribute(key, 'pz')
    pz_lantern = pz - lantern_i
    print(pz, pz_lantern)
    form.vertex_attribute(key, 'pz', pz_lantern)

for key in boundary:
    zlb = form.vertex_attribute(key, 'lb')
    form.vertex_attribute(key, 'lb', zlb - help_)
    zub = form.vertex_attribute(key, 'ub')
    form.vertex_attribute(key, 'ub', zub + help_)

# plotter = MeshPlotter(form, figsize=(5, 5))
# plotter.draw_edges()
# # plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})])
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz')) for key in form.vertices()})
# plotter.draw_faces()
# plotter.show()

# view_meshes([vault.intrados, vault.extrados, vault.fill]).show()

# --------------------- 3. Create Starting point with TNA ---------------------

from compas_tno.problems import initialize_loadpath
from compas_tno.algorithms import form_update_with_parallelisation

form_lp = os.path.join(folder, file_name + '-lp.json')

# try:
#     print('Found and loaded LP file:', form_lp)
#     form = FormDiagram.from_json(form_lp)
# except:
# form_update_with_parallelisation(form, plot=True)
initialize_loadpath(form)
# form = form.form_update_with_parallelisation(plot=False)
# form.initialise_loadpath()
# print('back_here')
form.to_json(form_lp)
print('Computed and saved LP file:', form_lp)
# plot_form(form, show_q=True).show()
# view_mesh(form).show()

# --------------------- 4. Create Minimisation Optimiser ---------------------

optimiser = Optimiser()
optimiser.settings['library'] = solver
optimiser.settings['solver'] = solver
optimiser.settings['constraints'] = constraints
optimiser.settings['variables'] = variables
optimiser.settings['objective'] = objective
optimiser.settings['features'] = features
optimiser.settings['printout'] = True
optimiser.settings['plot'] = False
optimiser.settings['find_inds'] = True
optimiser.settings['starting_point'] = starting_point
optimiser.settings['qmax'] = 10000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['max_iter'] = max_iter

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

if objective == 't':
    n_reduction = - 1 * analysis.optimiser.fopt
    thk_min = thk - 2 * n_reduction
    print('Approx. Minimum THK:', thk_min)
    # data_shape['thk'] = thk_min

plot_form(form, show_q=False, cracks=True).show()

# analytical_shape = Shape.from_library(data_shape)
# form.envelope_from_shape(analytical_shape)

# plot_form(form, show_q=False, cracks=True).show()

file_solution = os.path.join(folder, file_name + '_lantern_form_' + objective + '.json')
form.to_json(file_solution)

file_shape = os.path.join(folder, file_name + '_lantern_shape_' + objective + '.json')
vault.to_json(file_shape)

print('Saved the solutions in:')
print(file_solution)
print(file_shape)

view_solution(form, vault).show()
