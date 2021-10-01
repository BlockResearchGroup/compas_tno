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

thk = 0.25
error = 0.0
x0 = -4.431
xf = 1.313
y0 = 14.968
yf = 18.312

# small vault
# x0 = -4.190
# y0 = 15.135
# xf = 1.073
# yf = 18.145

k = 1.0
n = 2
ro_fill = 14.0

objective = 'max'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['ind', 'zb']
features = ['fixed']  # , 'symmetry']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'current'
gradients = True
max_iter = 500
help_ = 0.035

# ------------------- FormDiagram --------------------

folder = '/Users/mricardo/compas_dev/me/anagni/'

file_formdiagram = os.path.join(folder, 'top_vault_less_discr.json')

mesh = Mesh.from_json(file_formdiagram)
vertices, faces = mesh.to_vertices_and_faces()
form = FormDiagram.from_vertices_and_faces(vertices, faces)

# Find the supports on the corners
for key in form.vertices():
    deg = form.vertex_degree(key)
    x, y, _ = form.vertex_coordinates(key)
    a = abs(x - x0) < 0.1
    b = abs(y - y0) < 0.1
    c = abs(x - xf) < 0.1
    d = abs(y - yf) < 0.1
    test = [a, b, c, d]
    if sum(test) == 2:
        form.vertex_attribute(key, 'is_fixed', True)

plotter = MeshPlotter(form, figsize=(5, 5))
plotter.draw_edges()
plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})])
plotter.draw_faces()
plotter.show()

# ----------------------- Point Cloud -----------------------

file_name = 'top_vault'  # This point cloud considers the extrados with 25 cm only
# Note: The fill loads must be added
pointcloud = folder + file_name + '.json'

points_ub = []
points_lb = []
points_fill = []
tol = 10e-4

with open(pointcloud) as json_file:
    data = json.load(json_file)
    for key, pt in data['UB'].items():
        points_ub.append(pt)
    for key, pt in data['LB'].items():
        points_lb.append(pt)
    for key, pt in data['FILL'].items():
        points_fill.append(pt)

print(points_ub)
print(points_lb)
print(points_fill)

# triangulated_shape = Shape.from_pointcloud(points_lb, points_ub)
# view_shapes_pointcloud(triangulated_shape).show()

form_base = form

# structured_shape = Shape.from_pointcloud_and_topology(form_base, points_lb, points_ub)
# view_shapes_pointcloud(structured_shape).show()

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
vault = Shape.from_pointcloud_and_topology(form, points_lb, points_ub, fill_pts=points_fill, data={'type': 'general', 't': 0.0, 'thk': thk})
# more improved, considers the real middle
# # vault = Shape.from_meshes_and_formdiagram(form, structured_shape.intrados, structured_shape.extrados, data={'type': 'general', 't': 0.0, 'thk': thk})
vault.store_normals(plot=False)
# view_shapes_pointcloud(vault).show()
# view_normals(vault).show()

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

pzt0 = 0
for key in form.vertices():
    pzi = form.vertex_attribute(key, 'pz')
    pzt0 += pzi

print('Selfweight not considering the fill is:', pzt0)

# Add the fill load to the nodes affected

# kinks = [13, 21, 121, 129]
# kinks2 = [17, 28, 114, 125]

# plotter = MeshPlotter(form, figsize=(5, 5))
# plotter.draw_edges()
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 2) for key in form.vertices()})
# plotter.draw_faces()
# plotter.show()

for key in vault.fill.vertices():
    pz = form.vertex_attribute(key, 'pz')
    zub = form.vertex_attribute(key, 'ub')
    z_fill = vault.fill.vertex_attribute(key, 'z')
    z_diff = (z_fill - zub)
    if z_diff > 0.01:
        print('add load to:', key)
        ai = vault.middle.vertex_area(key)
        pz_fill = -1 * ai * z_diff * ro_fill  # downwards loads are negative
        pzt = pz + pz_fill
        form.vertex_attribute(key, 'pz', pzt)
        # print('improve height of:', key)
        # zubnew = zub + help_
        # form.vertex_attribute(key, 'ub', zubnew)

# plotter = MeshPlotter(form, figsize=(5, 5))
# plotter.draw_edges()
# plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 2) for key in form.vertices()})
# plotter.draw_faces()
# plotter.show()

pztfinal = 0
for key in form.vertices():
    pzi = form.vertex_attribute(key, 'pz')
    pztfinal += pzi

print('Selfweight with the considering the fill is:', pztfinal)
print('Fill weight is:', pztfinal - pzt0)

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
plot_form(form, show_q=True).show()
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
optimiser.settings['qmax'] = 1000.0
optimiser.settings['gradient'] = gradients
optimiser.settings['jacobian'] = gradients
optimiser.settings['max_iter'] = max_iter

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, form, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

if objective == 't':
    n_reduction = - 1 * analysis.optimiser.fopt
    thk_min = thk - 2 * n_reduction
    print('Approx. Minimum THK:', thk_min)
    # data_shape['thk'] = thk_min

plot_form(form, show_q=False, cracks=True).show()

file_solution = os.path.join(folder, file_name + '_less_discr_form_' + objective + '.json')
form.to_json(file_solution)

file_shape = os.path.join(folder, file_name + '_less_discr_shape_' + objective + '.json')
vault.to_json(file_shape)

view_solution(form, vault).show()
