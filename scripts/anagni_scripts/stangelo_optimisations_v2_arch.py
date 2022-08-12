# import compas_tno
import math
from compas_tno.algorithms import form_update_with_parallelisation
from compas_tno.problems import initialize_loadpath
from compas_tno.problems import initialize_tna
from compas_tno.algorithms import apply_sag
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_bounds_tub_tlb
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
# from compas_tno.plotters import plot_form
from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
# from compas_tno.viewers import view_shapes_pointcloud
# from compas_tno.viewers import view_shapes
# from compas_tno.viewers import view_meshes
# from compas_tno.viewers import view_mesh
from compas.datastructures import Mesh
from compas_tno.plotters import TNOPlotter
from compas_tno.utilities import move_pattern_to_origin

from numpy import array
from compas_plotters import Plotter
import json
import os

# Basic parameters

thk = 0.25
additional_thk = 0.00
error = 0.0

# #eassuress from Rhinoceros
xspan = [-4.431, 1.313]
yspan = [14.968, 18.312]

fill_loads = True
diagram_name = 'mesh-B3'

x0 = 0.0
xf = xspan[1] - xspan[0]
y0 = 0.0
yf = yspan[1] - yspan[0]
xc = (xf + x0)/2
yc = (yf + y0)/2
discretisation = 16

corners = [[x0, y0], [xf, y0], [xf, yf], [x0, yf]]
print('Corners:', corners)

k = 1.0
n = 2
ro_fill = 14.0

objective = 'min'
solver = 'IPOPT'
constraints = ['funicular', 'envelope']
variables = ['ind', 'zb']
features = ['fixed']  # , 'symmetry']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'current'
gradients = True
max_iter = 1000
help_ = 0.03

# ------------------- FormDiagram --------------------

folder = '/Users/mricardo/compas_dev/me/anagni/'

if diagram_name == 'cross_fd':

    form = FormDiagram.create_cross_form(discretisation=discretisation)
    diagram_name = diagram_name + '-' + str(discretisation)

elif diagram_name == 'fan_fd':

    form = FormDiagram.create_fan_form(discretisation=discretisation)
    diagram_name = diagram_name + '-' + str(discretisation)

elif diagram_name == 'arch':

    form = FormDiagram.create_cross_form(discretisation=discretisation)

    arch = FormDiagram.create_linear_form_diagram(L=xf, x0=x0,discretisation=discretisation+1)
    diagram_name = diagram_name + '-' + str(discretisation)

else:

    file_formdiagram = os.path.join(folder, 'meshes', 'CISM', diagram_name + '.json')

    mesh = Mesh.from_json(file_formdiagram)

    vertices, faces = mesh.to_vertices_and_faces()
    form = FormDiagram.from_vertices_and_faces(vertices, faces)

    delfaces = []
    for fkey in form.faces():
        a = form.face_area(fkey)
        if a > 3.0:
            print('delete:', fkey)
            delfaces.append(fkey)

    for fkey in delfaces:
        form.delete_face(fkey)

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
plotter.show()

trans, factors = move_pattern_to_origin(form, corners)

lines = []
vector_supports = []
for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)
    if x < xc:
        dXbi = [-1, 0, 0]
    else:
        dXbi = [1, 0, 0]
    vector_supports.append(dXbi)
    lines.append([[x, y, z], [x + dXbi[0], y + dXbi[1], z + dXbi[2]]])

dXb = array(vector_supports)

plotter = TNOPlotter(form)
plotter.draw_form(scale_width=False)
# plotter.draw_lines(lines)
plotter.draw_supports()
plotter.zoom_extents()  # why zoom extents does not work?
plotter.show()

# ----------------------- Point Cloud -----------------------

file_name = 'sangelo_vault_top_final'  # This point cloud considers the extrados with 25 cm only
# Note: The fill loads must be added
pointcloud = folder + file_name + '.json'

print('Point-cloud_location:', pointcloud)

points = {'UB': [], 'LB': [], 'FILL': []}

tol = 1e-3

with open(pointcloud) as json_file:
    data = json.load(json_file)

for tp in ['UB', 'LB', 'FILL']:
    for key, pt in data[tp].items():
        x, y, z = pt
        pt = [x - xspan[0], y - yspan[0], z]

        if abs(pt[0] - x0) < tol:
            pt[0] = pt[0] - tol * 2
        elif abs(pt[0] - xf) < tol:
            pt[0] = pt[0] + tol * 2

        if abs(pt[1] - y0) < tol:
            pt[1] = pt[1] - tol * 2
        elif abs(pt[1] - yf) < tol:
            pt[1] = pt[1] + tol * 2

        points[tp].append(pt)


form_base = form

# ------- Create shape given a topology and a point cloud --------

# roots - not considering real middle
vault = Shape.from_pointcloud_and_topology(form, points['LB'], points['UB'], fill_pts=points['FILL'], data={'type': 'general', 't': 0.0, 'thk': thk + additional_thk})
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

apply_selfweight_from_shape(form, vault)
apply_envelope_from_shape(form, vault)
apply_bounds_tub_tlb(form)

pzt0 = 0
for key in form.vertices():
    pzi = form.vertex_attribute(key, 'pz')
    pzt0 += pzi

print('Selfweight NOT considering the fill is:', pzt0)

view = Viewer(form, vault)
view.draw_shape()
view.draw_middle_shape()
view.draw_mesh(vault.fill)
view.show()

# Add the fill load to the nodes affected

# kinks = [13, 21, 121, 129]
# kinks2 = [17, 28, 114, 125]

plotter = Plotter()
artist = plotter.add(form)
artist.draw_edges()
artist.draw_vertices()
artist.draw_vertexlabels()
plotter.show()

if fill_loads:
    pzfilldic = {}
    at = 0.0
    for key in vault.fill.vertices():
        pz = form.vertex_attribute(key, 'pz')
        zub = form.vertex_attribute(key, 'ub')
        z_fill = vault.fill.vertex_attribute(key, 'z')
        z_diff = (z_fill - zub)
        if z_diff > 0.01:
            # print('add load to:', key)
            ai = form.vertex_area(key)  # since the form is planar at this point this corresponds to the projected area of the node
            pz_fill = -1 * ai * z_diff * ro_fill  # downwards loads are negative
            pzt = pz + pz_fill
            pzfilldic[key] = round(pz_fill, 1)
            form.vertex_attribute(key, 'pz', pzt)
            print('improve height of:', key)
            at += ai
            print('Node | height | area | fill load:', key, z_diff, ai, pz_fill)
        # zubnew = zub + help_
        # form.vertex_attribute(key, 'ub', zubnew)

    print('total area fill:', at)
    plotter = Plotter()
    artist = plotter.add(form)
    artist.draw_edges()
    # artist.draw_vertexlabels(text={key: round(form.vertex_attribute(key, 'pz'), 2) for key in form.vertices()})
    artist.draw_vertexlabels(text=pzfilldic)
    plotter.zoom_extents()
    plotter.show()

pztfinal = 0
for key in form.vertices():
    pzi = form.vertex_attribute(key, 'pz')
    pztfinal += pzi

print('Selfweight considering the fill is:', pztfinal)
print('Fill weight is:', pztfinal - pzt0)

# Copy information to the arch based on the closest node

from compas.geometry import distance_point_point_xy
for key in arch.vertices():
    pt_arch = arch.vertex_coordinates(key)
    for pkey in form.vertices():
        pt_pattern = form.vertex_coordinates(pkey)
        if distance_point_point_xy(pt_arch, pt_pattern) < 0.01:
            zub = form.vertex_attribute(pkey, 'ub')
            zlb = form.vertex_attribute(pkey, 'lb')
            target = form.vertex_attribute(pkey, 'target')
            pz = form.vertex_attribute(pkey, 'pz')

            arch.vertex_attribute(key, 'ub', zub)
            arch.vertex_attribute(key, 'lb', zlb)
            arch.vertex_attribute(key, 'target', target)
            arch.vertex_attribute(key, 'pz', pz)

            print('copied {} to {}', key, pkey, pt_arch, pt_pattern)
            break

pzt = 0
for key in arch.vertices():
    pzi = arch.vertex_attribute(key, 'pz')
    pzt += pzi

print('Self weight in the arch is:', pzt)

# --------------------- 3. Create Starting point with TNA ---------------------

initialize_loadpath(arch)

view = Viewer(arch, show_grid=False)
view.draw_thrust()
view.show()

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
optimiser.settings['support_displacement'] = dXb

# --------------------- 5. Set up and run analysis ---------------------

analysis = Analysis.from_elements(vault, arch, optimiser)
# analysis.apply_selfweight()
# analysis.apply_envelope()
# analysis.apply_reaction_bounds()
analysis.set_up_optimiser()
analysis.run()

fopt = optimiser.fopt
print('Thrust/weight:', fopt/pzt)

# plot_form(form, show_q=False, cracks=True).show()

folder = os.path.join(folder, 'revision')

if optimiser.exitflag == 0:
    os.makedirs(os.path.join(folder, diagram_name), exist_ok=True)
    file_solution = os.path.join(folder, diagram_name, file_name + '_' + diagram_name + '_' + objective + '.json')
    if fill_loads:
        file_solution = os.path.join(folder, diagram_name, file_name + '_' + diagram_name + '_with_fill_' + objective + '.json')

    form.to_json(file_solution)

    print('Saved solution to:', file_solution)

viewer = Viewer(arch, shape=vault, show_grid=False)
viewer.settings['camera.target'] = [2.8, 1.6, 12]
viewer.settings['camera.distance'] = 18
viewer.settings['scale.reactions'] = 0.005 * 4
viewer.settings['opacity.shapes'] =  0.3
viewer.draw_thrust()
viewer.draw_cracks()
viewer.draw_shape()
viewer.draw_reactions()
viewer.show()
