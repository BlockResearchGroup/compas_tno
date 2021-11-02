import compas_tno
from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import plot_form

from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_envelope_on_xy

from compas_tno.plotters import plot_superimposed_diagrams

from compas_view2.shapes import Arrow
from numpy import array

from compas.geometry import norm_vector

from compas_tno.analysis import Analysis
from compas_tno.viewers import Viewer
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
ro_fill = 10.0
fill = False
c = 0.01

objectives = ['Ecomp-linear']
solver = 'SLSQP'
constraints = ['funicular', 'envelope', 'envelopexy']
variables = ['ind', 'zb']
features = []  # , 'symmetry']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
# qmax = 10e+6
starting_point = 'loadpath'
type_movement = 'corner_pull'  # corner_pull
gradients = True
max_iter = 500
save = False
sign = +1.0

# ------------------- FormDiagram --------------------

folder = '/Users/mricardo/compas_dev/me/anagni/'

file_formdiagram = os.path.join(folder, 'top_vault_form.json')

mesh = Mesh.from_json(file_formdiagram)
vertices, faces = mesh.to_vertices_and_faces()
form = FormDiagram.from_vertices_and_faces(vertices, faces)
left_supports = []

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
        if a:
            left_supports.append(key)

view = Viewer(form)
view.settings['size.edge.max_thickness'] = 15.0
vector_supports = []
lines = []

for key in form.vertices_where({'is_fixed': True}):
    x, y, z = form.vertex_coordinates(key)

    if key in left_supports:
        if type_movement == 'corner_settlement':
            dXbi = [0, 0, sign*-1]
        if type_movement == 'corner_pull':
            dXbi = [sign*-1, 0, 0]
    else:
        dXbi = [0, 0, 0]

    vector_supports.append(dXbi)
    lines.append({
        'start': [x, y, z],
        'end': [x + dXbi[0], y + dXbi[1], z + dXbi[2]],
        'width': 3
    })

    if norm_vector(dXbi) > 0.01:
        arrow = Arrow([x, y, z], dXbi)
        view.app.add(arrow, color=(0, 0, 0))

dXb = array(vector_supports)
print(dXb)

view.view_shape()
view.show()
view.clear()

# plotter = MeshPlotter(form, figsize=(5, 5))
# plotter.draw_edges()
# plotter.draw_vertices(keys=[key for key in form.vertices_where({'is_fixed': True})])
# plotter.draw_faces()
# plotter.show()

# ----------------------- Point Cloud -----------------------

file_name = 'sangelo_vault_top_final'  # This point cloud considers the extrados with 25 cm only
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

form_base = form

solutions = {}
i = 0

for objective in objectives:

    solutions[objective] = {}

    for n in [0.0]:  # [0.05, 0.06, 0.07, 0.08, 0.09]

        solutions[objective][n] = {}

        # ------- Create shape given a topology and a point cloud --------

        vault = Shape.from_pointcloud_and_topology(form, points_lb, points_ub, fill_pts=points_fill, data={'type': 'general', 't': 0.0, 'thk': thk})
        vault.store_normals(plot=False)

        area = vault.middle.area()
        swt = vault.compute_selfweight()

        print('Interpolated Volume Data:')
        print('Self-weight is: {0:.2f}'.format(swt))
        print('Area is: {0:.2f}'.format(area))

        if i == 0:
            apply_selfweight_from_shape(form, vault)
        i += 1

        for key in vault.intrados.vertices():
            normal = vault.intrados.vertex_attribute(key, 'n')
            nflip = [-normal[0], -normal[1], -normal[2]]
            vault.intrados.vertex_attribute(key, 'n', nflip)

        for key in vault.extrados.vertices():
            normal = vault.extrados.vertex_attribute(key, 'n')
            nflip = [-normal[0], -normal[1], -normal[2]]
            vault.extrados.vertex_attribute(key, 'n', nflip)

        apply_envelope_from_shape(form, vault)

        pzt0 = 0
        for key in form.vertices():
            pzi = form.vertex_attribute(key, 'pz')
            pzt0 += pzi

        print('Selfweight not considering the fill is:', pzt0)

        if 'envelopexy' in constraints:
            apply_envelope_on_xy(form, c=c)

        form_base = form.copy()




        if fill:

            for key in vault.fill.vertices():
                pz = form.vertex_attribute(key, 'pz')
                zub = form.vertex_attribute(key, 'ub')
                z_fill = vault.fill.vertex_attribute(key, 'z')
                z_diff = (z_fill - zub)
                if z_diff > 0.01:
                    ai = vault.middle.vertex_area(key)
                    pz_fill = -1 * ai * z_diff * ro_fill  # downwards loads are negative
                    pzt = pz + pz_fill
                    form.vertex_attribute(key, 'pz', pzt)

            pztfinal = 0
            for key in form.vertices():
                pzi = form.vertex_attribute(key, 'pz')
                pztfinal += pzi

            print('Selfweight with the considering the fill is:', pztfinal)
            print('Fill weight is:', pztfinal - pzt0)

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
        optimiser.settings['axis_symmetry'] = axis_sym
        optimiser.settings['support_displacement'] = dXb
        optimiser.settings['Ecomp_method'] = 'complete'

        # --------------------- 5. Set up and run analysis ---------------------

        analysis = Analysis.from_elements(vault, form, optimiser)
        analysis.set_up_optimiser()
        analysis.run()

        weight = 0
        for key in form.vertices():
            weight += form.vertex_attribute(key, 'pz')

        thrust = form.thrust()

        T_over_W = abs(thrust/weight)
        print('Ratio Thrust/Weight:', T_over_W)

        solutions[objective][n] = T_over_W

        # plot_form(form, show_q=False, cracks=True).show()

        folder = os.path.join(folder, 'compl-energy')
        os.makedirs(folder, exist_ok=True)

        if save:
            file_solution = os.path.join(folder, file_name + '_solution_form_' + type_movement + '_sign_' + str(sign) + '_' + objective + '.json')
            form.to_json(file_solution)

            file_shape = os.path.join(folder, file_name + '_solution_shape_type_' + type_movement + '_sign_' + str(sign) + '_' + objective + '.json')
            vault.to_json(file_shape)

        print(solutions)

        plot_superimposed_diagrams(form, form_base).show()

        view.show_solution()

print(solutions)
