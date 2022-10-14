from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
import compas_tno
from compas.geometry import Line
from compas_view2.shapes import Arrow
from compas.geometry import normalize_vector
from compas.geometry import norm_vector
from compas.geometry import distance_point_point_xy
from numpy import array
from numpy import zeros
import math
# from compas_tno.utilities import apply_envelope_from_shape

thk = 0.50
spr_angle = 30
discr = 16
xf = 10.0
x0 = 0.0

xc = yc = (x0 + xf)/2
xyspan = [[x0, xf], [x0, xf]]
alpha = 1/math.cos(math.radians(spr_angle))
L = xf * alpha
Ldiff = L - xf
xyspan_shape = [[-Ldiff/2, xf + Ldiff/2], [-Ldiff/2, xf + Ldiff/2]]

form = FormDiagram.create_cross_form(xy_span=xyspan, discretisation=discr)
vault = Shape.create_crossvault(xy_span=xyspan_shape, discretisation=discr*2)

objective = 'min'  # try 'max' 'Ecomp-linear'
solver = 'SLSQP'  # try SLSQP
constraints = ['funicular', 'envelope']
variables = ['q', 'zb']  # 'lambdv', 't'
features = ['fixed']
axis_sym = None  # [[0.0, 5.0], [10.0, 5.0]]
starting_point = 'loadpath'

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

if objective == 'max_load':
    xp, yp = (5.0, 2.5)
    tol = 1e-3
    for key in form.vertices():
        xi, yi, _ = form.vertex_coordinates(key)
        if distance_point_point_xy([xi, yi], [xp, yp]) < tol:
            loaded_node = key
    max_load_mult = 2000.0
    n = form.number_of_vertices()
    pzv = zeros((n, 1))
    pzv[loaded_node] = -1.0

    optimiser.settings['max_lambd'] = max_load_mult
    optimiser.settings['load_direction'] = pzv

if objective == 'Ecomp-linear':
    Xc = [5.0, 5.0, 0.0]
    vector_supports = []
    sign = +1  # +1 for outwards / -1 for inwards

    for key in form.vertices_where({'is_fixed': True}):
        x, y, z = form.vertex_coordinates(key)
        xi, yi, zi = vault.intrados.vertex_coordinates(key)

        dXbi = normalize_vector([sign*(x - Xc[0]), sign*(y - Xc[1]), sign*(z - Xc[2])])  # 4 corners

        if x > Xc[0] and y < Xc[1]:             # vertical settlement
            # dXbi = normalize_vector([0, 0, -1])
            pass
        else:
            dXbi = [0, 0, 0]

        vector_supports.append(dXbi)
        form.vertex_attribute(key, 'dXb', dXbi)

    dXb = array(vector_supports)
    print(dXb)

    optimiser.settings['support_displacement'] = dXb

analysis = Analysis.from_elements(vault, form, optimiser)
analysis.apply_selfweight()
analysis.apply_envelope()
# analysis.set_up_optimiser()

vault0 = Shape.from_formdiagram_and_attributes(form)

# analysis.run()

# path = compas_tno.get('')
# form_path = path + '/form-' + objective + '.json'
# analysis_path = path + '/analysis-' + objective + '.json'
# form.to_json(form_path)
# print('Form saved to:', form_path)
# analysis.to_json(analysis_path)

# this is the solution for visualise the vault with a diagonal pull
form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-general1-Ecomp-linear-10-thk-0.5-corner-diagonal.json')

# this is the solution for visualise the
form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-direct_path-max_load-16-thk-0.5.json')

# form = FormDiagram.from_json('/Users/mricardo/compas_dev/compas_tno/data/CISM/form-appliedload.json')

vault2 = Shape.from_formdiagram_and_attributes(form)
viewer = Viewer(form, vault2)

viewer.draw_thrust()
# viewer.draw_middle_shape()

if objective == 't':
    # viewer.settings['color.mesh.intrados'] = (255, 150, 150)
    # viewer.settings['color.mesh.extrados'] = (255, 150, 150)
    viewer.draw_mesh(vault0.intrados, show_edges=False, color=(125, 125, 125))
    viewer.draw_mesh(vault0.extrados, show_edges=False, color=(125, 125, 125))
else:
    viewer.draw_cracks()

viewer.draw_shape()

if objective == 'max_load':
    length = 2.0
    x, y, z = form.vertex_coordinates(loaded_node)
    z += length + 0.1
    arrow = Arrow([x, y, z], [0, 0, -length])
    viewer.app.add(arrow, color=(0, 0, 0))

if objective == 'Ecomp-linear':
    for key in form.vertices_where({'is_fixed': True}):
        dXbi = form.vertex_attribute(key, 'dXb')
        x, y, _ = form.vertex_coordinates(key)
        z = form.vertex_attribute(key, 'lb') - 0.1
        if norm_vector(dXbi) > 0.01:
            arrow = Arrow([x, y, z], dXbi)
            viewer.app.add(arrow, color=(0, 0, 0))

viewer.show()

vault0 = Shape.from_formdiagram_and_attributes(form)
viewer = Viewer(form, vault0)

viewer.settings['scale.reactions'] = 0.005/8

viewer.draw_thrust()
viewer.draw_cracks()
viewer.draw_shape()
viewer.draw_reactions()

viewer.show()


points = []
supports = []
edges = []

for u, v in form.edges_where({'_is_edge': True}):
    Xu = form.vertex_coordinates(u)
    Xv = form.vertex_coordinates(v)

    # if abs(Xu[0] - xc) < 1e-3 and abs(Xv[0] - xc) < 1e-3:
    # if abs(Xu[1] - Xu[0]) < 1e-3 and abs(Xv[1] - Xv[0]) < 1e-3:
    if abs(Xu[1] + Xu[0] - 10.0) < 1e-3 and abs(Xv[1] + Xv[0] - 10) < 1e-3:
        edges.append((u, v))
        if u not in points:
            points.append(u)
        if v not in points:
            points.append(v)
        if form.vertex_attributes(u, 'is_fixed'):
            supports.append(u)
        if form.vertex_attributes(v, 'is_fixed'):
            supports.append(v)

# delete_faces = []
# for face in pavillion.intrados.faces():
#     Xf = pavillion.intrados.face_centroid(face)
#     if Xf[1] < yc:
#         delete_faces.append(face)

# for face in delete_faces:
#     pavillion.intrados.delete_face(face)
#     pavillion.extrados.delete_face(face)

view: Viewer = Viewer(form, shape=vault0)
view.settings['camera.show.grid'] = False
view.settings['camera.distance'] = 35
view.settings['camera.target'] = [5, 5, 0]
view.settings['camera.rz'] = 45
view.settings['camera.rx'] = 90.0
view.settings['camera.ry'] = 0.0
viewer.settings['scale.reactions'] = 0.005/3
view.draw_form(edges=edges, cull_negative=True)
view.draw_cracks(points=points, cull_negative=True)
# view.draw_reactions(supports=supports, extend_reactions=False)
view.draw_shape()
view.show()

