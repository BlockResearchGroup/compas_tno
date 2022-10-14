from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import TNOPlotter
from compas_tno.algorithms import apply_sag
from compas.geometry import intersection_segment_segment_xy
from compas.geometry import distance_point_point_xy
import compas_tno
from compas.geometry import Line
from compas.geometry import Point
from compas_view2.shapes import Arrow
from compas.geometry import normalize_vector
from compas.geometry import norm_vector
from numpy import array
from numpy import zeros
from compas_plotters import Plotter
from compas.colors import Color
import math
# from compas_tno.utilities import apply_envelope_from_shape
solution = {}

def form_add_lines_support(form, loaded_node, supports):

    text = {}
    new_lines = []
    text[loaded_node] = loaded_node

    xp, yp, _ = form.vertex_coordinates(loaded_node)
    fixed_coords = [form.vertex_coordinates(vertex) for vertex in form.vertices_where({'is_fixed': True})]

    lines = form.to_lines()
    support_lines = []
    points = []

    for support in supports:
        xs, ys, _ = form.vertex_coordinates(support)
        support_line = [[xs, ys, 0.0], [xp, yp, 0.0]]
        support_lines.append(support_line)

    for line in lines:
        int_pt1 = intersection_segment_segment_xy(support_lines[0], line)
        int_pt2 = intersection_segment_segment_xy(support_lines[1], line)
        if int_pt1 or int_pt2:
            int_pt = int_pt1 or int_pt2
            points.append(int_pt)
            if int_pt == line[0] or int_pt == line[1]:
                new_lines.append(line)  # necessary?
                pass
            else:
                cutline1 = [line[0], int_pt]
                cutline2 = [int_pt, line[1]]
                new_lines.append(cutline1)
                new_lines.append(cutline2)
        else:
            new_lines.append(line)

    for i in range(len(points) - 1):
        j = i + 1
        if distance_point_point_xy(points[i], points[j]) > 1e-3:
            new_lines.append([points[i], points[j]])

    form = FormDiagram.from_lines(new_lines)

    # find new vertices
    for vertex in form.vertices():
        coord = form.vertex_coordinates(vertex)
        dist = distance_point_point_xy(coord, [xp, yp])
        if dist < 1e-3:
            loaded_node = vertex
            text[loaded_node] = loaded_node
            continue
        for coord_fix in fixed_coords:
            dist = distance_point_point_xy(coord, coord_fix)
            if dist < 1e-3:
                form.vertex_attribute(vertex, 'is_fixed', True)

    return form


def find_loaded_node(form, load_ycoord):

    for key in form.vertices():
        pt = form.vertex_coordinates(key)
        if distance_point_point_xy(pt, [xc, load_ycoord, 0.0]) <1e-3:
            loaded_node = key
            break

    return loaded_node

solver = 'IPOPT'  # try SLSQP
starting_point = 'loadpath'
max_lambd = 200.0

discr = 14

sag = True
add_support = True
save = False

# folder = '/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/original/'
# folder = '/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/sag/'
# folder = '/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/sag_line/'
# folder = '/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/line/'
# folder = '/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/lambda=50/'
# folder = '/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/diagonals/'

# form = FormDiagram.create_cross_with_diagonal(discretisation=discr)

solutions = {}
count_points = int(discr/2 + 1)  # apply load to all points in a ridge
for j in range(count_points):
# for j in [4]:
# for j in [0]:

    solutions[j] = {}

    span = 10.0
    xc = yc = span/2
    spr_angle = 30.0
    alpha = 1/math.cos(math.radians(spr_angle))
    L = span * alpha
    Ldiff = L - span
    xyspan_shape = [[-Ldiff/2, span + Ldiff/2], [-Ldiff/2, span + Ldiff/2]]

    vault = Shape.create_crossvault(xyspan_shape)
    vault.ro = 1.0
    form = FormDiagram.create_cross_form(discretisation=discr)
    # form = FormDiagram.create_cross_with_diagonal(discretisation=discr)
    # form = FormDiagram.create_parametric_form(lambd=0.5)
    # form = FormDiagram.from_json('/Users/mricardo/compas_dev/me/max_load/crossvault/moving_load/mixed_form.json')
    type_fd = form.parameters['type']

    load_ycoord = yc - j/(count_points - 1) * yc
    loaded_node = find_loaded_node(form, load_ycoord)
    print('loaded node is:', loaded_node)

    apply_envelope_from_shape(form, vault)
    shape_original = Shape.from_formdiagram_and_attributes(form)

    if sag:
        apply_sag(form, boundary_force=25.0)
        load_ycoord = form.vertex_coordinates(loaded_node)[1]

    if add_support:
        supports = []
        for key in form.vertices_where({'is_fixed': True}):
            x, y, z = form.vertex_coordinates(key)
            if y < yc:
                supports.append(key)

        form = form_add_lines_support(form, loaded_node, supports)
        loaded_node = find_loaded_node(form, load_ycoord)  # new loaded node

    n = form.number_of_vertices()
    load_direction = zeros((n, 1))
    load_direction[loaded_node] = -1.0
    print('New Loaded Node:', loaded_node)

    load_file = folder + 'crossvault_maxload_{}_discr_{}_pos_{}.json'.format(type_fd, discr, j)
    form = FormDiagram.from_json(load_file)

    viewer = Viewer(form, shape_original)
    viewer.scale_edge_thickness(2.0)
    viewer.settings['camera.show.grid'] = False
    viewer.settings['camera.distance'] = 35
    viewer.draw_thrust()
    viewer.draw_cracks()
    viewer.draw_shape()

    length = 2.0
    x, y, z = form.vertex_coordinates(loaded_node)
    z += length + 0.1
    arrow = Arrow([x, y, z], [0, 0, -length])
    viewer.app.add(arrow, color=(0, 0, 0))

    viewer.show()
