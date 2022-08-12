from compas_tno.shapes import Shape
from compas_tno.viewers import Viewer
from compas_tno.diagrams import FormDiagram
from compas_tno.analysis import Analysis
from compas_tno.optimisers import Optimiser
from compas_tno.plotters import TNOPlotter
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
# from compas_tno.utilities import apply_envelope_from_shape
sols = {}

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

for j in [1, 2, 3, 4]:  # range(6):

    delta = 1.0
    span = 10.0
    xspan = yspan = [0.0, span]
    xspan_vault = yspan_vault = [- delta, span + delta]
    thk = 0.5

    data_vault = {
        'type': 'crossvault',
        'xy_span': [xspan_vault, yspan_vault],
        'thk': thk,
        'discretisation': 10,
        't': 0.0
    }

    data_diagram = {
        'type': 'cross_fd',
        'xy_span': [xspan, yspan],
        'thk': thk,
        'discretisation': 10,
        'fix': 'corners'
    }

    vault = Shape.from_library(data_vault)
    form = FormDiagram.from_library(data_diagram)

    # from compas_tno.algorithms import apply_sag
    # apply_sag(form, boundary_force=25.0)

    objective = 'max_load'  # try 'max' 'Ecomp-linear'
    solver = 'IPOPT'  # try SLSQP
    constraints = ['funicular', 'envelope']
    variables = ['q', 'zb', 'lambdv']  # 't'
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
    optimiser.settings['max_iter'] = 5000
    optimiser.settings['gradient'] = True
    optimiser.settings['jacobian'] = True
    optimiser.settings['printout'] = True
    optimiser.settings['starting_point'] = starting_point
    optimiser.settings['sym_loads'] = False

    loaded_node = 59 - j

    # adding a line to the supports

    supports = [0, 109]

    form = form_add_lines_support(form, loaded_node, supports)

    max_load_mult = 2000.0
    n = form.number_of_vertices()
    pzv = zeros((n, 1))
    pzv[loaded_node] = -1.0
    print('\n----- Loading Node:', loaded_node)

    optimiser.settings['max_lambd'] = max_load_mult
    optimiser.settings['load_direction'] = pzv

    plotter = Plotter()
    plotter.fontsize = 6
    artist = plotter.add(form)
    # for line in support_lines:
    #     plotter.add(Line(line[0], line[1]), zorder=100000, width=5.0, draw_as_segment=True)

    # for line in new_lines:
    #     plotter.add(Line(line[0], line[1]), zorder=100000, width=5.0, draw_as_segment=True)
    # for pt in points:
    #     pt = Point(pt[0], pt[1], pt[2])
    #     plotter.add(pt, zorder=100001, width=5.0)
    artist.draw_vertexlabels()
    artist.draw_facelabels()
    plotter.zoom_extents()
    plotter.show()

    plotter = TNOPlotter(form)
    plotter.draw_form(scale_width=False)
    plotter.draw_supports()
    plotter.show()

    analysis = Analysis.from_elements(vault, form, optimiser)
    analysis.apply_selfweight()
    analysis.apply_envelope()
    analysis.set_up_optimiser()

    vault0 = Shape.from_formdiagram_and_attributes(form)

    pzt = 0
    for key in form.vertices():
        pz = form.vertex_attribute(key, 'pz')
        pzt += pz
    print('Total load of:', pzt)

    analysis.run()

    plotter = TNOPlotter(form, vault)
    plotter.show_solution()

    viewer = Viewer(form)
    viewer.settings['camera.show.grid'] = False
    viewer.settings['camera.distance'] = 35
    viewer.draw_thrust()
    viewer.draw_cracks()
    viewer.draw_shape()

    load_applied = analysis.optimiser.fopt
    length = 2.0
    x, y, z = form.vertex_coordinates(loaded_node)
    z += length + 0.1
    arrow = Arrow([x, y, z], [0, 0, -length])
    viewer.app.add(arrow, color=(0, 0, 0))
    pct = load_applied/pzt
    print('Percentage of load:', pct)
    sols[j] = pct

    viewer.show()

    print(sols)
