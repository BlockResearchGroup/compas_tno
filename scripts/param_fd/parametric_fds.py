from multiprocessing.sharedctypes import Value
from compas.datastructures import Mesh
from compas_plotters import Plotter
from compas.geometry import Line
from compas.geometry import Point

from compas.geometry import mirror_points_line
from compas.geometry import rotate_points_xy

from compas.geometry import intersection_segment_segment_xy
from compas.geometry import distance_point_point_xy
from compas.geometry import closest_point_in_cloud_xy

from compas_tno.utilities import split_intersection_lines
from compas_tno.utilities.form import split_intersection_linesv0
from compas_tno.utilities import move_pattern_to_origin

xy_span = [[0, 10.0], [0, 8.0]]
xy_span = [[0.0, 5.744], [0.0, 3.344]]
discretisation = 10
lambd = 0.5


def create_parametric_form(xy_span, discretisation, lambd):

    if 0.0 > lambd or lambd > 1.0:
        raise ValueError('Lambda should be between 0.0 and 1.0')

    x_span = xy_span[0][1] - xy_span[0][0]
    y_span = xy_span[1][1] - xy_span[1][0]
    is_retangular = False

    if abs(x_span - y_span) > 1e-6:
        is_retangular = True
        y_span = x_span = 10.0
        x0, x1 = y0, y1 = 0.0, 10.0
    else:
        x0, x1 = xy_span[0][0], xy_span[0][1]
        y0, y1 = xy_span[1][0], xy_span[1][1]

    xc = (x1 + x0)/2
    yc = (y1 + y0)/2

    xc0 = x0 + x_span/2
    yc0 = y0 + y_span/2
    division_x = discretisation
    division_y = discretisation
    dx = float(x_span/division_x)
    dy = float(y_span/division_y)
    nx = int(division_x/2)
    line_hor = [[x0, yc0, 0.0], [xc0, yc0, 0.0]]
    line_ver = [[xc0, y0, 0.0], [xc0, yc0, 0.0]]

    def append_mirrored_lines(line, list_):
        mirror_a = mirror_points_line(line, line_hor)
        mirror_b = mirror_points_line(mirror_a, line_ver)
        mirror_c = mirror_points_line(line, line_ver)
        list_.append(mirror_a)
        list_.append(mirror_b)
        list_.append(mirror_c)

    lines = []

    for i in range(nx + 1):
        j = i

        xa = xc
        ya = xc - dy*j
        xb = (ya - y0) * (1 - lambd)
        yb = (ya - y0) * (1 - lambd)

        if distance_point_point_xy([xa, ya], [xb, yb]):
            lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            append_mirrored_lines([[xa, ya, 0.0], [xb, yb, 0.0]], lines)

        if i == 0:
            xa = x0
            ya = y0

            if distance_point_point_xy([xa, ya], [xb, yb]):
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

                append_mirrored_lines([[xa, ya, 0.0], [xb, yb, 0.0]], lines)

        xa = xc - dx*j
        xb = xc - dx*j
        ya = yc - dy*j
        yb = y0

        if distance_point_point_xy([xa, ya], [xb, yb]):
            lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            append_mirrored_lines([[xa, ya, 0.0], [xb, yb, 0.0]], lines)

        xa = yc - dy*j
        ya = yc
        xb = (xa - x0) * (1 - lambd)
        yb = xb * (y1 - y0) / (x1 - x0) + y0

        if distance_point_point_xy([xa, ya], [xb, yb]):
            lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            append_mirrored_lines([[xa, ya, 0.0], [xb, yb, 0.0]], lines)

        xa = xc - dx*j
        xb = x0
        ya = yc - dy*j
        yb = yc - dy*j

        if distance_point_point_xy([xa, ya], [xb, yb]):
            lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            append_mirrored_lines([[xa, ya, 0.0], [xb, yb, 0.0]], lines)

    clean_lines = split_intersection_linesv0(lines)

    mesh = Mesh.from_lines(clean_lines, delete_boundary_face=True)

    if is_retangular:
        x0, x1 = xy_span[0][0], xy_span[0][1]
        y0, y1 = xy_span[1][0], xy_span[1][1]
        move_pattern_to_origin(mesh, corners=[[x0, y0], [x0, y1], [x1, y1], [x0, y1]])

    return mesh


mesh = create_parametric_form(xy_span, discretisation, lambd)

plotter = Plotter()
artist = plotter.add(mesh)
# artist.draw_facelabels()
plotter.show()
