
import math
from compas.utilities import geometric_key
from compas.geometry import mirror_points_line
from compas.datastructures import Mesh


def create_cross_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=10, fix='corners'):
    """ Helper to construct a FormDiagram based on cross discretiastion with orthogonal arrangement and diagonal.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    discretisation: int
        Set the density of the grid in x and y directions.

    fix : string
        Option to select the constrained nodes: 'corners', 'all' are accepted.

    Returns
    -------
    obj
        FormDiagram.

    """

    if isinstance(discretisation, list):
        discretisation = discretisation[0]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]
    x_span = x1 - x0
    y_span = y1 - y0
    dx = x_span/discretisation
    dy = y_span/discretisation

    lines = []

    for i in range(discretisation+1):
        for j in range(discretisation+1):
            if i < discretisation and j < discretisation:
                # Vertical Members:
                xa = x0 + dx*i
                ya = y0 + dy*j
                xb = x0 + dx*(i + 1)
                yb = y0 + dy*j
                # Horizontal Members:
                xc = x0 + dx*i
                yc = y0 + dy*j
                xd = x0 + dx*i
                yd = y0 + dy*(j + 1)
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])
                lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])
                if i == j:
                    # Diagonal Members in + Direction:
                    xc = x0 + dx*i
                    yc = y0 + dy*j
                    xd = x0 + dx*(i + 1)
                    yd = y0 + dy*(j + 1)
                    lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])
                if i + j == discretisation:
                    # Diagonal Members in - Direction:
                    xc = x0 + dx*i
                    yc = y0 + dy*j
                    xd = x0 + dx*(i - 1)
                    yd = y0 + dy*(j + 1)
                    lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])
                    if i == (discretisation - 1):
                        xc = x0 + dx*i
                        yc = y0 + dy*j
                        xd = x0 + dx*(i + 1)
                        yd = y0 + dy*(j - 1)
                        lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])
            else:
                if i == discretisation and j < discretisation:
                    # Vertical Members on last column:
                    xa = x0 + dx*j
                    ya = y0 + dy*i
                    xb = x0 + dx*(j + 1)
                    yb = y0 + dy*i
                    # Horizontal Members:
                    xc = x0 + dx*i
                    yc = y0 + dy*j
                    xd = x0 + dx*i
                    yd = y0 + dy*(j + 1)
                    lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])
                    lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])

    mesh = Mesh.from_lines(lines, delete_boundary_face=True)
    form = cls.from_mesh(mesh)
    # form = cls.from_lines(lines, delete_boundary_face=True)
    gkey_key = form.gkey_key()

    if fix == 'corners':
        form.vertex_attribute(gkey_key[geometric_key([x0, y0, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x0, y1, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x1, y0, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x1, y1, 0.0])], 'is_fixed', True)
    else:
        [bnds] = form.vertices_on_boundaries()
        for key in bnds:
            form.vertex_attribute(key, 'is_fixed', True)

    return form


def create_fan_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10, 10], fix='corners'):
    """ Helper to construct a FormDiagram based on fan discretiastion with straight lines to the corners.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    discretisation: int
        Set the density of the grid in x and y directions.

    fix : string
        Option to select the constrained nodes: 'corners', 'all' are accepted.

    Returns
    -------
    obj
        FormDiagram.

    """

    if isinstance(discretisation, int):
        discretisation = [discretisation, discretisation]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    x_span = x1 - x0
    y_span = y1 - y0
    xc0 = x0 + x_span/2
    yc0 = y0 + y_span/2
    division_x = discretisation[0]
    division_y = discretisation[1]
    dx = float(x_span/division_x)
    dy = float(y_span/division_y)
    nx = int(division_x/2)
    ny = int(division_y/2)
    line_hor = [[x0, yc0, 0.0], [xc0, yc0, 0.0]]
    line_ver = [[xc0, y0, 0.0], [xc0, yc0, 0.0]]

    lines = []

    for i in range(nx):
        for j in range(ny+1):
            # Diagonal Members:
            xa = x0 + dx * i
            ya = y0 + dy * j * i / nx
            xb = x0 + dx * (i + 1)
            yb = y0 + dy * j * (i + 1) / nx
            lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

            a_mirror, b_mirror = mirror_points_line([[xa, ya, 0.0], [xb, yb, 0.0]], line_hor)
            lines.append([a_mirror, b_mirror])
            a_mirror, b_mirror = mirror_points_line([a_mirror, b_mirror], line_ver)
            lines.append([a_mirror, b_mirror])
            a_mirror, b_mirror = mirror_points_line([[xa, ya, 0.0],[xb, yb, 0.0]], line_ver)
            lines.append([a_mirror, b_mirror])

            xa_ = x0 + dy * j * i / nx
            ya_ = y0 + dx * i
            xb_ = x0 + dy * j * (i + 1) / nx
            yb_ = y0 + dx * (i + 1)
            lines.append([[xa_, ya_, 0.0], [xb_, yb_, 0.0]])

            a_mirror, b_mirror = mirror_points_line([[xa_, ya_, 0.0], [xb_, yb_, 0.0]], line_hor)
            lines.append([a_mirror, b_mirror])
            a_mirror, b_mirror = mirror_points_line([a_mirror, b_mirror], line_ver)
            lines.append([a_mirror, b_mirror])
            a_mirror, b_mirror = mirror_points_line([[xa_, ya_, 0.0], [xb_, yb_, 0.0]], line_ver)
            lines.append([a_mirror, b_mirror])

            if j < ny:
                # Vertical or Horizontal Members:
                xc = x0 + dx * (i + 1)
                yc = y0 + dy * j * (i + 1) / nx
                xd = x0 + dx * (i + 1)
                yd = y0 + dy * (j + 1) * (i + 1) / nx
                lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])

                c_mirror, d_mirror = mirror_points_line([[xc, yc, 0.0], [xd, yd, 0.0]], line_hor)
                lines.append([c_mirror, d_mirror])
                c_mirror, d_mirror = mirror_points_line([c_mirror, d_mirror], line_ver)
                lines.append([c_mirror, d_mirror])
                c_mirror, d_mirror = mirror_points_line([[xc, yc, 0.0], [xd, yd, 0.0]], line_ver)
                lines.append([c_mirror, d_mirror])

                xc_ = x0 + dy * j * (i + 1) / nx
                yc_ = y0 + dx * (i + 1)
                xd_ = x0 + dy * (j + 1) * (i + 1) / nx
                yd_ = y0 + dx * (i + 1)
                lines.append([[xc_, yc_, 0.0], [xd_, yd_, 0.0]])

                c_mirror, d_mirror = mirror_points_line([[xc_, yc_, 0.0], [xd_, yd_, 0.0]], line_hor)
                lines.append([c_mirror, d_mirror])
                c_mirror, d_mirror = mirror_points_line([c_mirror, d_mirror], line_ver)
                lines.append([c_mirror, d_mirror])
                c_mirror, d_mirror = mirror_points_line([[xc_, yc_, 0.0], [xd_, yd_, 0.0]], line_ver)
                lines.append([c_mirror, d_mirror])

    form = cls.from_lines(lines)

    gkey_key = form.gkey_key()

    if fix == 'corners':
        form.vertex_attribute(gkey_key[geometric_key([x0, y0, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x0, y1, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x1, y0, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x1, y1, 0.0])], 'is_fixed', True)
    else:
        [bnds] = form.vertices_on_boundaries()
        for key in bnds:
            form.vertex_attribute(key, 'is_fixed', True)

    return form


def create_ortho_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10,10], fix='corners'):
    """ Helper to construct a FormDiagram based on cross discretiastion with orthogonal arrangement and diagonal.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    discretisation: int
        Set the density of the grid in x and y directions.

    fix : string
        Option to select the constrained nodes: 'corners', 'all' are accepted.

    Returns
    -------
    obj
        FormDiagram.

    """

    if isinstance(discretisation, int):
        discretisation = [discretisation, discretisation]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]
    x_span = x1 - x0
    y_span = y1 - y0
    division_x = discretisation[0]
    division_y = discretisation[1]
    dx = float(x_span/division_x)
    dy = float(y_span/division_y)

    lines = []
    points = []
    pt_i = 0

    for i in range(division_x+1):
        for j in range(division_y+1):
            if i < division_x and j < division_y:
                # Vertical Members:
                xa = x0 + dx*i
                ya = y0 + dy*j
                xb = x0 + dx*(i + 1)
                yb = y0 + dy*j
                # Horizontal Members:
                xc = x0 + dx*i
                yc = y0 + dy*j
                xd = x0 + dx*i
                yd = y0 + dy*(j + 1)
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])
                lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])
            elif i == division_x and j == division_y - 1:
                # Horizontal Members:
                xc = x0 + dx*i
                yc = y0 + dy*j
                xd = x0 + dx*i
                yd = y0 + dy*(j + 1)
                lines.append([[xc, yc, 0.0], [xd, yd, 0.0]])
            elif j == division_y and i == division_x - 1:
                # Vertical Members:
                xa = x0 + dx*(i - 1)
                ya = y0 + dy*j
                xb = x0 + dx*i
                yb = y0 + dy*j
                lines.append([[xa, ya, 0.0], [xb, yb, 0.0]])

    # mesh = Mesh.from_lines(lines)
    # # mesh.plot()
    # form = cls.from_mesh(mesh)
    # form.plot()
    form = cls.from_lines(lines, delete_boundary_face=True)
    gkey_key = form.gkey_key()

    if fix == 'corners':
        form.vertex_attribute(gkey_key[geometric_key([x0, y0, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x0, y1, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x1, y0, 0.0])], 'is_fixed', True)
        form.vertex_attribute(gkey_key[geometric_key([x1, y1, 0.0])], 'is_fixed', True)
    else:
        [bnds] = form.vertices_on_boundaries()
        for key in bnds:
            form.vertex_attribute(key, 'is_fixed', True)

    return form
