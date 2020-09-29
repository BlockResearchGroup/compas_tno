from numpy import arange
import math
from compas.datastructures import Mesh
from numpy import array
from numpy import ones
from numpy import zeros


def cross_vault_highfields(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=None, tol=10e-6, t=10.0, discretisation=[100, 100]):
    """ Set Cross-Vault heights.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface

    tol : float (optional)
        Approximates the equations avoiding negative square-roots.

    discretisation: list
        Density of the grid that approximates the surfaces in x- and y- directions.

    t: float
        Negative lower-bound for the reactions position.

    Returns
    -------
    obj
        interp2d.

    Notes
    ----------------------
    Position of the quadrants is as in the schema below:

        Q3
    Q2      Q1
        Q4

    """

    if isinstance(discretisation, int):
        discretisation = [discretisation, discretisation]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    dx = x1 - x0
    dy = y1 - y0
    density_x = discretisation[0]
    density_y = discretisation[1]
    x = arange(x0, x1 + dx/density_x, dx/density_x)
    y = arange(y0, y1 + dy/density_y, dy/density_y)

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2
    hc = max(rx, ry)
    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]
            xd = x0 + (x1 - x0)/(y1 - y0) * (yi - y0)
            yd = y0 + (y1 - y0)/(x1 - x0) * (xi - x0)
            hxd = (math.sqrt((rx)**2 - ((xd - x0) - rx)**2))
            hyd = (math.sqrt((ry)**2 - ((yd - y0) - ry)**2))
            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
                z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
                z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
                z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
                z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
            else:
                print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi, yi))
                z = 0.0

            x1d.append(xi)
            y1d.append(yi)
            z1d.append(z)

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                face = [p1, p2, p4, p3]
                faces.append(face)
            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    middle = Mesh.from_vertices_and_faces(xyz, faces_i)

    extrados = cross_vault_highfields_ub(xy_span=xy_span, thk=thk, tol=tol, discretisation=discretisation)
    intrados = cross_vault_highfields_lb(xy_span=xy_span, thk=thk, tol=tol, t=t, discretisation=discretisation)

    return intrados, extrados, middle


def cross_vault_highfields_ub(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=None, tol=10e-6, discretisation=[100, 100]):
    """ Helper function to set the extrados of a parametric Cross-Vault.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    thk : float
        Thickness of the vault - perpendicular to the middle surface

    tol : float
        Approximates the equations avoiding negative square-roots.

    discretisation: int
        Density of the grid that approximates the surfaces

    Returns
    -------
    obj
        interp2d.

    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]
    dx = x1 - x0
    dy = y1 - y0
    density_x = discretisation[0]
    density_y = discretisation[1]
    x = arange(x0, x1 + dx/density_x, dx/density_x)
    y = arange(y0, y1 + dy/density_y, dy/density_y)

    y1 = xy_span[1][1] + thk/2
    y0 = xy_span[1][0] - thk/2
    x1 = xy_span[0][1] + thk/2
    x0 = xy_span[0][0] - thk/2
    rx = (x1 - x0)/2
    ry = (y1 - y0)/2
    hc = max(rx, ry)
    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]
            xd = x0 + (x1 - x0)/(y1 - y0) * (yi - y0)
            yd = y0 + (y1 - y0)/(x1 - x0) * (xi - x0)
            hxd = (math.sqrt((rx)**2 - ((xd - x0) - rx)**2))
            hyd = (math.sqrt((ry)**2 - ((yd - y0) - ry)**2))
            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
                z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
                z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
                z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
                z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
            else:
                print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi, yi))
                z = 0.0

            x1d.append(xi)
            y1d.append(yi)
            z1d.append(z)

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                face = [p1, p2, p4, p3]
                faces.append(face)

            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    extrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return extrados


def cross_vault_highfields_lb(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=None, tol=10e-6, t=10.0, discretisation=[100, 100]):
    """ Helper function to set the intrados of Cross-Vaults.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface

    tol : float (optional)
        Approximates the equations avoiding negative square-roots.

    discretisation: int
        Density of the grid that approximates the surfaces

    t: float
        Negative lower-bound for the reactions position.

    Returns
    -------
    obj
        interp2d.

    """

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]
    dx = x1 - x0
    dy = y1 - y0
    density_x = discretisation[0]
    density_y = discretisation[1]
    x = arange(x0, x1 + dx/density_x, dx/density_x)
    y = arange(y0, y1 + dy/density_y, dy/density_y)

    y1 = xy_span[1][1] - thk/2
    y0 = xy_span[1][0] + thk/2
    x1 = xy_span[0][1] - thk/2
    x0 = xy_span[0][0] + thk/2

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2
    hc = max(rx, ry)
    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]

            x1d.append(xi)
            y1d.append(yi)

            if ((yi) > y1 and ((xi) > x1 or (xi) < x0)) or ((yi) < y0 and ((xi) > x1 or (xi) < x0)):
                z = - 1*t
            else:
                if yi > y1:
                    yi = y1
                elif yi < y0:
                    yi = y0
                elif xi > x1:
                    xi = x1
                elif xi < x0:
                    xi = x0
                xd = x0 + (x1 - x0)/(y1 - y0) * (yi - y0)
                yd = y0 + (y1 - y0)/(x1 - x0) * (xi - x0)
                hxd = math.sqrt(abs((rx)**2 - ((xd - x0) - rx)**2))
                hyd = math.sqrt(abs((ry)**2 - ((yd - y0) - ry)**2))
                if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
                    z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
                elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
                    z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
                elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
                    z = hc*(hxd + math.sqrt((ry)**2 - ((yi - y0) - ry)**2))/(rx + ry)
                elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
                    z = hc*(hyd + math.sqrt((rx)**2 - ((xi - x0) - rx)**2))/(rx + ry)
                else:
                    print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi, yi))
                    z = 0.0

            z1d.append(z)

            if i < len(x) - 1 and j < len(y) - 1:
                p1 = (i, j)
                p2 = (i, j+1)
                p3 = (i+1, j)
                p4 = (i+1, j+1)
                face = [p1, p2, p4, p3]
                faces.append(face)

            index = index + 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    intrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return intrados


def crossvault_ub_lb_update(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], tol=1e-6):  # Only for Square Cross Vault

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    y1_ub = y1 + thk/2
    y0_ub = y0 - thk/2
    x1_ub = x1 + thk/2
    x0_ub = x0 - thk/2

    y1_lb = y1 - thk/2
    y0_lb = y0 + thk/2
    x1_lb = x1 - thk/2
    x0_lb = x0 + thk/2

    rx_ub = (x1_ub - x0_ub)/2
    ry_ub = (y1_ub - y0_ub)/2
    rx_lb = (x1_lb - x0_lb)/2
    ry_lb = (y1_lb - y0_lb)/2

    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t

    for i in range(len(x)):
        xi, yi = x[i], y[i]
        intrados_null = False
        if (yi > y1_lb and (xi > x1_lb or xi < x0_lb)) or (yi < y0_lb and (xi > x1_lb or xi < x0_lb)):
            intrados_null = True
        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub) - ry_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb) - ry_lb)**2)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub) - rx_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb) - rx_lb)**2)
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub) - ry_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb) - ry_lb)**2)
        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub) - rx_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb) - rx_lb)**2)
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return ub, lb


def crossvault_dub_dlb(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], tol=1e-6):  # Only for Square Cross Vault

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    y1_ub = y1 + thk/2
    y0_ub = y0 - thk/2
    x1_ub = x1 + thk/2
    x0_ub = x0 - thk/2

    y1_lb = y1 - thk/2
    y0_lb = y0 + thk/2
    x1_lb = x1 - thk/2
    x0_lb = x0 + thk/2

    rx_ub = (x1_ub - x0_ub)/2
    ry_ub = (y1_ub - y0_ub)/2
    rx_lb = (x1_lb - x0_lb)/2
    ry_lb = (y1_lb - y0_lb)/2

    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t
    dub = zeros((len(x), 1))
    dlb = zeros((len(x), 1))

    for i in range(len(x)):
        xi, yi = x[i], y[i]
        intrados_null = False
        if (yi > y1_lb and (xi > x1_lb or xi < x0_lb)) or (yi < y0_lb and (xi > x1_lb or xi < x0_lb)):
            intrados_null = True
        if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q1
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub) - ry_ub)**2)
            dub[i] = 1/2 * ry_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb) - ry_lb)**2)
                dlb[i] = - 1/2 * ry_lb/lb[i]
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol:  # Q3
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub) - rx_ub)**2)
            dub[i] = 1/2 * rx_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb) - rx_lb)**2)
                dlb[i] = - 1/2 * rx_lb/lb[i]
        elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q2
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub) - ry_ub)**2)
            dub[i] = 1/2 * ry_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb) - ry_lb)**2)
                dlb[i] = - 1/2 * ry_lb/lb[i]
        elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol:  # Q4
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub) - rx_ub)**2)
            dub[i] = 1/2 * rx_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb) - rx_lb)**2)
                dlb[i] = - 1/2 * rx_lb/lb[i]
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return dub, dlb  # ub, lb

