

from numpy import zeros
from numpy import ones
from numpy import array
from numpy import linspace

from compas_tno.shapes import MeshDos
from compas_tno.shapes import rectangular_topology

import math


__all__ = ['pavillion_vault_highfields',
           'pavillionvault_ub_lb_update',
           'pavillionvault_dub_dlb',
           'pavillionvault_b_update',
           'pavillionvault_db'
           ]


def pavillion_vault_highfields(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=None, tol=10e-6, t=10.0, discretisation=[100, 100]):
    """ Set Pavillion-Vault heights.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface

    tol : float (optional)
        Approximates the equations avoiding negative square-roots.

    discretisation: int
        discretisation of the grid that approximates the surfaces

    t: float
        Negative lower-bound for the reactions position.

    Returns
    -------
    obj
        Mesh.

    Notes
    ----------------------
    Position of the quadrants is as in the schema below:

        Q3
    Q2      Q1
        Q4

    * W.I.P to update this function to work on rectangular pavillion vaults

    """

    # Uodate this function to work on rectangular Pavillion-Vaults

    if isinstance(discretisation, int):
        discretisation = [discretisation, discretisation]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    density_x = discretisation[0]
    density_y = discretisation[1]
    x = linspace(x0, x1, num=density_x+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
    y = linspace(y0, y1, num=density_y+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)

    xi, yi, faces_i = rectangular_topology(x, y)

    zt = pavillionvault_middle_update(xi, yi, t, xy_span=xy_span, tol=1e-6)
    xyzt = array([xi, yi, zt.flatten()]).transpose()
    middle = MeshDos.from_vertices_and_faces(xyzt, faces_i)

    zub, zlb = pavillionvault_ub_lb_update(xi, yi, thk, t, xy_span=xy_span, tol=1e-6)

    xyzub = array([xi, yi, zub.flatten()]).transpose()
    xyzlb = array([xi, yi, zlb.flatten()]).transpose()

    extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)
    intrados = MeshDos.from_vertices_and_faces(xyzlb, faces_i)

    return intrados, extrados, middle


def pavillionvault_middle_update(x, y, xy_span=[[0.0, 10.0], [0.0, 10.0]], tol=1e-6):

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2

    z = zeros((len(x), 1))

    for i in range(len(x)):
        xi, yi = x[i], y[i]
        if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q1
            z[i] = math.sqrt((ry)**2 - ((yi - y0)-ry)**2)
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q3
            z[i] = math.sqrt((ry)**2 - ((yi - y0)-ry)**2)
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q2
            z[i] = math.sqrt((rx)**2 - ((xi - x0)-rx)**2)
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q4
            z[i] = math.sqrt((rx)**2 - ((xi - x0)-rx)**2)
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return z


def pavillionvault_ub_lb_update(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], tol=1e-6):

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]

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
        if yi > y1_lb or xi > x1_lb or xi < x0_lb or yi < y0_lb:
            intrados_null = True
        if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q1
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub)-ry_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb)-ry_lb)**2)
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q3
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub)-ry_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb)-ry_lb)**2)
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q2
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub)-rx_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb)-rx_lb)**2)
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q4
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub)-rx_ub)**2)
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb)-rx_lb)**2)
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return ub, lb


def pavillionvault_dub_dlb(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], tol=1e-6):

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]

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
        if yi > y1_lb or xi > x1_lb or xi < x0_lb or yi < y0_lb:
            intrados_null = True
        if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q1
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub)-ry_ub)**2)
            dub[i] = 1/2 * ry_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb)-ry_lb)**2)
                dlb[i] = - 1/2 * ry_lb/lb[i]
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q3
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub)-ry_ub)**2)
            dub[i] = 1/2 * ry_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb)-ry_lb)**2)
                dlb[i] = - 1/2 * ry_lb/lb[i]
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q2
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub)-rx_ub)**2)
            dub[i] = 1/2 * rx_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb)-rx_lb)**2)
                dlb[i] = - 1/2 * rx_lb/lb[i]
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q4
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub)-rx_ub)**2)
            dub[i] = 1/2 * rx_ub/ub[i]
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb)-rx_lb)**2)
                dlb[i] = - 1/2 * rx_lb/lb[i]
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return dub, dlb  # ub, lb


def pavillionvault_b_update(x, y, thk, fixed, xy_span=[[0.0, 10.0], [0.0, 10.0]]):

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]
    b = zeros((len(fixed), 2))

    for i in range(len(fixed)):
        index = fixed[i]
        xi, yi = x[index], y[index]
        if xi == x0:
            b[[i], :] += [-thk/2, 0]
        elif xi == x1:
            b[i, :] += [+thk/2, 0]
        if yi == y0:
            b[i, :] += [0, -thk/2]
        elif yi == y1:
            b[i, :] += [0, +thk/2]

    return abs(b)


def pavillionvault_db(x, y, thk, fixed, xy_span=[[0.0, 10.0], [0.0, 10.0]]):

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]
    db = zeros((len(fixed), 2))

    for i in range(len(fixed)):
        index = fixed[i]
        xi, yi = x[index], y[index]
        if xi == x0:
            db[i, :] += [-1/2, 0]
        elif xi == x1:
            db[i, :] += [+1/2, 0]
        if yi == y0:
            db[i, :] += [0, -1/2]
        elif yi == y1:
            db[i, :] += [0, +1/2]

    return abs(db)
