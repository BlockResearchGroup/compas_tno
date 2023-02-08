

from numpy import zeros
from numpy import ones
from numpy import array
from numpy import linspace

from compas_tno.shapes import MeshDos
from compas_tno.shapes import rectangular_topology

import math


def pavillion_vault_highfields_proxy(xy_span, thk=0.5, discretisation=10, t=0.0):
    """Get Pavillionvault from proxy.

    Parameters
    ----------
    xy_span : [list, list], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        The thickness of the vault, by default 0.5
    tol : float, optional
        Tolerance, by default 10e-6
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0
    discretisation : list|int, optional
        Level of discretisation of the shape, by default [100, 100]

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape
    """

    intrados, extrados, middle = pavillion_vault_highfields(xy_span=xy_span, thk=thk, discretisation=discretisation, t=t)
    return intrados.to_data(), extrados.to_data(), middle.to_data()


def pavillion_vault_highfields(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=None, tol=10e-6, t=0.0, discretisation=[100, 100], spr_angle=0.0, expanded=False):
    """Set Pavillion vault heights

    Parameters
    ----------
    xy_span : [list, list], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        The thickness of the vault, by default 0.5
    tol : float, optional
        Tolerance, by default 10e-6
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0
    discretisation : list|int, optional
        Level of discretisation of the shape, by default [100, 100]
    spr_angle : float, optional
        Springing angle, by default 0.0
    expanded : bool, optional
        If the extrados should extend beyond the floor plan, by default False

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape
    """

    # TODO: Update this function to work on rectangular Pavillion-Vaults

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

    zt = pavillionvault_middle_update(xi, yi, xy_span=xy_span, spr_angle=spr_angle, tol=1e-6)
    xyzt = array([xi, yi, zt.flatten()]).transpose()
    middle = MeshDos.from_vertices_and_faces(xyzt, faces_i)

    zub, zlb = pavillionvault_ub_lb_update(xi, yi, thk, t, xy_span=xy_span, spr_angle=spr_angle, tol=1e-6)

    xyzub = array([xi, yi, zub.flatten()]).transpose()
    xyzlb = array([xi, yi, zlb.flatten()]).transpose()

    extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)
    intrados = MeshDos.from_vertices_and_faces(xyzlb, faces_i)

    if expanded:
        x = linspace(x0 - thk/2 / math.cos(math.radians(spr_angle)), x1 + thk/2 / math.cos(math.radians(spr_angle)),
                     num=density_x+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
        y = linspace(y0 - thk/2 / math.cos(math.radians(spr_angle)), y1 + thk/2 / math.cos(math.radians(spr_angle)),
                     num=density_y+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)
        xi, yi, faces_i = rectangular_topology(x, y)
        zub, zlb = pavillionvault_ub_lb_update(xi, yi, thk, t, xy_span=xy_span, spr_angle=spr_angle, tol=1e-6)
        xyzub = array([xi, yi, zub.flatten()]).transpose()
        extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)

    return intrados, extrados, middle


def pavillionvault_middle_update(x, y, xy_span=[[0.0, 10.0], [0.0, 10.0]], spr_angle=0.0, tol=1e-6):
    """Update middle of a pavillion vault based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the vault
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    spr_angle : float, optional
        Springing angle, by default 0.0
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    z : array
        Values of the middle surface in the points
    """

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]

    if spr_angle == 0.0:
        z_ = 0.0
    else:
        alpha = 1/math.cos(math.radians(spr_angle))
        z_ = (x1 - x0)/2 * math.tan(math.radians(spr_angle))
        L = x1 * alpha
        Ldiff = L - x1
        x0, x1 = -Ldiff/2, x1 + Ldiff/2
        y0, y1 = -Ldiff/2, y1 + Ldiff/2

    rx = (x1 - x0)/2
    ry = (y1 - y0)/2

    z = zeros((len(x), 1))

    for i in range(len(x)):
        xi, yi = x[i], y[i]
        if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q1
            z[i] = math.sqrt((ry)**2 - ((yi - y0)-ry)**2) - z_
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q3
            z[i] = math.sqrt((ry)**2 - ((yi - y0)-ry)**2) - z_
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q2
            z[i] = math.sqrt((rx)**2 - ((xi - x0)-rx)**2) - z_
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q4
            z[i] = math.sqrt((rx)**2 - ((xi - x0)-rx)**2) - z_
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return z


def pavillionvault_ub_lb_update(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], spr_angle=0.0, tol=1e-6):
    """Update upper and lower bounds of a pavillionvault based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the vault
    t : float
        Parameter for lower bound in nodes in the boundary
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """

    x0, x1 = xy_span[0]
    y0, y1 = xy_span[1]

    if spr_angle == 0.0:
        z_ = 0.0
    else:
        alpha = 1/math.cos(math.radians(spr_angle))
        z_ = (x1 - x0)/2 * math.tan(math.radians(spr_angle))
        L = x1 * alpha
        Ldiff = L - x1
        x0, x1 = -Ldiff/2, x1 + Ldiff/2
        y0, y1 = -Ldiff/2, y1 + Ldiff/2

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
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub)-ry_ub)**2) - z_
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb)-ry_lb)**2) - z_
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q3
            ub[i] = math.sqrt((ry_ub)**2 - ((yi - y0_ub)-ry_ub)**2) - z_
            if not intrados_null:
                lb[i] = math.sqrt((ry_lb)**2 - ((yi - y0_lb)-ry_lb)**2) - z_
        elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol:  # Q2
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub)-rx_ub)**2) - z_
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb)-rx_lb)**2) - z_
        elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol:  # Q4
            ub[i] = math.sqrt((rx_ub)**2 - ((xi - x0_ub)-rx_ub)**2) - z_
            if not intrados_null:
                lb[i] = math.sqrt((rx_lb)**2 - ((xi - x0_lb)-rx_lb)**2) - z_
        else:
            print('Error Q. (x,y) = ({0},{1})'.format(xi, yi))

    return ub, lb


def pavillionvault_dub_dlb(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], tol=1e-6):
    """Computes the sensitivities of upper and lower bounds in the x, y coordinates and thickness specified.

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the vault
    t : float
        Parameter for lower bound in nodes in the boundary
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    dub : array
        Values of the sensitivities for the upper bound in the points
    dlb : array
        Values of the sensitivities for the lower bound in the points
    """

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
    """Computes the ``b`` of parameter x, y coordinates and thickness specified.

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the vault
    fixed : list
        The list with indexes of the fixed vertices
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]

    Returns
    -------
    dub : array
        Values of the sensitivities for the upper bound in the points
    dlb : array
        Values of the sensitivities for the lower bound in the points
    """

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
    """Computes the sensitivities of the ``b`` parameter in the x, y coordinates and thickness specified.

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the vault
    fixed : list
        The list with indexes of the fixed vertices
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]

    Returns
    -------
    db : array
        Values of the sensitivities of the ``b`` parameter in the points
    """

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
