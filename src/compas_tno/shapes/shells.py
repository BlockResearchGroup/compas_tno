from numpy import arange
from numpy import array
from numpy import linspace

from compas_tno.shapes import MeshDos
from compas_tno.shapes import rectangular_topology

import math


def domical_vault(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.5, radius=None, center=None, tol=10e-6, t=0.0, discretisation=[100, 100]):
    """Set domical vault upper and lower meshes

    Parameters
    ----------
    xy_span : list, optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        [description], by default 0.5
    thk : float
        Thickness of the vault
    center : [float, float], optional
        Center of the dome in which the geometry of the vault is based, by default None in which the center point of the pattern is taken
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

    if isinstance(discretisation, int):
        discretisation = [discretisation, discretisation]

    y1 = xy_span[1][1]
    y0 = xy_span[1][0]
    x1 = xy_span[0][1]
    x0 = xy_span[0][0]

    dx = x1 - x0
    dy = y1 - y0
    d = math.sqrt(dx**2 + dy**2)
    r_diagonal = d/2

    if radius is None:
        radius = r_diagonal
    elif radius < r_diagonal:
        print('Check your radius!')
        raise Exception

    re = radius + thk/2
    ri = radius - thk/2

    if center is None:
        xc = (x1 + x0)/2
        yc = (y1 + y0)/2
    else:
        xc = center[0]
        yc = center[1]
        print('Center provided by the user:', center)

    density_x = discretisation[0]
    density_y = discretisation[1]
    x = arange(x0, x1 + dx/density_x, dx/density_x)
    y = arange(y0, y1 + dy/density_y, dy/density_y)

    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []
    zi1d = []
    ze1d = []

    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]
            z2 = radius**2 - (xi - xc)**2 - (yi - yc)**2
            zi2 = ri**2 - (xi - xc)**2 - (yi - yc)**2
            ze2 = re**2 - (xi - xc)**2 - (yi - yc)**2
            z = math.sqrt(z2)
            ze = math.sqrt(ze2)
            if zi2 < 0:
                zi = -t
            else:
                zi = math.sqrt(zi2)

            x1d.append(xi)
            y1d.append(yi)
            z1d.append(z)
            zi1d.append(zi)
            ze1d.append(ze)

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
    middle = MeshDos.from_vertices_and_faces(xyz, faces_i)

    xyz = array([x1d, y1d, zi1d]).transpose()
    intrados = MeshDos.from_vertices_and_faces(xyz, faces_i)

    xyz = array([x1d, y1d, ze1d]).transpose()
    extrados = MeshDos.from_vertices_and_faces(xyz, faces_i)

    return intrados, extrados, middle


def parabolic_shell_highfields(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.5, hc=5.0, tol=10e-6, t=0.0, discretisation=[100, 100]):
    """ Set parabolic vault heights.

    Parameters
    ----------
    xy_span : list, optional
        [description], by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        [description], by default 0.5
    discretisation : list, optional
        [description], by default [10, 10]
    hc : float, optional
        Height in the middle point of the vault, by default 5.0
    tol : float, optional
        Tolerance, by default 10e-6
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0
    discretisation : list, optional
        [description], by default [100, 100]

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape

    Notes
    -------------
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

    density_x = discretisation[0]
    density_y = discretisation[1]
    x = linspace(x0, x1, num=density_x+1, endpoint=True)  # arange(x0, x1 + dx/density_x, dx/density_x)
    y = linspace(y0, y1, num=density_y+1, endpoint=True)  # arange(y0, y1 + dy/density_y, dy/density_y)

    xi, yi, faces_i = rectangular_topology(x, y)

    zt = parabolic_shell_middle_update(xi, yi, t, xy_span=xy_span, hc=hc, tol=tol)
    xyzt = array([xi, yi, zt.flatten()]).transpose()
    middle = MeshDos.from_vertices_and_faces(xyzt, faces_i)

    zub, zlb = parabolic_shell_ub_lb_update(xi, yi, thk, t, xy_span=xy_span, hc=hc, tol=tol)
    xyzub = array([xi, yi, zub.flatten()]).transpose()
    xyzlb = array([xi, yi, zlb.flatten()]).transpose()
    extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)
    intrados = MeshDos.from_vertices_and_faces(xyzlb, faces_i)

    return intrados, extrados, middle


def parabolic_shell_middle_update(x, y, thk, xy_span=[[0.0, 10.0], [0.0, 10.0]], hc=5.0, tol=10e-6):
    """Update middle of a parabolic vault based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the arch
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    zt : array
        Values of the middle surface in the points
    """
    print('WIP')
    zt = 0
    return zt


def parabolic_shell_ub_lb_update(x, y, thk, t, xy_span=[[0.0, 10.0], [0.0, 10.0]], hc=5.0, tol=10e-6):
    """Update upper and lower bounds of a parabolic vault based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the arch
    t : float
        Parameter for lower bound in nodes in the boundary
    xy_span : [list[float], list[float]], optional
        xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
    hc : float, optional
        Height in the middle point of the vault, by default 8.0
    tol : float, optional
        Tolerance, by default 10e-6

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """
    print('WIP')
    ub = 0
    lb = 0
    return ub, lb
