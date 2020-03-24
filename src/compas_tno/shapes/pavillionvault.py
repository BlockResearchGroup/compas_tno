

from scipy.interpolate import interp2d
from numpy import arange
from numpy import zeros
from numpy import array
from compas.datastructures import Mesh
import math


def pavillion_vault_highfields(xy_span = [[0.0,10.0],[0.0,10.0]], thk = None, tol = 10e-6, t = 10.0, discretisation = [100,100]):
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
            if yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q1
                z = math.sqrt((rx)**2 - (xi-rx)**2)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi >= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) - tol: #Q3
                z = math.sqrt((rx)**2 - (yi-rx)**2)
            elif yi >= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) - tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q2
                z = math.sqrt((ry)**2 - (xi-ry)**2)
            elif yi <= y0 + (y1 - y0)/(x1 - x0) * (xi - x0) + tol and yi <= y1 - (y1 - y0)/(x1 - x0) * (xi - x0) + tol: #Q4
                z = math.sqrt((ry)**2 - (yi-ry)**2)
            else:
                print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi,yi))
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
            u,v = uv
            i = uv_i[(u,v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    middle = Mesh.from_vertices_and_faces(xyz, faces_i)

    extrados = pavillion_vault_highfields_ub(xy_span=xy_span, thk=thk, tol=tol, discretisation=discretisation)
    intrados = pavillion_vault_highfields_lb(xy_span=xy_span, thk=thk, tol=tol, discretisation=discretisation, t=t)

    return intrados, extrados, middle

def pavillion_vault_highfields_ub(xy_span=[[0.0,10.0],[0.0,10.0]], thk = None, tol = 10e-6, discretisation = [100, 100]):

    """ Helper function to set the extrados of a parametric Pavillion-Vault.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

    thk : float
        Thickness of the vault - perpendicular to the middle surface

    tol : float
        Approximates the equations avoiding negative square-roots.

    discretisation: int
        discretisation of the grid that approximates the surfaces

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

    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] + thk)/2.0
    for i in range(len(x)):
        for j in range(len(y)):
            uv_i[(i, j)] = index
            xi, yi = x[i], y[j]
            if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q1
                z = math.sqrt((ry)**2 - ((yi - y0)-ry)**2)
            elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q3
                z = math.sqrt((ry)**2 - ((yi - y0)-ry)**2)
            elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q2
                z = math.sqrt((rx)**2 - ((xi - x0)-rx)**2)
            elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q4
                z = math.sqrt((rx)**2 - ((xi - x0)-rx)**2)
            else:
                print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi,yi))
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
            u,v = uv
            i = uv_i[(u,v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    extrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return extrados

def pavillion_vault_highfields_lb(xy_span = [[0.0,10.0],[0.0,10.0]], thk = None, tol = 10e-6, t = 10.0, discretisation = [100, 100]):

    """ Helper function to set the intrados of Pavillion-Vaults.

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

    index = 0
    uv_i = {}
    faces = []
    faces_i = []
    x1d = []
    y1d = []
    z1d = []

    if xy_span[0] == xy_span[1]:
        rx = ry = (xy_span[0][1] - xy_span[0][0] - thk)/2.0

    for i in range(len(x)):
        for j in range(len(y)):
            xi, yi = x[i], y[j]
            uv_i[(i, j)] = index
            if (yi) > y1 or (xi) > x1 or (xi) < x0 or (yi) < y0:
                z = -1.0 * t
            else:
                if (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q1
                    z = math.sqrt(abs((ry)**2 - ((yi - y0)-ry)**2))
                elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q3
                    z = math.sqrt(abs((ry)**2 - ((yi - y0)-ry)**2))
                elif (yi - y0) <= y1/x1 * (xi - x0) + tol and (yi - y0) >= (y1 - y0) - (xi - x0) - tol: #Q2
                    z = math.sqrt(abs((rx)**2 - ((xi - x0)-rx)**2))
                elif (yi - y0) >= y1/x1 * (xi - x0) - tol and (yi - y0) <= (y1 - y0) - (xi - x0) + tol: #Q4
                    z = math.sqrt(abs((rx)**2 - ((xi - x0)-rx)**2))
                else:
                    print('Vertex did not belong to any Q. (x,y) = ({0},{1})'.format(xi,yi))
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
            u,v = uv
            i = uv_i[(u,v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz = array([x1d, y1d, z1d]).transpose()
    intrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return intrados