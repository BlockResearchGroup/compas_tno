from numpy import arange
import math
from compas.datastructures import Mesh
from numpy import array

def domical_vault(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.5, radius=None, center=None, tol=10e-6, t=0.0, discretisation=[100, 100]):
    """ Set domical vault upper and lower meshes.

    Parameters
    ----------
    xy_span : list
        List with initial- and end-points of the vault [(x0,x1),(y0,y1)].
    thk : float (optional)
        Thickness of the vault - perpendicular to the middle surface.
    r : float (optional)
        Radius of the middle surface of the vault. If None, it will be calculated as half diagonal
    center: list
        List with coordinates (x,y) of the center of the sphere [xc,yc]
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
    middle = Mesh.from_vertices_and_faces(xyz, faces_i)

    xyz = array([x1d, y1d, zi1d]).transpose()
    intrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    xyz = array([x1d, y1d, ze1d]).transpose()
    extrados = Mesh.from_vertices_and_faces(xyz, faces_i)

    return intrados, extrados, middle



def paraboloid_vault(xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=None, tol=10e-6, t=10.0, discretisation=[100, 100]):


    return
