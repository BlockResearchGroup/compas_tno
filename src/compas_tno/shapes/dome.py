from numpy import zeros
from numpy import array
from numpy import ones
import math

from compas_tno.shapes import MeshDos
from compas.datastructures import mesh_delete_duplicate_vertices


def dome_heightfields_proxy(center=[5.0, 5.0], radius=5.0, thk=0.30, t=5.0, discretisation=[8, 20], *args, **kwargs):
    """Function that computes the highfield of a dome through the proxy

    Parameters
    ----------
    center : [float, float], optional
        xy-span of the pattern, by default [[0.0, 10.0], [0.0, 10.0]]
    thk : float, optional
        The thickness of the vault, by default 0.5
    tol : float, optional
        Tolerance, by default 10e-6
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0
    discretisation : list|int, optional
        LEvel of discretisation of the shape, by default [100, 100]

    Returns
    -------
    intradosdata
        Data to a Mesh for the intrados of the pattern
    extradosdata
        Data to a Mesh for the extrados of the pattern
    middledata
        Data to a Mesh for the middle of the pattern
    """
    intrados, extrados, middle = set_dome_heighfield(center=center, radius=radius, thk=thk, t=t, discretisation=discretisation)
    return intrados.to_data(), extrados.to_data(), middle.to_data()


def set_dome_heighfield(center=[5.0, 5.0], radius=5.0, thk=0.30, t=0.0, discretisation=[8, 20], expanded=False):
    """Height of the dome heighfield

    Parameters
    ----------
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0
    thk : float, optional
        The thickness of the dome, by default 0.30
    t : float, optional
        Parameter to the lowerbound of the nodes in the boundary, by default 0.0
    discretisation : [float, float], optional
        Discretisation level for the dome (hoops and meridians), by default [8, 20]
    expanded : bool, optional
        If should expand 1 row to "close" the dome, by default False

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape
    """

    tol = 10e-3
    xc = center[0]
    yc = center[1]
    ri = radius - thk/2
    re = radius + thk/2
    n_radial = discretisation[0]
    n_spikes = discretisation[1]
    r_div = radius/n_radial
    theta = 2*math.pi/n_spikes

    x1d = []
    y1d = []
    zt1d = []
    zi1d = []
    ze1d = []
    faces = []
    faces_i = []
    i = 0
    uv_i = {}

    for nr in range(n_radial+1):
        for nc in range(n_spikes):
            uv_i[(nr, nc)] = i
            xi = xc + nr * r_div * math.cos(theta * nc)
            yi = yc + nr * r_div * math.sin(theta * nc)
            zt2 = radius**2 - (xi - xc)**2 - (yi - yc)**2
            zi2 = ri**2 - (xi - xc)**2 - (yi - yc)**2
            ze2 = re**2 - (xi - xc)**2 - (yi - yc)**2
            ze = math.sqrt(ze2)
            if zt2 < 0.0:
                zt = 0.0
            else:
                zt = math.sqrt(zt2)
            if zi2 < 0.0:
                if zi2 < -tol:
                    zi = -t
                else:
                    zi = 0.0
            else:
                zi = math.sqrt(zi2)

            x1d.append(xi)
            y1d.append(yi)
            zt1d.append(zt)
            zi1d.append(zi)
            ze1d.append(ze)

            if nr < n_radial:
                p1 = (nr, nc)
                p2 = (nr, nc+1)
                p3 = (nr+1, nc)
                p4 = (nr+1, nc+1)
                if nc < n_spikes - 1:  # General Case
                    # face = [p1, p2, p4, p3]
                    face = [p3, p4, p2, p1]
                else:
                    p2 = (nr, 0)
                    p4 = (nr+1, 0)
                    # face = [p1, p2, p4, p3]
                    face = [p3, p4, p2, p1]
                faces.append(face)

            i += 1

    if expanded:
        nr = n_radial + 1
        faces_exp = faces[:]
        x1d_exp = x1d[:]
        y1d_exp = y1d[:]
        for nc in range(n_spikes):
            uv_i[(nr, nc)] = i
            ze = -t
            zi = -t
            xi = xc + re * math.cos(theta * nc)
            yi = yc + re * math.sin(theta * nc)
            x1d_exp.append(xi)
            y1d_exp.append(yi)
            zi1d.append(zi)
            ze1d.append(ze)

            p1 = (nr-1, nc)
            p2 = (nr-1, nc+1)
            p3 = (nr, nc)
            p4 = (nr, nc+1)
            if nc < n_spikes - 1:  # General Case
                face = [p3, p4, p2, p1]
            else:
                p2 = (nr-1, 0)
                p4 = (nr, 0)
                face = [p3, p4, p2, p1]
            faces_exp.append(face)

            i += 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz_middle = array([x1d, y1d, zt1d]).transpose()
    middle = MeshDos.from_vertices_and_faces(xyz_middle, faces_i)

    if expanded:
        faces_i_exp = []
        for face in faces_exp:
            face_i = []
            for uv in face:
                u, v = uv
                i = uv_i[(u, v)]
                face_i.append(i)
            faces_i_exp.append(face_i)

        xyz_middle = array([x1d, y1d, zt1d]).transpose()
        xyz_intrados = array([x1d_exp, y1d_exp, zi1d]).transpose()
        xyz_extrados = array([x1d_exp, y1d_exp, ze1d]).transpose()
        middle = MeshDos.from_vertices_and_faces(xyz_middle, faces_i)
        intrados = MeshDos.from_vertices_and_faces(xyz_intrados, faces_i_exp)
        extrados = MeshDos.from_vertices_and_faces(xyz_extrados, faces_i_exp)
    else:
        xyz_middle = array([x1d, y1d, zt1d]).transpose()
        xyz_intrados = array([x1d, y1d, zi1d]).transpose()
        xyz_extrados = array([x1d, y1d, ze1d]).transpose()
        middle = MeshDos.from_vertices_and_faces(xyz_middle, faces_i)
        intrados = MeshDos.from_vertices_and_faces(xyz_intrados, faces_i)
        extrados = MeshDos.from_vertices_and_faces(xyz_extrados, faces_i)

    mesh_delete_duplicate_vertices(middle)  # this is a hack. Do it better
    mesh_delete_duplicate_vertices(intrados)
    mesh_delete_duplicate_vertices(extrados)

    return intrados, extrados, middle


def set_dome_with_spr(center=[5.0, 5.0], radius=5.0, thk=0.30, theta=[0, math.pi/2], t=5.0, discretisation=[8, 20]):
    """Height of the dome heighfield with spring angle

    Parameters
    ----------
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0
    thk : float, optional
        The thickness of the dome, by default 0.30
    theta : [float, float], optional
        The range to cover, by default [0, pi/2] which is a hemispheric dome
    t : float, optional
        Parameter to the lowerbound of the nodes in the boundary, by default 0.0
    discretisation : [float, float], optional
        Discretisation level for the dome (hoops and meridians), by default [8, 20]

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape
    """

    center.append(0.0)
    ri = radius - thk/2
    re = radius + thk/2
    phi_lower = 0
    phi_upper = 2 * math.pi
    phi_length = discretisation[1]
    phi_range = [phi_lower + x * (phi_upper - phi_lower) / phi_length for x in range(phi_length + 1)]

    theta_length = discretisation[0]
    theta_lower = theta[0]
    theta_upper = theta[1]
    theta_range = [theta_lower + x * (theta_upper - theta_lower) / theta_length for x in range(theta_length + 1)]

    [xmin, _, _] = geom_dome(center, radius, theta_range[0], phi_range[0])
    [xmax, _, _] = geom_dome(center, radius, theta_range[theta_length], phi_range[phi_length])
    theta_upper_i = (math.pi/2 - math.acos(min(1, (xmax-center[0])/ri)))
    theta_lower_e = (math.pi/2 - math.acos((xmin-center[0])/re))
    # print('min/max', xmin, xmax)
    # print('lower/lower_i/upper/upper_e', theta_lower/(math.pi/2), theta_lower_e/(math.pi/2), theta_upper/(math.pi/2), theta_upper_i/(math.pi/2))
    theta_range_i = [theta_lower + x * (theta_upper_i - theta_lower) / theta_length for x in range(theta_length + 1)]
    theta_range_e = [theta_lower_e + x * (theta_upper - theta_lower_e) / theta_length for x in range(theta_length + 1)]

    faces = []
    faces_i = []
    count = 0
    uv_i = {}
    xyz_middle = []
    xyz_intrados = []
    xyz_extrados = []

    for i in range(phi_length + 1):
        for j in range(theta_length + 1):
            uv_i[(i, j)] = count
            p1 = (i, j)
            p2 = (i, j+1)
            p3 = (i+1, j)
            p4 = (i+1, j+1)
            face = []

            ptm = geom_dome(center, radius, theta_range[j], phi_range[i])
            pti = geom_dome(center, ri, theta_range_i[j], phi_range[i])
            pte = geom_dome(center, re, theta_range_e[j], phi_range[i])

            xyz_middle.append(ptm)
            xyz_intrados.append(pti)
            xyz_extrados.append(pte)

            if i < phi_length and j < theta_length:
                face = [p3, p4, p2, p1]
                faces.append(face)

            count += 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz_middle = array(xyz_middle)
    xyz_extrados = array(xyz_extrados)
    xyz_intrados = array(xyz_intrados)

    middle = MeshDos.from_vertices_and_faces(xyz_middle, faces_i)
    extrados = MeshDos.from_vertices_and_faces(xyz_extrados, faces_i)
    intrados = MeshDos.from_vertices_and_faces(xyz_intrados, faces_i)

    mesh_delete_duplicate_vertices(middle)  # this is a hack. Do it better
    mesh_delete_duplicate_vertices(intrados)
    mesh_delete_duplicate_vertices(extrados)

    return intrados, extrados, middle


def set_dome_polar_coord(center=[5.0, 5.0], radius=5.0, thk=0.30, theta=[0, math.pi/2], t=0.0, discretisation=[8, 20]):
    """Height of the dome heighfield assuming polar coordinates

    Parameters
    ----------
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0
    thk : float, optional
        The thickness of the dome, by default 0.30
    theta : [float, float], optional
        The range to cover, by default [0, pi/2] which is a hemispheric dome
    t : float, optional
        Parameter to the lowerbound of the nodes in the boundary, by default 0.0
    discretisation : [float, float], optional
        Discretisation level for the dome (hoops and meridians), by default [8, 20]

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape
    """

    center0 = center
    center.append(0.0)
    ri = radius - thk/2
    re = radius + thk/2
    phi_lower = 0
    phi_upper = 2 * math.pi
    phi_length = discretisation[1]
    # step_phi = (phi_upper - phi_lower)
    phi_range = [phi_lower + x * (phi_upper - phi_lower) / phi_length for x in range(phi_length + 1)]

    theta_length = discretisation[0]
    theta_lower = theta[0]
    theta_upper = theta[1]
    theta_range = [theta_lower + x * (theta_upper - theta_lower) / theta_length for x in range(theta_length + 1)]

    [xmin, _, _] = geom_dome(center, radius, theta_range[0], phi_range[0])
    [xmax, _, _] = geom_dome(center, radius, theta_range[theta_length], phi_range[phi_length])
    theta_upper_i = (math.pi/2 - math.acos(min(1, (xmax-center[0])/ri)))
    theta_lower_e = (math.pi/2 - math.acos((xmin-center[0])/re))
    # print('min/max', xmin, xmax)
    # print('lower/lower_i/upper/upper_e', theta_lower/(math.pi/2), theta_lower_e/(math.pi/2), theta_upper/(math.pi/2), theta_upper_i/(math.pi/2))
    theta_range_i = [theta_lower + x * (theta_upper_i - theta_lower) / theta_length for x in range(theta_length + 1)]
    theta_range_e = [theta_lower_e + x * (theta_upper - theta_lower_e) / theta_length for x in range(theta_length + 1)]

    faces = []
    faces_i = []
    count = 0
    uv_i = {}
    xyz_middle = []
    xyz_intrados = []
    xyz_extrados = []

    for i in range(phi_length + 1):
        for j in range(theta_length + 2):
            uv_i[(i, j)] = count
            p1 = (i, j)
            p2 = (i, j+1)
            p3 = (i+1, j)
            p4 = (i+1, j+1)
            face = []

            if j <= theta_length:
                ptm = geom_dome(center, radius, theta_range[j], phi_range[i])
                pte = geom_dome(center, re, theta_range_e[j], phi_range[i])
                pti = geom_dome(center, ri, theta_range_i[j], phi_range[i])
            else:
                ptm = geom_dome(center, radius, theta_range[-1], phi_range[i])
                pte = geom_dome(center, re, theta_range_e[-1], phi_range[i])
                pti[:] = ptm
                pti[2] = -t

            xyz_middle.append(ptm)
            xyz_intrados.append(pti)
            xyz_extrados.append(pte)

            if i < phi_length and j < theta_length + 1:
                face = [p3, p4, p2, p1]
                faces.append(face)

            count += 1

    for face in faces:
        face_i = []
        for uv in face:
            u, v = uv
            i = uv_i[(u, v)]
            face_i.append(i)
        faces_i.append(face_i)

    xyz_middle = array(xyz_middle)
    xyz_extrados = array(xyz_extrados)
    xyz_intrados = array(xyz_intrados)

    middle = MeshDos.from_vertices_and_faces(xyz_middle, faces_i)
    extrados = MeshDos.from_vertices_and_faces(xyz_extrados, faces_i)
    intrados = MeshDos.from_vertices_and_faces(xyz_intrados, faces_i)

    for mesh in [middle, intrados, extrados]:
        mesh_delete_duplicate_vertices(mesh)  # this is a hack. Do it better

    for mesh in [middle, intrados, extrados]:
        delete_fkey = []
        for fkey in mesh.faces():
            if len(mesh.face_coordinates(fkey)) == 2:
                delete_fkey.append(fkey)
        for fkey in delete_fkey:
            mesh.delete_face(fkey)

    center = center0

    return intrados, extrados, middle


def geom_dome(p0, ro, theta, phi):
    """Geometry of the dome in polar coordinates

    Parameters
    ----------
    p0 : list
        center of the dome
    ro : float
        ro parameter
    theta : float
        theta parameter
    phi : float
        phi parameter

    Returns
    -------
    point : list
        List of cartesian coordinates
    """
    x = ro * math.sin(theta) * math.cos(phi)
    y = ro * math.sin(theta) * math.sin(phi)
    z = ro * math.cos(theta)
    point = [p0[0] + x, p0[1] + y, p0[2] + z]
    return point


def dome_zt_update(x, y, radius, t, center=[5.0, 5.0]):
    """Update middle of the dome based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    radius : float, optional
        The radius of the dome, by default 5.0
    t : float
        Parameter for lower bound in nodes in the boundary
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """

    xc = center[0]
    yc = center[1]
    zt = ones((len(x), 1))

    for i in range(len(x)):
        zt2 = radius**2 - (x[i] - xc)**2 - (y[i] - yc)**2
        if zt2 > 0:
            zt[i] = math.sqrt(zt2)
        else:
            zt[i] = 0.0

    return zt


def dome_ub_lb_update(x, y, thk, t, center=[5.0, 5.0], radius=5.0):
    """Update upper and lower bounds of the dome based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the dome
    t : float
        Parameter for lower bound in nodes in the boundary
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """

    xc = center[0]
    yc = center[1]
    ri = radius - thk/2
    re = radius + thk/2
    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * -t

    for i in range(len(x)):
        zi2 = ri**2 - (x[i] - xc)**2 - (y[i] - yc)**2
        ze2 = re**2 - (x[i] - xc)**2 - (y[i] - yc)**2
        ub[i] = math.sqrt(ze2)
        if zi2 > 0.0:
            lb[i] = math.sqrt(zi2)

    return ub, lb


def dome_dub_dlb(x, y, thk, t, center=[5.0, 5.0], radius=5.0):
    """Update sensitivites of upper and lower bounds of the dome based in the parameters

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the dome
    t : float
        Parameter for lower bound in nodes in the boundary
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0

    Returns
    -------
    dub : array
        Values of the sensitivites of upper bound in the points
    dlb : array
        Values of the sensitivites of lower bound in the points
    """

    xc = center[0]
    yc = center[1]
    ri = radius - thk/2
    re = radius + thk/2
    dub = zeros((len(x), 1))
    dlb = zeros((len(x), 1))
    dubdx = zeros((len(x), len(x)))
    dubdy = zeros((len(x), len(x)))
    dlbdx = zeros((len(x), len(x)))
    dlbdy = zeros((len(x), len(x)))

    for i in range(len(x)):
        zi2 = ri**2 - (x[i] - xc)**2 - (y[i] - yc)**2
        ze2 = re**2 - (x[i] - xc)**2 - (y[i] - yc)**2
        ze = math.sqrt(ze2)
        dub[i] = 1/2 * re/ze
        dubdx[i, i] = 1/2/ze * -2 * (x[i] - xc)
        dubdy[i, i] = 1/2/ze * -2 * (y[i] - yc)
        if zi2 > 0.0:
            zi = math.sqrt(zi2)
            dlb[i] = -1/2 * ri/zi
            dlbdx[i, i] = 1/2/zi * -2 * (x[i] - xc)
            dlbdy[i, i] = 1/2/zi * -2 * (y[i] - yc)

    return dub, dlb, dubdx, dubdy, dlbdx, dlbdy


def dome_b_update(x, y, thk, fixed, center=[5.0, 5.0], radius=5.0):
    """Updates the ``b`` parameter of a dome for a given thickness

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the dome
    fixed : list
        List with indexes of the fixed vertices
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0

    Returns
    -------
    b : array
        The ``b`` parameter
    """

    [xc, yc] = center[:2]
    b = zeros((len(fixed), 2))

    for i in range(len(fixed)):
        i_ = fixed[i]
        theta = math.atan2((y[i_] - yc), (x[i_] - xc))
        x_ = abs(thk/2 * math.cos(theta))
        y_ = abs(thk/2 * math.sin(theta))
        b[i, :] = [x_, y_]

    return b


def dome_db(x, y, thk, fixed, center=[5.0, 5.0], radius=5.0):
    """Updates the ``db`` parameter of a dome for a given thickness

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    thk : float
        Thickness of the dome
    fixed : list
        List with indexes of the fixed vertices
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]
    radius : float, optional
        The radius of the dome, by default 5.0

    Returns
    -------
    b : array
        The sensitivity of the ``b`` parameter
    """

    [xc, yc] = center[:2]
    db = zeros((len(fixed), 2))

    for i in range(len(fixed)):
        i_ = fixed[i]
        theta = math.atan2((y[i_] - yc), (x[i_] - xc))
        x_ = abs(1/2 * math.cos(theta))
        y_ = abs(1/2 * math.sin(theta))
        db[i, :] = [x_, y_]

    return db


def dome_b_update_with_n(x, y, n, fixed, b, center=[5.0, 5.0]):
    """Updates the ``b`` parameter of a dome for a given thickness considering as variable the offset from the surface

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    n : float
        Offset distance
    fixed : list
        List with indexes of the fixed vertices
    b : list
        Current ``b`` list
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]

    Returns
    -------
    new_b : array
        The new ``b`` parameter
    """

    [xc, yc] = center[:2]
    new_b = zeros((len(fixed), 2))

    for i in range(len(fixed)):
        i_ = fixed[i]
        theta = math.atan2((y[i_] - yc), (x[i_] - xc))
        new_b[i, 0] = b[i, 0] - n * abs(math.cos(theta))
        new_b[i, 1] = b[i, 1] - n * abs(math.sin(theta))

    return new_b


def dome_db_with_n(x, y, fixed, center=[5.0, 5.0]):
    """Updates the sensitivites of the ``b`` parameter of a dome for a given thickness considering as variable the offset from the surface

    Parameters
    ----------
    x : list
        x-coordinates of the points
    y : list
        y-coordinates of the points
    n : float
        Offset distance
    fixed : list
        List with indexes of the fixed vertices
    b : list
        Current ``b`` list
    center : [float, float], optional
        x, y coordinates of the center of the dome, by default [5.0, 5.0]

    Returns
    -------
    db : array
        The sensitivity of the ``b`` parameter
    """

    [xc, yc] = center[:2]
    db = zeros((len(fixed), 2))

    for i in range(len(fixed)):
        i_ = fixed[i]
        theta = math.atan2((y[i_] - yc), (x[i_] - xc))
        db[i, 0] = - abs(math.cos(theta))
        db[i, 1] = - abs(math.sin(theta))

    return db
