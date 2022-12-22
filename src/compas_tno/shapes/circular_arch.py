from compas_tno.shapes import MeshDos
from compas.datastructures import mesh_weld
from numpy import ones
from numpy import zeros
from numpy import array
import math


def arch_shape(H=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, t=0.0, total_nodes=100):
    """Helper to create meshes to define upper and lower bounds of 2D arch

    Parameters
    ----------
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0
    thk : float, optional
        Thickness of the arch, by default 0.20
    b : float, optional
        Out-of-plane measure of the arch, by default 0.5
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 0.0
    total_nodes : int, optional
        Density of the shape, by default 100

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the pattern
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the pattern
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the pattern

    """
    # Add option for starting from Hi and Li for a given thk.

    radius = H / 2 + (L**2 / (8 * H))
    ri = radius - thk/2
    re = radius + thk/2

    spr = math.atan2((L/2), (radius - H))
    tot_angle = 2*spr
    angle_init = (math.pi - tot_angle)/2
    an = tot_angle / (total_nodes - 1)
    zc = radius - H
    xc = L/2 + x0
    i = 0

    xs = []
    ys = []
    zts = []
    zis = []
    zes = []
    faces = []

    for i in range(total_nodes-1):
        angle_i = angle_init + i * an
        angle_f = angle_init + (i + 1) * an
        xi = xc - radius * math.cos(angle_i)
        xf = xc - radius * math.cos(angle_f)
        zi = radius * math.sin(angle_i) - zc
        zf = radius * math.sin(angle_f) - zc
        zei = math.sqrt(re**2 - (xi - xc)**2) - zc
        zef = math.sqrt(re**2 - (xf - xc)**2) - zc
        zii2 = ri**2 - (xi - xc)**2
        zif2 = ri**2 - (xf - xc)**2
        if zii2 > 0:
            zii = math.sqrt(ri**2 - (xi - xc)**2) - zc
        else:
            zii = -t
        if zif2 > 0:
            zif = math.sqrt(ri**2 - (xf - xc)**2) - zc
        else:
            zif = -t

        xs += [xi, xi, xi, xf, xf, xf]
        ys += [-b/2, 0, b/2, -b/2, 0, b/2]
        zts += [zi, zi, zi, zf, zf, zf]
        zis += [zii, zii, zii, zif, zif, zif]
        zes += [zei, zei, zei, zef, zef, zef]
        faces.append([6*i, 6*i + 1, 6*i + 4, 6*i + 3])
        faces.append([6*i + 1, 6*i + 2, 6*i + 5, 6*i + 4])

        i = i + 1

    intrados = mesh_weld(MeshDos.from_vertices_and_faces(array([xs, ys, zis]).transpose(), faces))
    extrados = mesh_weld(MeshDos.from_vertices_and_faces(array([xs, ys, zes]).transpose(), faces))
    middle = mesh_weld(MeshDos.from_vertices_and_faces(array([xs, ys, zts]).transpose(), faces))

    return intrados, extrados, middle


def arch_shape_polar(H=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, total_nodes=100):
    """Helper to create meshes to define upper and lower bounds of 2D arch with polar coordinates
    Note: The points in this mesh will not result in a direct projection for a typical form diagram

    Parameters
    ----------
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0
    thk : float, optional
        Thickness of the arch, by default 0.20
    b : float, optional
        Out-of-plane measure of the arch, by default 0.5
    total_nodes : int, optional
        Density of the shape, by default 100

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the pattern
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the pattern
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the pattern

    """
    # Add option for starting from Hi and Li for a given thk.

    radius = H / 2 + (L**2 / (8 * H))
    ri = radius - thk/2
    re = radius + thk/2

    spr = math.atan2((L/2), (radius - H))
    tot_angle = 2*spr
    angle_init = (math.pi - tot_angle)/2
    an = tot_angle / (total_nodes - 1)
    zc = radius - H
    xc = L/2 + x0
    i = 0

    xis = []
    xes = []
    xts = []
    ys = []
    zts = []
    zis = []
    zes = []
    faces = []

    for i in range(total_nodes-1):
        angle_i = angle_init + i * an
        angle_f = angle_init + (i + 1) * an

        xii = xc - ri * math.cos(angle_i)  # intrados initial
        xif = xc - ri * math.cos(angle_f)  # intrados final
        xei = xc - re * math.cos(angle_i)  # extrados initial
        xef = xc - re * math.cos(angle_f)  # extrados final
        xti = xc - radius * math.cos(angle_i)  # target initial
        xtf = xc - radius * math.cos(angle_f)  # target final

        zii = zc + ri * math.sin(angle_i)  # intrados initial
        zif = zc + ri * math.sin(angle_f)  # intrados final
        zei = zc + re * math.sin(angle_i)  # extrados initial
        zef = zc + re * math.sin(angle_f)  # extrados final
        zti = zc + radius * math.sin(angle_i)  # target initial
        ztf = zc + radius * math.sin(angle_f)  # target final

        xis += [xii, xii, xii, xif, xif, xif]
        xes += [xei, xei, xei, xef, xef, xef]
        xts += [xti, xti, xti, xtf, xtf, xtf]
        ys += [-b/2, 0, b/2, -b/2, 0, b/2]
        zts += [zti, zti, zti, ztf, ztf, ztf]
        zis += [zii, zii, zii, zif, zif, zif]
        zes += [zei, zei, zei, zef, zef, zef]
        faces.append([6*i, 6*i + 1, 6*i + 4, 6*i + 3])
        faces.append([6*i + 1, 6*i + 2, 6*i + 5, 6*i + 4])

        i = i + 1

    intrados = mesh_weld(MeshDos.from_vertices_and_faces(array([xis, ys, zis]).transpose(), faces))
    extrados = mesh_weld(MeshDos.from_vertices_and_faces(array([xes, ys, zes]).transpose(), faces))
    middle = mesh_weld(MeshDos.from_vertices_and_faces(array([xts, ys, zts]).transpose(), faces))

    return intrados, extrados, middle


def arch_ub_lb_update(x, y, thk, t, H=1.00, L=2.0, x0=0.0):
    """Update upper and lower bounds of an arch based in the parameters

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
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """

    radius = H / 2 + (L**2 / (8 * H))
    ri = radius - thk/2
    re = radius + thk/2
    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t
    zc = radius - H
    xc = L/2 + x0

    for i in range(len(x)):
        ub[i] = math.sqrt(re**2 - (x[i] - xc)**2) - zc
        lb2 = ri**2 - (x[i] - xc)**2
        if lb2 > 0:
            lb[i] = math.sqrt(lb2) - zc

    return ub, lb


def arch_dub_dlb(x, y, thk, t, H=1.00, L=2.0, x0=0.0):
    """Computes the sensitivities of upper and lower bounds in the x, y coordinates and thickness specified.

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
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0

    Returns
    -------
    dub : array
        Values of the sensitivities for the upper bound in the points
    dlb : array
        Values of the sensitivities for the lower bound in the points
    """

    radius = H / 2 + (L**2 / (8 * H))
    ri = radius - thk/2
    re = radius + thk/2
    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t
    dub = zeros((len(x), 1))
    dlb = zeros((len(x), 1))
    zc = radius - H
    xc = L/2 + x0

    for i in range(len(x)):
        ze = math.sqrt(re**2 - (x[i] - xc)**2) - zc
        ub[i] = ze
        dub[i] = re/(2*ze)
        zi2 = ri**2 - (x[i] - xc)**2
        if zi2 > 0:
            zi = math.sqrt(zi2) - zc
            lb[i] = zi
            dlb[i] = - ri/(2*zi)

    return dub, dlb  # ub, lb


def arch_b_update(x, y, thk, fixed, H=1.00, L=2.0, x0=0.0):
    """Computes the sensitivities of upper and lower bounds in the x, y coordinates and thickness specified.

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
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0

    Returns
    -------
    b : array
        Value of the ``b`` parameter.
    """

    radius = H / 2 + (L**2 / (8 * H))
    re = radius + thk/2
    zc = radius - H
    b = zeros((len(fixed), 2))
    x_lim = math.sqrt(re**2 - zc**2)

    for i in range(len(fixed)):
        b[i, :] = [x_lim - L/2, 0.0]

    return b


def arch_db(x, y, thk, fixed, H=1.00, L=2.0, x0=0.0):
    """Computes the sensitivities of upper and lower bounds in the x, y coordinates and thickness specified.

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
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.0
    x0 : float, optional
        Starting coordinate of the arch , by default 0.0

    Returns
    -------
    db : array
        Sensitivity of the ``b`` parameter.
    """

    radius = H / 2 + (L**2 / (8 * H))
    re = radius + thk/2
    zc = radius - H
    db = zeros((len(fixed), 2))
    x_lim = math.sqrt(re**2 - zc**2)

    for i in range(len(fixed)):
        db[i, :] = [1/2 * re/x_lim, 0.0]

    return db
