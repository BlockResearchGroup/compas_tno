from numpy import ones
from numpy import zeros
from numpy import array
from numpy import linspace

from compas_tno.shapes import rectangular_topology
from compas_tno.shapes import MeshDos

import math


def pointed_arch_shape(hc=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, t=5.0, total_nodes=100):
    """Function to create meshes to define upper and lower bounds of 2D pointed arch.

    Parameters
    ----------
    hc : float, optional
        Rise of the arch measured with regards to the center line, by default 1.00
    L : float, optional
        Span of the arch considered as center, to center. (L <= 2*H), by default 2.0
    x0 : float, optional
        Beginning of the linear form diagram, by default 0.0
    thk : float, optional
        Thickness of the arch, by default 0.20
    b : float, optional
        Out of plane dimension of the arch, by default 0.5
    t : float, optional
        Parameter for lower bound in nodes in the boundary, by default 5.0
    total_nodes : int, optional
        Numbers of nodes to be considered in the form diagram, by default 100

    Returns
    -------
    intrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the intrados of the shape
    extrados : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the extrados of the shape
    middle : :class:`~compas_tno.shapes.MeshDos`
        A MeshDos for the middle of the shape
    """

    radius = 1/L * (hc**2 + L**2/4)
    ri = radius - thk/2
    re = radius + thk/2
    # print('radius/ ri / re =', radius, ri, re)
    zc = 0.0
    xc1 = x0 + radius
    xc2 = x0 + L - radius
    # print('Centers x =', xc1, xc2)

    x = linspace(x0, x0 + L, num=total_nodes, endpoint=True)
    xs, ys, faces_i = rectangular_topology(x, [-b/2, 0, b/2])

    i = 0

    zts = []
    zis = []
    zes = []

    for xi in x:
        if xi <= x0 + L/2:
            dx = xi - xc1
        else:
            dx = xc2 - xi

        zi2 = radius**2 - dx**2
        if zi2 > 0:
            zi = math.sqrt(zi2) - zc
        else:
            zi = 0
        zei = math.sqrt(re**2 - dx**2) + zc
        zii2 = ri**2 - dx**2
        if zii2 > 0:
            zii = math.sqrt(zii2) - zc
        else:
            zii = -t

        zts += [zi, zi, zi]
        zis += [zii, zii, zii]
        zes += [zei, zei, zei]

        i = i + 1

    xyzlb = array([xs, ys, zis]).transpose()
    xyzub = array([xs, ys, zes]).transpose()
    xyzt = array([xs, ys, zts]).transpose()

    intrados = MeshDos.from_vertices_and_faces(xyzlb, faces_i)
    extrados = MeshDos.from_vertices_and_faces(xyzub, faces_i)
    middle = MeshDos.from_vertices_and_faces(xyzt, faces_i)

    return intrados, extrados, middle


def pointed_arch_ub_lb_update(x, y, thk, t, hc, L, x0):
    """Update upper and lower bounds of a pointed arch based in the parameters

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
    hc : float
        Rise of the arch measured with regards to the center line
    L : float
        Span of the arch considered as center, to center. (L <= 2*H)
    x0 : float
        Beginning of the linear form diagram

    Returns
    -------
    ub : array
        Values of the upper bound in the points
    lb : array
        Values of the lower bound in the points
    """

    radius = 1/L * (hc**2 + L**2/4)
    ri = radius - thk/2
    re = radius + thk/2
    zc = 0.0
    xc1 = x0 + radius
    xc2 = x0 + L - radius

    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t

    for i in range(len(x)):
        if x[i] <= x0 + L/2:
            dx = x[i] - xc1
        else:
            dx = xc2 - x[i]
        ub[i] = math.sqrt(re**2 - dx**2) - zc
        lb2 = ri**2 - dx**2
        if lb2 > 0:
            lb[i] = math.sqrt(lb2) - zc

    return ub, lb


def pointed_arch_dub_dlb(x, y, thk, t, hc, L, x0):
    """Computes the sensitivities of upper and lower bounds of a pointed arch based in the parameters

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
    hc : float
        Rise of the arch measured with regards to the center line
    L : float
        Span of the arch considered as center, to center. (L <= 2*H)
    x0 : float
        Beginning of the linear form diagram

    Returns
    -------
    dub : array
        Values of the sensitivities for the upper bound in the points
    dlb : array
        Values of the sensitivities for the lower bound in the points
    """

    radius = 1/L * (hc**2 + L**2/4)
    ri = radius - thk/2
    re = radius + thk/2
    xc1 = x0 + radius
    xc2 = x0 + L - radius
    zc = 0.0

    dub = zeros((len(x), 1))
    dlb = zeros((len(x), 1))

    for i in range(len(x)):
        if x[i] <= x0 + L/2:
            dx = x[i] - xc1
        else:
            dx = xc2 - x[i]
        ze = math.sqrt(re**2 - dx**2) - zc
        dub[i] = re/(2*ze)
        zi2 = ri**2 - dx**2
        if zi2 > 0:
            zi = math.sqrt(zi2) - zc
            dlb[i] = - ri/(2*zi)

    return dub, dlb  # ub, lb


# def pointed_arch_b_update(x, y, thk, fixed, H=1.00, L=2.0, x0=0.0):

#     radius = H / 2 + (L**2 / (8 * H))
#     re = radius + thk/2
#     zc = radius - H
#     b = zeros((len(fixed), 2))
#     x_lim = math.sqrt(re**2 - zc**2)

#     for i in range(len(fixed)):
#         b[i, :] = [x_lim - L/2, 0.0]

#     return b


# def pointed_arch_db(x, y, thk, fixed, H=1.00, L=2.0, x0=0.0):  # This does not work for spring angles

#     radius = H / 2 + (L**2 / (8 * H))
#     re = radius + thk/2
#     zc = radius - H
#     db = zeros((len(fixed), 2))
#     x_lim = math.sqrt(re**2 - zc**2)

#     for i in range(len(fixed)):
#         db[i, :] = [1/2 * re/x_lim, 0.0]

#     return db
