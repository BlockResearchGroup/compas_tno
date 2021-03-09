from compas_tno.datastructures import MeshDos
from compas.utilities import geometric_key
from numpy import ones
from numpy import zeros
from numpy import array
import math
from compas_tno.shapes.crossvault import rectangular_topology
from numpy import linspace


def pointed_arch_shape(hc=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, t=5.0, total_nodes=100):
    """ Helper to create meshes to define upper and lower bounds of 2D pointed arch.

    Parameters
    ----------
    hc : float
        Rise of the arch measured with regards to the center line.
    L : float
        Span of the arch considered as center, to center. (L <= 2*H).
    x0: float
        Beginning of the linear form diagram.
    thk: float
        Thickness of the arch
    b: float
        Out of plane dimension of the arch

    total_nodes : int
        Numbers of nodes to be considered in the form diagram.

    Returns
    -------
    obj
        FormDiagram.

    """

    # Add option for starting from Hi and Li for a given thk.

    radius = 1/L * (hc**2 + L**2/4)
    ri = radius - thk/2
    re = radius + thk/2
    print('radius/ ri / re =', radius, ri, re)
    zc = 0.0
    xc1 = x0 + radius
    xc2 = x0 + L - radius
    print('Centers x =', xc1, xc2)

    x = linspace(x0, x0 + L, num=total_nodes, endpoint=True)
    xs, ys, faces_i = rectangular_topology(x, [-b/2, 0, b/2])

    i = 0

    zts = []
    zis = []
    zes = []

    for xi in x:
        if xi <= L/2:
            dx = xi - xc1
        else:
            dx = xc2 - xi

        zi = math.sqrt(radius**2 - dx**2) + zc
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

    radius = 1/L * (hc**2 + L**2/4)
    ri = radius - thk/2
    re = radius + thk/2
    zc = 0.0
    xc1 = x0 + radius
    xc2 = x0 + L - radius

    ub = ones((len(x), 1))
    lb = ones((len(x), 1)) * - t

    for i in range(len(x)):
        if x[i] <= L/2:
            dx = x[i] - xc1
        else:
            dx = xc2 - x[i]
        ub[i] = math.sqrt(re**2 - dx**2) - zc
        lb2 = ri**2 - dx**2
        if lb2 > 0:
            lb[i] = math.sqrt(lb2) - zc

    return ub, lb


def pointed_arch_dub_dlb(x, y, thk, t, hc, L, x0):

    radius = 1/L * (hc**2 + L**2/4)
    ri = radius - thk/2
    re = radius + thk/2
    xc1 = x0 + radius
    xc2 = x0 + L - radius
    zc = 0.0

    dub = zeros((len(x), 1))
    dlb = zeros((len(x), 1))

    for i in range(len(x)):
        if x[i] <= L/2:
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


def pointed_arch_b_update(x, y, thk, fixed, H=1.00, L=2.0, x0=0.0):

    radius = H / 2 + (L**2 / (8 * H))
    re = radius + thk/2
    zc = radius - H
    b = zeros((len(fixed), 2))
    x_lim = math.sqrt(re**2 - zc**2)

    for i in range(len(fixed)):
        b[i, :] = [x_lim - L/2, 0.0]

    return b


def pointed_arch_db(x, y, thk, fixed, H=1.00, L=2.0, x0=0.0):  # This does not work for spring angles

    radius = H / 2 + (L**2 / (8 * H))
    re = radius + thk/2
    zc = radius - H
    db = zeros((len(fixed), 2))
    x_lim = math.sqrt(re**2 - zc**2)

    for i in range(len(fixed)):
        db[i, :] = [1/2 * re/x_lim, 0.0]

    return db
