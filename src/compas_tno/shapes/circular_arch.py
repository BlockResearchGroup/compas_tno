from compas.datastructures import Mesh
from compas.utilities import geometric_key
from numpy import array
import math


def arch_shape(H=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, t=5.0, total_nodes=100):
    """ Helper to create meshes to define upper and lower bounds of 2D arch.

    Parameters
    ----------
    H : float
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

    radius = radius = H / 2 + (L**2 / (8 * H))
    ri = radius - thk/2
    re = radius + thk/2
    print('radius/ ri / re =', radius, ri, re)
    spr = math.atan2((L/2),(radius - H))
    print('springing angle =', math.degrees(spr), spr)
    tot_angle = 2*spr
    angle_init = (math.pi - tot_angle)/2
    print('init angle =', math.degrees(angle_init))
    an = tot_angle / (total_nodes - 1)
    zc = radius - H
    xc = L/2
    lines = []
    gkey_fix = []
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
        zei = (re**2 - (xi - xc)**2)**(1/2) - zc
        zef = (re**2 - (xf - xc)**2)**(1/2) - zc
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

        if i == 0:
            gkey_fix.append(geometric_key([xi, 0.0, 0.0], precision=6))
        elif i == total_nodes - 2:
            gkey_fix.append(geometric_key([xf, 0.0, 0.0], precision=6))

        i = i + 1

    intrados = Mesh.from_vertices_and_faces(array([xs, ys, zis]).transpose(), faces)
    extrados = Mesh.from_vertices_and_faces(array([xs, ys, zes]).transpose(), faces)
    middle = Mesh.from_vertices_and_faces(array([xs, ys, zts]).transpose(), faces)

    # form = FormDiag.from_lines(lines, delete_boundary_face=False)
    # gkey_key = form.gkey_key(precision=6)
    # form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    # form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    # intrados = None
    # extrados = None
    # middle = None

    return intrados, extrados, middle
