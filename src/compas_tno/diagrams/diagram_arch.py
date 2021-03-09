from compas.datastructures import Mesh
from compas.utilities import geometric_key
from compas_tno.datastructures import MeshDos
from numpy import linspace
import math


def create_arch(cls, H=1.00, L=2.0, x0=0.0, total_nodes=100):
    """ Helper to create a arch linear form-diagram with equaly spaced (in 3D) nodes.

    Parameters
    ----------
    H : float
        Rise of the arch measured with regards to the center line.

    L : float
        Span of the arch considered as center, to center. (L <= 2*H).

    x0: float
        Beginning of the linear form diagram.

    total_nodes : int
        Numbers of nodes to be considered in the form diagram.

    Returns
    -------
    obj
        FormDiagram.

    """

    # Add option for starting from Hi and Li for a given thk.

    radius = H / 2 + (L**2 / (8 * H))
    print('radius =', radius)
    spr = math.atan2((L/2), (radius - H))
    print('springing angle =', math.degrees(spr))
    tot_angle = 2*spr
    angle_init = (math.pi - tot_angle)/2
    an = tot_angle / (total_nodes - 1)
    lines = []
    gkey_fix = []

    for i in range(total_nodes-1):
        angle_i = angle_init + i * an
        angle_f = angle_init + (i + 1) * an
        xi = L/2 - radius * math.cos(angle_i)
        xf = L/2 - radius * math.cos(angle_f)
        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi, 0.0, 0.0], precision=6))
        elif i == total_nodes - 2:
            gkey_fix.append(geometric_key([xf, 0.0, 0.0], precision=6))

    mesh = Mesh.from_lines(lines)
    form = cls.from_mesh(mesh)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form


def create_linear_form_diagram(cls, L=2.0, x0=0.0, total_nodes=100):
    """ Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes.

    Parameters
    ----------
    L : float
        Span of the arch considered as center, to center.
    x0: float
        Beginning of the linear form diagram.
    total_nodes : int
        Numbers of nodes to be considered in the form diagram.

    Returns
    -------
    obj
        FormDiagram.

    """

    x = linspace(x0, x0 + L, num=total_nodes, endpoint=True)  # Continue this
    lines = []
    gkey_fix = []

    for i in range(total_nodes-1):
        xi = x[i]
        xf = x[i + 1]
        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi, 0.0, 0.0], precision=6))
        elif i == total_nodes - 2:
            gkey_fix.append(geometric_key([xf, 0.0, 0.0], precision=6))

    mesh = MeshDos.from_lines(lines)
    form = cls.from_mesh(mesh)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form
