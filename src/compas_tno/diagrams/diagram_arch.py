
from compas.utilities import geometric_key
import math


def create_arch(FormDiag, D=2.00, x0=0.0, total_nodes=100):
    """ Helper to create a arch linear form-diagram.

    Parameters
    ----------
    D : float
        Central diameter of the arch.

    x0: float
        Beginning of the linear form diagram.

    total_nodes : int
        Numbers of nodes to be considered in the form diagram.

    Returns
    -------
    obj
        FormDiagram.

    """
    r = D/2
    xc = x0 + r
    lines = []
    total_edges = total_nodes - 1
    gkey_fix = []

    for i in range(total_edges):
        xi = xc - r*math.cos(i/total_edges*math.pi)
        xf = xc - r*math.cos((i+1)/total_edges*math.pi)
        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi, 0.0, 0.0], precision=6))
        elif i == total_edges - 1:
            gkey_fix.append(geometric_key([xf, 0.0, 0.0], precision=6))

    form = FormDiag.from_lines(lines, delete_boundary_face=False)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form
