from compas.datastructures import Mesh
from compas.utilities import geometric_key
import math


def create_arch_form_diagram(cls, H=1.0, L=2.0, x0=0.0, discretisation=100):
    """Construct a FormDiagram based on an arch linear discretisation.
    Note: The nodes of the form diagram are spaced following a projection in a semicircular arch.

    Parameters
    ----------
    H : float, optional
        Height of the arch, by default 1.00
    L : float, optional
        Span of the arch, by default 2.00
    x0 : float, optional
        Initial coordiante of the arch, by default 0.0
    discretisation : int, optional
        Numbers of nodes to be considered in the form diagram, by default 100

    Returns
    -------
    :class:`~compas_tno.diagrams.FormDiagram`
        The FormDiagram created.
    """

    # Add option for starting from Hi and Li for a given thk.

    radius = H / 2 + (L**2 / (8 * H))
    print('radius =', radius)
    spr = math.atan2((L/2), (radius - H))
    print('springing angle =', math.degrees(spr))
    tot_angle = 2*spr
    angle_init = (math.pi - tot_angle)/2
    an = tot_angle / (discretisation - 1)
    lines = []
    gkey_fix = []

    for i in range(discretisation-1):
        angle_i = angle_init + i * an
        angle_f = angle_init + (i + 1) * an
        xi = L/2 - radius * math.cos(angle_i) + x0
        xf = L/2 - radius * math.cos(angle_f) + x0
        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi, 0.0, 0.0], precision=6))
        elif i == discretisation - 2:
            gkey_fix.append(geometric_key([xf, 0.0, 0.0], precision=6))

    mesh = Mesh.from_lines(lines)
    form = cls.from_mesh(mesh)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form


def create_linear_form_diagram(cls, L=2.0, x0=0.0, discretisation=100):
    """ Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes.

    Parameters
    ----------
    L : float, optional
        Span of the arch, by default 2.00
    x0 : float, optional
        Initial coordiante of the arch, by default 0.0
    discretisation : int, optional
        Numbers of nodes to be considered in the form diagram, by default 100

    Returns
    -------
    :class:`~compas_tno.diagrams.FormDiagram`
        FormDiagram generated according to the parameters.

    """

    from numpy import linspace
    x = linspace(x0, x0 + L, discretisation)  # Continue this remove need of numpy in the future
    lines = []
    gkey_fix = []

    for i in range(discretisation-1):
        xi = x[i]
        xf = x[i + 1]
        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])
        if i == 0:
            gkey_fix.append(geometric_key([xi, 0.0, 0.0], precision=6))
        elif i == discretisation - 2:
            gkey_fix.append(geometric_key([xf, 0.0, 0.0], precision=6))

    mesh = Mesh.from_lines(lines)
    form = cls.from_mesh(mesh)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form


def create_linear_form_diagram_sp_ep(cls, sp=[0, 0, 0], ep=[2, 0, 0], discretisation=100):
    """ Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes based on starting and ending points

    Parameters
    ----------
    L : float, optional
        Span of the arch, by default 2.00
    x0 : float, optional
        Initial coordiante of the arch, by default 0.0
    discretisation : int, optional
        Numbers of nodes to be considered in the form diagram, by default 100

    Returns
    -------
    :class:`~compas_tno.diagrams.FormDiagram`
        FormDiagram generated according to the parameters.

    """

    x0, y0 = sp[:2]
    xf, yf = ep[:2]

    dx = (xf - x0)/discretisation
    dy = (yf - y0)/discretisation

    lines = []
    gkey_fix = []

    for i in range(discretisation):
        xi = x0 + i * dx
        xf = x0 + (i + 1) * dx
        yi = y0 + i * dy
        yf = y0 + (i + 1) * dy

        lines.append([[xi, yi, 0.0], [xf, yf, 0.0]])

        if i == 0:
            gkey_fix.append(geometric_key([xi, yi, 0.0], precision=6))
        elif i == discretisation - 1:
            gkey_fix.append(geometric_key([xf, yf, 0.0], precision=6))

    mesh = Mesh.from_lines(lines)
    form = cls.from_mesh(mesh)
    gkey_key = form.gkey_key(precision=6)

    form.vertex_attribute(gkey_key[gkey_fix[0]], 'is_fixed', True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], 'is_fixed', True)

    return form
