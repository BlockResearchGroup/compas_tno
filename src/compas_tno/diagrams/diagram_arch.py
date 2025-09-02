import math
from typing import TYPE_CHECKING
from typing import Annotated
from typing import Type
from typing import Union

from numpy import linspace

from compas.geometry import Point
from compas.tolerance import TOL

if TYPE_CHECKING:
    from compas_tna.diagrams import FormDiagram


def create_arch_form_diagram(
    cls: Type["FormDiagram"],
    H: float = 1.0,
    L: float = 2.0,
    x0: float = 0.0,
    discretisation: int = 100,
) -> "FormDiagram":
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
    spr = math.atan2((L / 2), (radius - H))
    tot_angle = 2 * spr
    angle_init = (math.pi - tot_angle) / 2
    an = tot_angle / (discretisation - 1)

    lines = []
    gkey_fix = []

    for i in range(discretisation - 1):
        angle_i = angle_init + i * an
        angle_f = angle_init + (i + 1) * an
        xi = L / 2 - radius * math.cos(angle_i) + x0
        xf = L / 2 - radius * math.cos(angle_f) + x0

        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])

        if i == 0:
            gkey_fix.append(TOL.geometric_key([xi, 0.0, 0.0], precision=6))

        elif i == discretisation - 2:
            gkey_fix.append(TOL.geometric_key([xf, 0.0, 0.0], precision=6))

    form = cls.from_lines(lines)
    gkey_key = {TOL.geometric_key(form.vertex_coordinates(vertex), precision=6): vertex for vertex in form.vertices()}

    form.vertex_attribute(gkey_key[gkey_fix[0]], "is_support", True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], "is_support", True)

    return form


def create_linear_form_diagram(
    cls: Type["FormDiagram"],
    L: float = 2.0,
    x0: float = 0.0,
    discretisation: int = 100,
) -> "FormDiagram":
    """Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes.

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
    x = linspace(x0, x0 + L, discretisation)  # Continue this remove need of numpy in the future
    lines = []
    gkey_fix = []

    for i in range(discretisation - 1):
        xi = x[i]
        xf = x[i + 1]

        lines.append([[xi, 0.0, 0.0], [xf, 0.0, 0.0]])

        if i == 0:
            gkey_fix.append(TOL.geometric_key([xi, 0.0, 0.0], precision=6))

        elif i == discretisation - 2:
            gkey_fix.append(TOL.geometric_key([xf, 0.0, 0.0], precision=6))

    form = cls.from_lines(lines)
    gkey_key = {TOL.geometric_key(form.vertex_coordinates(vertex), precision=6): vertex for vertex in form.vertices()}

    form.vertex_attribute(gkey_key[gkey_fix[0]], "is_support", True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], "is_support", True)

    return form


def create_linear_form_diagram_sp_ep(
    cls: Type["FormDiagram"],
    sp: Union[Point, Annotated[list[float], 3]] = [0, 0, 0],
    ep: Union[Point, Annotated[list[float], 3]] = [2, 0, 0],
    discretisation: int = 100,
) -> "FormDiagram":
    """Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes based on starting and ending points

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

    dx = (xf - x0) / discretisation
    dy = (yf - y0) / discretisation

    lines = []
    gkey_fix = []

    for i in range(discretisation):
        xi = x0 + i * dx
        xf = x0 + (i + 1) * dx
        yi = y0 + i * dy
        yf = y0 + (i + 1) * dy

        lines.append([[xi, yi, 0.0], [xf, yf, 0.0]])

        if i == 0:
            gkey_fix.append(TOL.geometric_key([xi, yi, 0.0], precision=6))

        elif i == discretisation - 1:
            gkey_fix.append(TOL.geometric_key([xf, yi, 0.0], precision=6))

    form = cls.from_lines(lines)
    gkey_key = {TOL.geometric_key(form.vertex_coordinates(vertex), precision=6): vertex for vertex in form.vertices()}

    form.vertex_attribute(gkey_key[gkey_fix[0]], "is_support", True)
    form.vertex_attribute(gkey_key[gkey_fix[1]], "is_support", True)

    return form
