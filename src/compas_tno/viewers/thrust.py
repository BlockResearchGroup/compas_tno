from compas_viewers.objectviewer import ObjectViewer
from compas.geometry import Line
from compas_tno.shapes import Shape

from compas.geometry import Point

__all__ = [
    'view_thrust',
    'view_thrusts',
    'view_solution',
    'view_thrust_as_lines',
    'view_bestfit_solution'
]


def view_thrust(form, settings_form=None, cracks=True):
    """ Viewer showing the thrust network.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    settings_form : dict (None)
        Settings to modify the plot
    cracks : bool (True)
        Option to show the cracks as points

    Returns
    ----------
    obj
        Viewer object.

    """

    viewer = ObjectViewer()

    intrad = 1
    extrad = 1

    for key in form.vertices():
        lb = form.vertex_attribute(key, 'lb')
        ub = form.vertex_attribute(key, 'ub')
        x, y, z = form.vertex_coordinates(key)
        if ub is not None and lb is not None:
            if cracks:
                if abs(ub - z) < 10e-4:
                    viewer.add(Point(x, y, z), name="Extrados (%s)" % extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
                    extrad += 1
                elif abs(lb - z) < 10e-4:
                    viewer.add(Point(x, y, z), name="Intrados (%s)" % intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
                    intrad += 1

    if not settings_form:
        settings_form = {
            'color': '#FF0000',
            'edges.color': '#FF0000',
            'edges.width': 2,
            'opacity': 0.8,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': False,
        }

    viewer.add(form, name="FormDiagram", settings=settings_form)

    return viewer


def view_thrust_as_lines(form, settings_form=None, cracks=True, viewer=None):
    """ Viewer showing the thrust network.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    settings_form : dict (None)
        Settings to modify the plot
    cracks : bool (True)
        Option to show the cracks as points

    Returns
    ----------
    obj
        Viewer object.

    """

    if not viewer:
        viewer = ObjectViewer()

    f = [abs(form.edge_attribute((u, v), 'q') * form.edge_length(u, v)) for u, v in form.edges_where({'_is_edge': True})]
    fmax = max(f)
    max_width = 0.5
    tol = 10e-4

    i = 0
    for u, v in form.edges_where({'_is_edge': True}):
        fi = f[i]
        if fi > tol:
            width = (fi / fmax) * max_width
            start, end = form.edge_coordinates(u, v)
            line = Line(start, end)
            viewer.add(line, name="i - f ({0} - {1:.2f})".format(i, fi), settings={'color': '#FF0000', 'width': width})

        i += 1

    # add reactions

    return viewer


def view_thrusts(forms, settings_form=None, cracks=True):
    """ Viewer showing the thrust network.

    Parameters
    ----------
    form : list of FormDiagram
        List og FormDiagrams to plot
    settings_form : dict (None)
        Settings to modify the plot
    cracks : bool (True)
        Option to show the cracks as points


    Returns
    ----------
    obj
        Viewer object.

    """

    viewer = ObjectViewer()

    for form in forms:

        intrad = 1
        extrad = 1

        for key in form.vertices():
            lb = form.vertex_attribute(key, 'lb')
            ub = form.vertex_attribute(key, 'ub')
            x, y, z = form.vertex_coordinates(key)
            if cracks:
                if abs(ub - z) < 10e-4:
                    viewer.add(Point(x, y, z), name="Extrados (%s)" % extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
                    extrad += 1
                elif abs(lb - z) < 10e-4:
                    viewer.add(Point(x, y, z), name="Intrados (%s)" % intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
                    intrad += 1

        if not settings_form:
            settings_form = {
                'color': '#FF0000',
                'edges.color': '#FF0000',
                'edges.width': 2,
                'opacity': 0.8,
                'vertices.size': 0,
                'vertices.on': False,
                'edges.on': True,
                'faces.on': False,
            }

        viewer.add(form, name="FormDiagram", settings=settings_form)

    return viewer


def view_solution(form, shape=None, settings_form=None, settings_bounds=None, cracks=True, reactions=True, thickness=False, outside=True):
    """ Viewer showing the thrust network together with intrados and extrados.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    shape : Shape, optional
        Shape to plot
        If no Shape is given, the shape is constructed from the form's attributes ``UB`` and ``LB``.

    Returns
    ----------
    obj
        Viewer object.

    """

    if not shape:
        shape = Shape.from_formdiagram_and_attributes(form)

    intrados = shape.intrados
    extrados = shape.extrados

    viewer = ObjectViewer()

    intrad = 1
    extrad = 1
    out = 1

    for key in form.vertices():
        lb = form.vertex_attribute(key, 'lb')
        ub = form.vertex_attribute(key, 'ub')
        x, y, z = form.vertex_coordinates(key)
        if cracks:
            if abs(ub - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Extrados (%s)" % extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
                extrad += 1
            elif abs(lb - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Intrados (%s)" % intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
                intrad += 1
        if outside:
            if z > ub:
                viewer.add(Point(x, y, z), name="Outside - Intra (%s)" % out, settings={'vertices.color': '#000000', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
                out += 1
            elif z < lb:
                viewer.add(Point(x, y, z), name="Outside - Extra (%s)" % out, settings={'vertices.color': '#000000', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
                out += 1

    if not settings_form:
        settings_form = {
            'color': '#FF0000',
            'edges.color': '#FF0000',
            'edges.width': 2,
            'opacity': 0.8,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': False,
        }

    if not settings_bounds:
        settings_bounds = {
            'color': '#999999',
            'edges.width': 3,
            'opacity': 0.5,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': False,
            'faces.on': True,
        }

    viewer.add(intrados, name="Intrados", settings=settings_bounds)
    viewer.add(extrados, name="Extrados", settings=settings_bounds)

    if not thickness:
        viewer.add(form, name="FormDiagram", settings=settings_form)
    else:
        # f_range = [form.edge_attribute((u, v), 'f') for u, v in form.edges()]
        for u, v in form.edges_where({'_is_edge': True}):
            return
            # Figure out how to make this plot.

    return viewer


def view_bestfit_solution(form, shape=None, settings_form=None, settings_bounds=None, reactions=True, thickness=False, outside=True):
    """ Viewer showing the thrust network together with intrados and extrados.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    shape : Shape, optional
        Shape to plot
        If no Shape is given, the shape is constructed from the form's attributes ``UB`` and ``LB``.

    Returns
    ----------
    obj
        Viewer object.

    """

    if not shape:
        shape = Shape.from_formdiagram_and_attributes(form)

    middle = shape.middle

    viewer = ObjectViewer()

    # for key in form.vertices():
    #     lb = form.vertex_attribute(key, 'lb')
    #     ub = form.vertex_attribute(key, 'ub')
    #     x, y, z = form.vertex_coordinates(key)

    if not settings_form:
        settings_form = {
            'color': '#FF0000',
            'edges.color': '#FF0000',
            'edges.width': 2,
            'opacity': 0.8,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': False,
        }

    if not settings_bounds:
        settings_bounds = {
            'color': '#999999',
            'edges.width': 3,
            'opacity': 0.5,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': False,
            'faces.on': True,
        }

    viewer.add(middle, name="Target", settings=settings_bounds)
    viewer.add(form, name="FormDiagram", settings=settings_form)

    return viewer
