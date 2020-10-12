from compas_viewers.meshviewer import MeshViewer
from compas_viewers.multimeshviewer import MultiMeshViewer
from compas_viewers.objectviewer import ObjectViewer
from compas.datastructures import Mesh

from compas.geometry import Point

__all__ = [
    'view_thrust',
    'view_thrusts',
    'view_solution',
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
                    viewer.add(Point(x, y, z), name="Extrados (%s)"%extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
                    extrad += 1
                elif abs(lb - z) < 10e-4:
                    viewer.add(Point(x, y, z), name="Intrados (%s)"%intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
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
                    viewer.add(Point(x, y, z), name="Extrados (%s)"%extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
                    extrad += 1
                elif abs(lb - z) < 10e-4:
                    viewer.add(Point(x, y, z), name="Intrados (%s)"%intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
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

def view_solution(form, shape, settings_form=None, settings_bounds=None, cracks=True):
    """ Viewer showing the thrust network together with intrados and extrados.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot
    shape : Shape
        Shape to plot

    Returns
    ----------
    obj
        Viewer object.

    """

    intrados = shape.intrados
    extrados = shape.extrados

    viewer = ObjectViewer()

    intrad = 1
    extrad = 1

    for key in form.vertices():
        lb = form.vertex_attribute(key, 'lb')
        ub = form.vertex_attribute(key, 'ub')
        x, y, z = form.vertex_coordinates(key)
        if cracks:
            if abs(ub - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Extrados (%s)"%extrad, settings={'vertices.color': '#008000', 'vertices.size': 15.0, 'vertices.on': True})  # green extrados
                extrad += 1
            elif abs(lb - z) < 10e-4:
                viewer.add(Point(x, y, z), name="Intrados (%s)"%intrad, settings={'vertices.color': '#0000FF', 'vertices.size': 15.0, 'vertices.on': True})  # blue intrados
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

    viewer.add(form, name="FormDiagram", settings=settings_form)
    viewer.add(intrados, name="Intrados", settings=settings_bounds)
    viewer.add(extrados, name="Extrados", settings=settings_bounds)

    return viewer
