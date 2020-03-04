from compas_viewers.meshviewer import MeshViewer
from compas_viewers.multimeshviewer import MultiMeshViewer
from compas.datastructures import Mesh

__all__ = [
    'view_thrust',
    'view_solution',
]

def view_thrust(form):
    """ Viewer showing the thrust network.

    Parameters
    ----------
    form : FormDiagram
        FormDiagram to plot

    Returns
    ----------
    obj
        Plotter object.

    """

    mesh = form
    viewer = MeshViewer()
    viewer.mesh = mesh

    return viewer

def view_solution(form, shape):
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
        Plotter object.

    """

    mesh = form
    intrados = shape.intrados
    extrados = shape.extrados

    viewer = MultiMeshViewer()
    viewer.meshes = [mesh, intrados, extrados]

    return viewer
