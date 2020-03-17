from compas_viewers.meshviewer import MeshViewer
from compas_viewers.multimeshviewer import MultiMeshViewer
from compas.datastructures import Mesh


__all__ = [
    'view_intrados',
    'view_extrados',
    'view_middle',
    'view_shapes'
]


def view_intrados(shape):
    """ Viewer showing the intrados of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Plotter object.

    """

    mesh = shape.intrados
    viewer = MeshViewer()
    viewer.mesh = mesh

    return viewer

def view_extrados(shape):
    """ Viewer showing the extrados of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Plotter object.

    """

    mesh = shape.extrados
    viewer = MeshViewer()
    viewer.mesh = mesh


    return viewer

def view_middle(shape):
    """ Viewer showing the middle (target) surface of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Plotter object.

    """

    mesh = shape.middle
    viewer = MeshViewer()
    viewer.mesh = mesh

    return viewer

def view_shapes(shape, cut_negatives=True):
    """ Viewer showing the middle (target) surface of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot
    cut_negatives : boll
        If true, it will hid the negative allowance on the close-to-support points (WIP).

    Returns
    ----------
    obj
        Plotter object.

    """

    middle = shape.middle
    intrados = shape.intrados
    extrados = shape.extrados

    viewer = MultiMeshViewer()
    viewer.meshes = [middle, intrados, extrados]

    return viewer
