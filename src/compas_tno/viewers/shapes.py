from compas_viewers.multimeshviewer import MultiMeshViewer
from compas_viewers.objectviewer import ObjectViewer
from compas.geometry import Line


__all__ = [
    'view_mesh',
    'view_meshes',
    'view_intrados',
    'view_extrados',
    'view_middle',
    'view_shapes_pointcloud',
    'view_shapes',
    'view_fill',
    'view_normals'
]


def view_mesh(mesh, name='mesh', normals=False, length=1.0, settings=None):
    """ Viewer showing the extrados of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Viewer object.

    """

    viewer = ObjectViewer()

    if not settings:
        settings = {
            'color': '#999999',
            'edges.width': 1,
            'opacity': 0.5,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': True,
        }

    viewer.add(mesh, name=name, settings=settings)

    if normals:
        for key in mesh.vertices():
            n = mesh.vertex_attribute(key, 'n')
            pt0 = mesh.vertex_coordinates(key)
            pt1 = [pt0[0] + n[0]*length, pt0[1] + n[1]*length, pt0[2] + n[2]*length]
            line = Line(pt0, pt1)
            viewer.add(line, name='normal-{}'.format(key))

    return viewer


def view_meshes(meshes, settings=None):

    viewer = ObjectViewer()

    if not settings:
        settings = {
            'color': '#999999',
            'edges.width': 1,
            'opacity': 0.5,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': True,
        }

    for mesh in meshes:
        viewer.add(mesh, settings=settings)

    return viewer


def view_intrados(shape, settings_bounds=None):
    """ Viewer showing the intrados of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Viewer object.

    """

    return view_mesh(shape.intrados, name="Intrados", settings=settings_bounds)


def view_extrados(shape, settings_bounds=None):
    """ Viewer showing the extrados of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Viewer object.

    """

    return view_mesh(shape.extrados, name="Extrados", settings=settings_bounds)


def view_middle(shape, settings_middle=None):
    """ Viewer showing the middle (target) surface of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Viewer object.

    """

    return view_mesh(shape.middle, name="Middle/Target", settings=settings_middle)


def view_shapes_pointcloud(shape, settings_bounds=None):
    """ Viewer showing the extrados of a shape

    Parameters
    ----------
    shape : obj
        Shape to plot

    Returns
    ----------
    obj
        Viewer object.

    """

    intrados = shape.intrados
    extrados = shape.extrados

    viewer = ObjectViewer()

    if not settings_bounds:
        settings_bounds = {
            'color': '#999999',
            'edges.width': 1,
            'opacity': 0.5,
            'vertices.size': 10,
            'vertices.on': True,
            'edges.on': False,
            'faces.on': True,
        }

    viewer.add(intrados, name="Intrados", settings=settings_bounds)
    viewer.add(extrados, name="Extrados", settings=settings_bounds)

    return viewer


def view_shapes(shape, show_middle=True, settings_middle=None, settings_bounds=None):
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
        Viewer object.

    """

    intrados = shape.intrados
    extrados = shape.extrados
    middle = shape.middle

    viewer = ObjectViewer()

    if not settings_middle:
        settings_middle = {
            'color': '#777777',
            'edges.width': 1,
            'opacity': 0.8,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': True,
        }

    if not settings_bounds:
        settings_bounds = {
            'color': '#999999',
            'edges.width': 1,
            'opacity': 0.5,
            'vertices.size': 0,
            'vertices.on': False,
            'edges.on': True,
            'faces.on': True,
        }

    if show_middle:
        viewer.add(middle, name="Middle/Target", settings=settings_middle)
    viewer.add(intrados, name="Intrados", settings=settings_bounds)
    viewer.add(extrados, name="Extrados", settings=settings_bounds)

    return viewer


def view_normals(shape, length=0.5):
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
        Viewer object.

    """

    viewer = view_shapes(shape, show_middle=False)
    intrados = shape.intrados
    extrados = shape.extrados

    for mesh in [intrados, extrados]:
        for key in mesh.vertices():
            n = mesh.vertex_attribute(key, 'n')
            pt0 = mesh.vertex_coordinates(key)
            pt1 = [pt0[0] + n[0]*length, pt0[1] + n[1]*length, pt0[2] + n[2]*length]
            line = Line(pt0, pt1)
            viewer.add(line, name='normal-{}'.format(key))

    return viewer


def view_fill(shape, cut_negatives=True):
    """ Viewer showing the fill

    Parameters
    ----------
    shape : obj
        Shape to plot
    cut_negatives : boll
        If true, it will hid the negative allowance on the close-to-support points (WIP).

    Returns
    ----------
    obj
        Viewer object.

    """

    middle = shape.middle
    intrados = shape.intrados
    extrados = shape.extrados
    try:
        extrados_fill = shape.extrados_fill
    except BaseException:
        print('No fill is assigned')
        extrados_fill = shape.extrados

    viewer = MultiMeshViewer()
    viewer.meshes = [middle, intrados, extrados, extrados_fill]

    return viewer
