from compas.datastructures import Mesh
from compas.utilities import geometric_key_xy
from numpy import array
from scipy import interpolate
from triangle import triangulate
import math


__all__ = [
    'interpolate_from_poincloud',
    'get_shape_ub',
    'get_shape_ub_pattern',
    'get_shape_ub_fill',
    'get_shape_lb',
    'get_shape_lb_pattern',
    'get_shape_middle',
    'get_shape_target',
    'get_shape_target_pattern',
    'delaunay_mesh_from_points',
    'delaunay_mesh_from_points',
    'create_mesh_from_topology_and_pointcloud',
    'create_mesh_from_topology_and_basemesh'
]


def interpolate_from_poincloud(pointcloud, XY, method='linear'):
    pointcloud_array = array(pointcloud)
    return interpolate.griddata(pointcloud_array[:, :2], pointcloud_array[:, 2], array(XY), method=method)


def get_shape_ub(shape, x, y):
    """Get the height of the extrados in the point.

    Parameters
    ----------
    x : float
        x-coordinate of the point to evaluate.
    y : float
        y-coordinate of the point to evaluate.
    Returns
    -------
    z : float
        The extrados evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.extrados.vertices_attributes('xyz'), [x, y], method=method)


def get_shape_ub_pattern(shape, XY):
    """Get the height of the extrados in a list of points.

    Parameters
    ----------
    XY : list or array
        list of the x-coordinate and y-coordinate of the points to evaluate.
    Returns
    -------
    z : float
        The extrados evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.extrados.vertices_attributes('xyz'), XY, method=method)


def get_shape_ub_fill(shape, x, y):
    """Get the height of the fill in the point.

    Parameters
    ----------
    x : float
        x-coordinate of the point to evaluate.
    y : float
        y-coordinate of the point to evaluate.
    Returns
    -------
    z : float
        The extrados evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.extrados_fill.vertices_attributes('xyz'), [x, y], method=method)


def get_shape_lb(shape, x, y):
    """Get the height of the intrados in the point.

    Parameters
    ----------
    x : float
        x-coordinate of the point to evaluate.
    y : float
        y-coordinate of the point to evaluate.
    Returns
    -------
    z : float
        The intrados evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.intrados.vertices_attributes('xyz'), [x, y], method=method)


def get_shape_lb_pattern(shape, XY):
    """Get the height of the intrados in a list of points.

    Parameters
    ----------
    XY : list or array
        list of the x-coordinate and y-coordinate of the points to evaluate.
    Returns
    -------
    z : float
        The extrados evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.intrados.vertices_attributes('xyz'), XY, method=method)


def get_shape_middle(shape, x, y):
    """Get the height of the target/middle surface in the point.

    Parameters
    ----------
    x : float
        x-coordinate of the point to evaluate.
    y : float
        y-coordinate of the point to evaluate.
    Returns
    -------
    z : float
        The middle surface evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.middle.vertices_attributes('xyz'), [x, y], method=method)


def get_shape_middle_pattern(shape, XY):
    """Get the height of the target/middle surface in a list of points.

    Parameters
    ----------
    XY : list or array
        list of the x-coordinate and y-coordinate of the points to evaluate.
    Returns
    -------
    z : float
        The extrados evaluated in the point.
    """
    method = shape.datashape.get('interpolation', 'linear')
    return interpolate_from_poincloud(shape.middle.vertices_attributes('xyz'), XY, method=method)


def get_shape_target(shape, x, y):
    """Get the height of the target/middle surface in the point."""
    return get_shape_middle(shape, x, y)


def get_shape_target_pattern(shape, XY):
    """Get the height of the target/middle surface in the point."""
    return get_shape_middle_pattern(shape, XY)


def delaunay_mesh_from_points(points):
    """Construct a Delaunay triangulation of set of vertices.

    Parameters
    ----------
    points : list
        XY(Z) coordinates of the points to triangulate.
    Returns
    -------
    obj
        MeshDos.

    """

    data = {'vertices': [point[0:2] for point in points]}
    result = triangulate(data, opts='c')
    vertices = []
    vertices_flat = []
    i = 0
    for x, y in result['vertices']:
        vertices.append([x, y, points[i][2]])
        vertices_flat.append([x, y, 0.0])
        i += 1
    faces = result['triangles']
    faces_meaningful = []

    mesh_flat = Mesh.from_vertices_and_faces(vertices_flat, faces)

    i = 0
    for fkey in mesh_flat.faces():
        if mesh_flat.face_area(fkey) > 0.001:
            faces_meaningful.append(faces[i])
        i = i+1

    mesh = Mesh.from_vertices_and_faces(vertices, faces_meaningful)  # Did this to avoid faces with area = 0, see if it is necessary

    return mesh


def create_mesh_from_topology_and_pointcloud(meshtopology, pointcloud, isnan_height=0.0):
    """
    Create a mesh based on a given topology and the heights based a pointcloud.

    Parameters
    ----------
    meshtopology : compas_tno.diagrams.FormDiagram
        Topology that is intended to keep
    mesh_base : mesh
        Mesh usually denser to base the heights on
    keep_normals : bool (True)
        Go through the process of
    Returns
    -------
    obj
        Mesh.

    """
    vertices, faces = meshtopology.to_vertices_and_faces()
    mesh = Mesh.from_vertices_and_faces(vertices, faces)
    XY = mesh.vertices_attributes('xy')
    z = interpolate_from_poincloud(pointcloud, XY)

    for i, key in enumerate(mesh.vertices()):
        if math.isnan(z[i]):
            print('Height (nan) for [x,y]:', XY[i])
            z[i] = isnan_height
        mesh.vertex_attribute(key, 'z', float(z[i]))

    return mesh


def create_mesh_from_topology_and_basemesh(meshtopology, mesh_base):
    """
    Create a mesh based on a given topology and the heights based in a base mesh (usually denser).

    Parameters
    ----------
    meshtopology : compas_tno.diagrams.FormDiagram
        Topology that is intended to keep
    mesh_base : mesh
        Mesh usually denser to base the heights on
    keep_normals : bool (True)
        Go through the process of
    Returns
    -------
    obj
        Mesh.

    """
    vertices, faces = meshtopology.to_vertices_and_faces()
    mesh = Mesh.from_vertices_and_faces(vertices, faces)
    XY = mesh.vertices_attributes('xy')
    XYZ_base = mesh_base.vertices_attributes('xyz')
    z = interpolate_from_poincloud(XYZ_base, XY)

    for i, key in enumerate(mesh.vertices()):
        mesh.vertex_attribute(key, 'z', float(z[i]))

    return mesh
