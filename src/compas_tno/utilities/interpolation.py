import math
from typing import Annotated
from typing import Optional
from typing import Type

from numpy import asarray
from scipy import interpolate

from compas.datastructures import Mesh
from compas.geometry import delaunay_triangulation


def interpolate_from_pointcloud(pointcloud, XY, method="linear"):
    pointcloud_array = asarray(pointcloud)
    return interpolate.griddata(pointcloud_array[:, :2], pointcloud_array[:, 2], asarray(XY), method=method)


def mesh_from_pointcloud(points, cls=None) -> Mesh:
    """Construct a Delaunay triangulation of set of vertices.

    Parameters
    ----------
    points : list
        XY(Z) coordinates of the points to triangulate.

    Returns
    -------
    mesh : MeshDos
        Delaunay mesh created.

    """
    mesh: Mesh = delaunay_triangulation(points)

    XY = mesh.vertices_attributes("xy")
    z = interpolate_from_pointcloud(points, XY)

    for i, key in enumerate(mesh.vertices()):
        mesh.vertex_attribute(key, "z", z[i])

    if cls:
        return mesh.copy(cls=cls)

    return mesh


def create_mesh_from_topology_and_pointcloud(
    topology: Mesh,
    pointcloud: list[Annotated[list[float], 3]],
    isnan_height: float = 0.0,
    cls: Optional[Type[Mesh]] = None,
) -> Mesh:
    """
    Create a mesh based on a given topology and the heights based in a pointcloud.

    Parameters
    ----------
    topology : FormDiagram
        Topology that is intended to keep
    pointcloud : list[point]
        Mesh usually denser to use as base to interpolate the heights
    isnan_height : bool, optional
        Value if the height is nan, by default is 0.0

    Returns
    -------
    obj
        Mesh.

    """
    cls = cls or Mesh

    vertices, faces = topology.to_vertices_and_faces()

    mesh = cls.from_vertices_and_faces(vertices, faces)

    XY = mesh.vertices_attributes("xy")
    z = interpolate_from_pointcloud(pointcloud, XY)

    for i, key in enumerate(mesh.vertices()):
        if math.isnan(z[i]):
            print("Height (nan) for [x,y]:", XY[i])
            z[i] = isnan_height

        mesh.vertex_attribute(key, "z", float(z[i]))

    return mesh


def create_mesh_from_topology_and_basemesh(
    meshtopology: Mesh,
    mesh_base: Mesh,
    cls: Optional[Type[Mesh]],
) -> Mesh:
    """
    Create a mesh based on a given topology and the heights based in a base mesh (usually denser).

    Parameters
    ----------
    meshtopology : compas_tno.diagrams.FormDiagram
        Topology that is intended to keep
    mesh_base : mesh
        Mesh usually denser to use as base to interpolate the heights

    Returns
    -------
    obj
        Mesh.

    """
    cls = cls or Mesh

    vertices, faces = meshtopology.to_vertices_and_faces()

    mesh = cls.from_vertices_and_faces(vertices, faces)
    XY = mesh.vertices_attributes("xy")
    XYZ_base = mesh_base.vertices_attributes("xyz")
    z = interpolate_from_pointcloud(XYZ_base, XY)

    for i, key in enumerate(mesh.vertices()):
        mesh.vertex_attribute(key, "z", float(z[i]))

    return mesh


# def delaunay_mesh_from_points(points):
#     """Construct a Delaunay triangulation of set of vertices.

#     Parameters
#     ----------
#     points : list
#         XY(Z) coordinates of the points to triangulate.

#     Returns
#     -------
#     mesh : MeshDos
#         Delaunay mesh created.

#     """
#     raise NotImplementedError

#     # from triangle import triangulate  # check to do it with scipy
#     # data = {'vertices': [point[0:2] for point in points]}
#     # result = triangulate(data, opts='c')
#     # vertices = []
#     # vertices_flat = []
#     # i = 0
#     # for x, y in result['vertices']:
#     #     vertices.append([x, y, points[i][2]])
#     #     vertices_flat.append([x, y, 0.0])
#     #     i += 1
#     # faces = result['triangles']
#     # faces_meaningful = []

#     # mesh_flat = MeshDos.from_vertices_and_faces(vertices_flat, faces)

#     # i = 0
#     # for fkey in mesh_flat.faces():  # Did this to avoid faces with area = 0, see if it is necessary
#     #     if mesh_flat.face_area(fkey) > 0.001:
#     #         faces_meaningful.append(faces[i])
#     #     i = i+1

#     # mesh = MeshDos.from_vertices_and_faces(vertices, faces_meaningful)

#     # return mesh
