from compas.geometry import subtract_vectors
from compas.geometry import length_vector
from compas.geometry import cross_vectors

from compas.datastructures import Mesh

from triangle import triangulate

from scipy import interpolate

from numpy import array

import math


__all__ = ['MeshDos']


class MeshDos(Mesh):

    """The ``MeshDos`` inherits compas ``Mesh`` object and add some features to create Intrados and Extrados.

    Notes
    -----
    A ``MeshDos`` has the following constructor functions

    *   ``from_library`` : Construct the shape from a dictionary with instructions, the library supports the creation of parametric arches, domes, and vaults.
    *   ``from_rhinomesh`` : Construct Extrados, Intrados and Middle surfaces from RhinoMeshes.
    *   ``from_rhinosurface`` : Construct Extrados, Intrados and Middle surfaces using the U and V isolines.
    *   ``from_pointcloud`` : Construct Extrados, Intrados and Middle surfaces interoplating a pointcloud.

    """

    __module__ = 'compas_tno.datastructures'

    def __init__(self):
        super(Mesh, self).__init__()

    @classmethod
    def from_mesh(cls, mesh):
        """
        Construct object from an existing COMPAS mesh.

        Parameters
        ----------
        obj : BaseMesh
            Base mesh.
        Returns
        -------
        obj
            MeshDos.

        """

        vertices, faces = mesh.to_vertices_and_faces()

        return cls().from_vertices_and_faces(vertices, faces)

    @classmethod
    def from_points_delaunay(cls, points):
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
        i = 0
        for x, y in result['vertices']:
            vertices.append([x, y, points[i][2]])
            i += 1
        faces = result['triangles']

        mesh = cls().from_vertices_and_faces(vertices, faces)

        return mesh


    @classmethod
    def from_formdiagram_attribute(cls, formdiagram, attribute='lb'):

        vertices, faces = formdiagram.to_vertices_and_faces()
        mesh = cls().from_vertices_and_faces(vertices, faces)

        for key in formdiagram.vertices():
            z_ = formdiagram.vertex_attribute(key, attribute)
            mesh.vertex_attribute(key, 'z', z_)

        return mesh

    @classmethod
    def from_topology_and_pointcloud(cls, formdiagram, pointcloud):

        vertices, faces = formdiagram.to_vertices_and_faces()
        mesh = cls().from_vertices_and_faces(vertices, faces)
        XY = array(mesh.vertices_attributes('xy'))
        z = interpolate.griddata(pointcloud[:, :2], pointcloud[:, 2], XY, method='linear')

        for i, key in enumerate(mesh.vertices()):
            mesh.vertex_attribute(key, 'z', float(z[i]))

        return mesh

    def offset_mesh(self, n=0.1, direction='up'):
        """
        Offset the mesh upwards considering it as intrados.

        Parameters
        ----------
        n : float
            Distance to offset the mesh.
        Returns
        -------
        obj
            MeshDos.

        """

        new_mesh = self.copy()
        new_offset = []

        for key in new_mesh.vertices():
            normal = new_mesh.vertex_normal(key)
            z = new_mesh.vertex_coordinates(key)[2]
            deviation = math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
            if direction == 'up':
                new_offset.append(z + n*1/deviation)
            else:
                new_offset.append(z - n*1/deviation)

        for i, key in enumerate(new_mesh.vertices()):
            new_mesh.vertex_attribute(key, 'z', new_offset[i])

        return new_mesh


    def vertex_projected_area(self, key):
        """Compute the projected tributary area of a vertex.

        Parameters
        ----------
        key : int
            The identifier of the vertex.

        Returns
        -------
        float
            The projected tributary area.

        """

        area = 0.

        p0 = self.vertex_coordinates(key)
        p0[2] = 0

        for nbr in self.halfedge[key]:
            p1 = self.vertex_coordinates(nbr)
            p1[2] = 0
            v1 = subtract_vectors(p1, p0)

            fkey = self.halfedge[key][nbr]
            if fkey is not None:
                p2 = self.face_centroid(fkey)
                p2[2] = 0
                v2 = subtract_vectors(p2, p0)
                area += length_vector(cross_vectors(v1, v2))

            fkey = self.halfedge[nbr][key]
            if fkey is not None:
                p3 = self.face_centroid(fkey)
                p3[2] = 0
                v3 = subtract_vectors(p3, p0)
                area += length_vector(cross_vectors(v1, v3))

        return 0.25 * area
