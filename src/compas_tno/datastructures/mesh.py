from compas.geometry import subtract_vectors
from compas.geometry import length_vector
from compas.geometry import cross_vectors
from compas.geometry import normalize_vector
from compas.geometry import is_point_in_convex_polygon_xy
from compas.geometry import centroid_points
from compas.geometry import normal_polygon
from compas.utilities import geometric_key_xy

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

    def __init__(self):  # add '_is_outside': False as default
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

    @classmethod
    def from_topology_and_mesh(cls, formdiagram, mesh_base, keep_normals=True):

        vertices, faces = formdiagram.to_vertices_and_faces()
        mesh = cls().from_vertices_and_faces(vertices, faces)
        XY = array(mesh.vertices_attributes('xy'))
        XYZ_base = array(mesh_base.vertices_attributes('xyz'))
        z = interpolate.griddata(XYZ_base[:, :2], XYZ_base[:, 2], XY, method='linear')
        count_vertex_normals = 0
        count_face_normals = 0
        if keep_normals:
            XY_project = []
            keys_project = []
            base_gkey_key_xy = mesh_base.geometric_key_xy_key()
            for i, key in enumerate(mesh.vertices()):
                x, y = XY[i]
                try:
                    key_base = base_gkey_key_xy[geometric_key_xy([x, y])]
                    n_vertex = mesh_base.vertex_normal(key_base)
                    mesh.vertex_attribute(key, 'n', n_vertex)
                    count_vertex_normals += 1
                except:
                    XY_project.append([x, y])
                    keys_project.append(key)
                    count_face_normals += 1
            if XY_project:
                normals_faces = mesh_base.get_xy_face_normals(XY_project)
                # print(len(XY_project), len(keys_project), len(normals_faces))
                for i in range(len(normals_faces)):
                    mesh.vertex_attribute(keys_project[i], 'n', normals_faces[i])

        for i, key in enumerate(mesh.vertices()):
            mesh.vertex_attribute(key, 'z', float(z[i]))

        # for key in mesh.vertices():
        #     print('n', mesh.vertex_attribute(key, 'n'))

        print('count normals {0} and faces {1}'.format(count_vertex_normals, count_face_normals))

        return mesh

    def offset_mesh(self, n=0.1, direction='up', t=0.0):
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
            normal = new_mesh.vertex_attribute(key, 'n')
            # normal = new_mesh.vertex_normal(key)
            z = new_mesh.vertex_coordinates(key)[2]
            if new_mesh.vertex_attribute(key, '_is_outside'):
                new_offset.append(t)
            else:
                deviation = math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
                if direction == 'up':
                    new_offset.append(z + n*1/deviation)
                else:
                    new_offset.append(z - n*1/deviation)

        for i, key in enumerate(new_mesh.vertices()):
            new_mesh.vertex_attribute(key, 'z', new_offset[i])

        return new_mesh

    def get_xy_face_normal(self, x, y):
        """
        Get normal based on the face normal in which the point encounters (please use vertex_normal if x, y = xi, yi where i is a vertex of the mesh) .

        Parameters
        ----------
        x : float
            X-Coordinate.
        y : float
            Y-Coordinate.
        Returns
        -------
        obj
            MeshDos.

        """

        for fkey in self.faces():
            xy_face = self.face_coordinates(fkey)
            print('xy_face', xy_face)
            if is_point_in_convex_polygon_xy([x, y], xy_face):
                return self.face_normal(fkey)

        return None

    def get_xy_face_normals(self, XY):
        """
        Get normal based on the face normal in which the points are in the plan.

        Parameters
        ----------
        XY : list
            List of Points
        Returns
        -------
        obj
            MeshDos.

        """

        face_coords = [self.face_coordinates(fkey) for fkey in self.faces()]
        normals = []
        for xy in XY:
            normals_for_pt = []
            for face_coord in face_coords:
                i = 0
                if is_point_in_convex_polygon_xy(xy, face_coord):
                    normal = normal_polygon(face_coord, unitized=False)
                    normals_for_pt.append(normal)
                    i = i+1
            # print('xy, normal', xy, normals_for_pt)
            if len(normals_for_pt) == 1:
                normal = normalize_vector(normals_for_pt[0])
            else:
                normal = normalize_vector(centroid_points(normals_for_pt))
                # print('sum vector', normal)
            normals.append(normal)
        # print(len(normals), len(XY))
        return normals

    def geometric_key_xy_key(self, precision=None):
        """Returns a dictionary that maps *geometric keys* of a certain precision
        to the keys of the corresponding vertices.

        Parameters
        ----------
        precision : str (3f)
            The float precision specifier used in string formatting.

        Returns
        -------
        dict
            A dictionary of geometric key-key pairs.

        """
        XY = self.vertices_attributes('xy')
        return {geometric_key_xy(XY[i], precision): key for i, key in enumerate(self.vertices())}


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
