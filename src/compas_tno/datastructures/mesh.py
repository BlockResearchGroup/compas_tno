from compas.geometry import subtract_vectors
from compas.geometry import length_vector
from compas.geometry import cross_vectors

from compas.datastructures import Mesh


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

        vertices, faces = mesh.to_vertices_and_faces()

        return cls().from_vertices_and_faces(vertices, faces)

    @classmethod
    def from_pointcloud(cls, mesh):

        # TODO: probably with CGAL

        return

    # --------------------------------------------------------------- #
    # ------------------------ULTILITIES----------------------------- #
    # --------------------------------------------------------------- #

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

        Example
        -------
        >>>

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
