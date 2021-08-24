from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino
from compas_rhino.artists import BaseArtist
from compas.utilities import pairwise
from compas.geometry import centroid_polygon


__all__ = ['ShapeArtist']


class ShapeArtist(BaseArtist):
    """Artist for shape in TNO.

    Parameters
    ----------
    shape: compas_tno.shapes.Shape
        The shape to draw.

    """

    def __init__(self, shape, scene=None, name=None, layer=None, visible=True, settings=None):
        super(ShapeArtist, self).__init__()
        self._shape = shape
        self.intrados_layer = "Shape::Intrados"
        self.extrados_layer = "Shape::Extrados"
        self.middle_layer = "Shape::Middle"

    @property
    def shape(self):
        """The shape associated with the artist."""
        return self._shape

    @shape.setter
    def shape(self, shape):
        self._shape = shape

    def draw_intrados(self, color=(0, 0, 0), displacement=None, layer=None):

        layer = layer or self.intrados_layer
        intrados = self.shape.intrados

        return self.draw_shape_mesh(intrados, color=color, layer=layer)

    def draw_extrados(self, color=(0, 0, 0), displacement=None, layer=None):

        layer = layer or self.extrados_layer
        extrados = self.shape.extrados

        return self.draw_shape_mesh(extrados, color=color, layer=layer)

    def draw_middle(self, color=(0, 0, 0), displacement=None, layer=None):

        layer = layer or self.middle_layer
        middle = self.shape.middle

        return self.draw_shape_mesh(middle, color=color, layer=layer)

    def draw_shape_mesh(self, mesh, color=(0, 0, 0), disjoint=False, layer=None):
        """Draw the meshes of the shape object as a consolidated RhinoMesh.

        Parameters
        ----------
        color : tuple, optional
            The color of the mesh.
            Default is black, ``(0, 0, 0)``.
        disjoint : bool, optional
            Draw the faces of the mesh with disjoint vertices.
            Default is ``False``.

        Returns
        -------
        list
            The GUIDs of the created Rhino objects.

        Notes
        -----
        The mesh should be a valid Rhino Mesh object, which means it should have only triangular or quadrilateral faces.
        Faces with more than 4 vertices will be triangulated on-the-fly.

        """
        vertex_index = mesh.key_index()
        vertices = [mesh.vertex_coordinates(vertex) for vertex in mesh.vertices()]
        faces = [[vertex_index[vertex] for vertex in mesh.face_vertices(face)] for face in mesh.faces()]
        new_faces = []
        for face in faces:
            f = len(face)
            if f == 3:
                new_faces.append(face + face[-1:])
            elif f == 4:
                new_faces.append(face)
            elif f > 4:
                centroid = len(vertices)
                vertices.append(centroid_polygon([vertices[index] for index in face]))
                for a, b in pairwise(face + face[0:1]):
                    new_faces.append([centroid, a, b, b])
            else:
                continue
        name = "{}".format(mesh.name)
        guid = compas_rhino.draw_mesh(vertices, new_faces, layer=layer, name=name, color=color, disjoint=disjoint)
        return [guid]
