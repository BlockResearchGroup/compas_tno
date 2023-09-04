from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino
import math
from functools import partial
from math import sqrt
from math import pi
from compas_tno.rhino.diagramartist import DiagramArtist
from compas.utilities import color_to_colordict
from compas.geometry import add_vectors
from compas.geometry import scale_vector
from compas.geometry import length_vector
from compas.geometry import norm_vector
# from compas.colors import Color  # Future

from compas.utilities import is_color_rgb

colordict = partial(color_to_colordict, colorformat='rgb', normalize=False)


__all__ = ['FormArtist']


class FormArtist(DiagramArtist):
    """Artist for form diagram in TNO.

    Parameters
    ----------
    form: compas_tno.diagrams.FormDiagram
        The form diagram to draw.

    Attributes
    ----------
    mesh : :class:`compas_tno.diagrams.FormDiagram`
        The form diagram associated with the artist.
    color_compression : 3-tuple
        Default color for compression.
    color_tension : 3-tuple
        Default color for tension.
    scale_forces : float
        Scale factor for the force pipes.
    tol_forces : float
        Tolerance for force magnitudes.
    vertices : list
        The vertices to include in the drawing.
        Default is all vertices.
    edges : list
        The edges to include in the drawing.
        Default is all edges.
    faces : list
        The faces to include in the drawing.
        Default is all faces.
    vertex_xyz : dict
        The view coordinates of the vertices.
        Default is to use the actual vertex coordinates.
    vertex_color : dict
        Mapping between vertices and colors.
        Default is to use the default color for vertices.
    edge_color : dict
        Mapping between edges and colors.
        Default is to use the default color for edges.
    face_color : dict
        Mapping between faces and colors.
        Default is to use the default color for faces.
    vertex_text : dict
        Mapping between vertices and text labels.
    edge_text : dict
        Mapping between edges and text labels.
    face_text : dict
        Mapping between faces and text labels.
    """

    def __init__(self, form, layer=None):
        super(FormArtist, self).__init__(form, layer=layer)
        self.color_compression = (255, 0, 0)
        self.color_tension = (0, 255, 0)
        self.color_reaction = (125, 125, 125)
        self.default_color = (0, 0, 0)
        self.color_mesh_thrust = (0, 0, 0)
        self.color_mesh_intrados = (125, 125, 125)
        self.color_mesh_extrados = (125, 125, 125)
        self.color_mesh_middle = (0, 0, 0)
        self.color_vertex_extrados = (0, 125, 0)
        self.color_vertex_intrados = (0, 0, 255)
        self.color_faces = (0, 0, 0)
        self.scale_forces = 0.001
        self.pipes_scale = 0.025
        self.scale_line = 0.1
        self.tol_forces = 0.001
        self.radius_sphere = 0.15
        self.layer = 'FormDiagram'

        self._mesh = None
        self._vertices = None
        self._edges = None
        self._faces = None
        self._color = None
        self._vertex_xyz = None
        self._vertex_color = None
        self._vertex_text = None
        self._vertex_size = None
        self._edge_color = None
        self._edge_text = None
        self._edge_width = None
        self._face_color = None
        self._face_text = None

        self._vertexcollection = None
        self._edgecollection = None
        self._facecollection = None
        self._vertexnormalcollection = None
        self._facenormalcollection = None
        self._vertexlabelcollection = None
        self._edgelabelcollection = None
        self._facelabelcollection = None

        self.mesh = form

    @property
    def mesh(self):
        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        self._mesh = mesh
        self._vertex_xyz = None

    @property
    def vertices(self):
        if self._vertices is None:
            self._vertices = list(self.mesh.vertices())
        return self._vertices

    @vertices.setter
    def vertices(self, vertices):
        self._vertices = vertices

    @property
    def edges(self):
        if self._edges is None:
            self._edges = list(self.mesh.edges())
        return self._edges

    @edges.setter
    def edges(self, edges):
        self._edges = edges

    @property
    def faces(self):
        if self._faces is None:
            self._faces = list(self.mesh.faces())
        return self._faces

    @faces.setter
    def faces(self, faces):
        self._faces = faces

    @property
    def color(self):
        if not self._color:
            self._color = self.default_color
        return self._color

    @color.setter
    def color(self, color):
        if is_color_rgb(color):
            self._color = color

    @property
    def vertex_xyz(self):
        if self._vertex_xyz is None:
            return {vertex: self.mesh.vertex_attributes(vertex, 'xyz') for vertex in self.mesh.vertices()}
        return self._vertex_xyz

    @vertex_xyz.setter
    def vertex_xyz(self, vertex_xyz):
        self._vertex_xyz = vertex_xyz

    @property
    def vertex_color(self):
        if self._vertex_color is None:
            self._vertex_color = {vertex: self.default_vertexcolor for vertex in self.mesh.vertices()}
        return self._vertex_color

    @vertex_color.setter
    def vertex_color(self, vertex_color):
        if isinstance(vertex_color, dict):
            self._vertex_color = vertex_color
        elif is_color_rgb(vertex_color):
            self._vertex_color = {vertex: vertex_color for vertex in self.mesh.vertices()}

    @property
    def vertex_text(self):
        if self._vertex_text is None:
            self._vertex_text = {vertex: str(vertex) for vertex in self.mesh.vertices()}
        return self._vertex_text

    @vertex_text.setter
    def vertex_text(self, text):
        if text == 'key':
            self._vertex_text = {vertex: str(vertex) for vertex in self.mesh.vertices()}
        elif text == 'index':
            self._vertex_text = {vertex: str(index) for index, vertex in enumerate(self.mesh.vertices())}
        elif isinstance(text, dict):
            self._vertex_text = text

    @property
    def vertex_size(self):
        if not self._vertex_size:
            self._vertex_size = {vertex: self.default_vertexsize for vertex in self.mesh.vertices()}
        return self._vertex_size

    @vertex_size.setter
    def vertex_size(self, vertexsize):
        if isinstance(vertexsize, dict):
            self._vertex_size = vertexsize
        elif isinstance(vertexsize, (int, float)):
            self._vertex_size = {vertex: vertexsize for vertex in self.mesh.vertices()}

    @property
    def edge_color(self):
        if self._edge_color is None:
            self._edge_color = {edge: self.default_edgecolor for edge in self.mesh.edges()}
        return self._edge_color

    @edge_color.setter
    def edge_color(self, edge_color):
        if isinstance(edge_color, dict):
            self._edge_color = edge_color
        elif is_color_rgb(edge_color):
            self._edge_color = {edge: edge_color for edge in self.mesh.edges()}

    @property
    def edge_text(self):
        if self._edge_text is None:
            self._edge_text = {edge: "{}-{}".format(*edge) for edge in self.mesh.edges()}
        return self._edge_text

    @edge_text.setter
    def edge_text(self, text):
        if text == 'key':
            self._edge_text = {edge: "{}-{}".format(*edge) for edge in self.mesh.edges()}
        elif text == 'index':
            self._edge_text = {edge: str(index) for index, edge in enumerate(self.mesh.edges())}
        elif isinstance(text, dict):
            self._edge_text = text

    @property
    def edge_width(self):
        if not self._edge_width:
            self._edge_width = {edge: self.default_edgewidth for edge in self.mesh.edges()}
        return self._edge_width

    @edge_width.setter
    def edge_width(self, edgewidth):
        if isinstance(edgewidth, dict):
            self._edge_width = edgewidth
        elif isinstance(edgewidth, (int, float)):
            self._edge_width = {edge: edgewidth for edge in self.mesh.edges()}

    @property
    def face_color(self):
        if self._face_color is None:
            self._face_color = {face: self.default_facecolor for face in self.mesh.faces()}
        return self._face_color

    @face_color.setter
    def face_color(self, face_color):
        if isinstance(face_color, dict):
            self._face_color = face_color
        elif is_color_rgb(face_color):
            self._face_color = {face: face_color for face in self.mesh.faces()}

    @property
    def face_text(self):
        if self._face_text is None:
            self._face_text = {face: str(face) for face in self.mesh.faces()}
        return self._face_text

    @face_text.setter
    def face_text(self, text):
        if text == 'key':
            self._face_text = {face: str(face) for face in self.mesh.faces()}
        elif text == 'index':
            self._face_text = {face: str(index) for index, face in enumerate(self.mesh.faces())}
        elif isinstance(text, dict):
            self._face_text = text

    # def draw_edges(self, edges=None, color=None, displacement=None, layer='Thrust'):
    #     """Draw a selection of edges.

    #     Parameters
    #     ----------
    #     edges : list, optional
    #         A selection of edges to draw.
    #         The default is ``None``, in which case all edges are drawn.
    #     color : tuple or dict of tuple, optional
    #         The color specififcation for the edges.
    #         The default color is black, ``(0, 0, 0)``.

    #     Returns
    #     -------
    #     list
    #         The GUIDs of the created Rhino objects.

    #     """
    #     edges = edges or list(self.diagram.edges_where({'_is_edge': True}))
    #     vertex_xyz = self.vertex_xyz
    #     if displacement:
    #         for key in vertex_xyz:
    #             vertex_xyz[key][0] += displacement[0]
    #             vertex_xyz[key][1] += displacement[1]
    #     edge_color = colordict(color, edges, default=self.color_edges)
    #     lines = []
    #     for edge in edges:
    #         lines.append({
    #             'start': vertex_xyz[edge[0]],
    #             'end': vertex_xyz[edge[1]],
    #             'color': edge_color[edge],
    #             'name': "{}.edge.{}-{}".format(self.diagram.name, *edge)})
    #     return compas_rhino.draw_lines(lines, layer=layer, clear=False, redraw=False)

    def redraw(self):

        compas_rhino.rs.EnableRedraw(True)

        return

    def draw_mesh(self, displacement=None):
        """Draw a mesh for the thrust network.

        Parameters
        ----------
        edges : list, optional
            A selection of edges to draw.
            The default is ``None``, in which case all edges are drawn.
        displacement : list, optional
            A displacement to add mesh to the scene.

        Returns
        -------
        list
            The GUIDs of the created Rhino object.

        """
        vertices, faces = self.diagram.to_vertices_and_faces()
        return compas_rhino.draw_mesh(vertices, faces, name="Thrust", color=self.color_mesh_thrust, disjoint=True, layer="Thrust-Mesh")

    def draw_thrust(self, scale_width=None, layer="FormDiagram::Thrust", tol=1e-2):
        """Draw the thrust network as a set of lines with width proportional to the square of the carried force.

        Parameters
        ----------
        edges : list, optional
            A selection of edges to draw.
            The default is ``None``, in which case all edges are drawn.
        displacement : list, optional
            A displacement to add mesh to the scene.

        Returns
        -------
        list
            The GUIDs of the created Rhino object.

        """

        layer = layer or self.layer
        color = self.color_compression
        scale_width = scale_width or self.scale_line
        lines = []

        for u, v in self.diagram.edges_where({'_is_edge': True}):
            Xu = self.diagram.vertex_coordinates(u)
            Xv = self.diagram.vertex_coordinates(v)

            q = self.diagram.edge_attribute((u, v), 'q')
            length = self.diagram.edge_length(u, v)
            force = abs(q*length)

            if force < tol:
                continue

            radius = sqrt(force)

            thk = radius * scale_width

            lines.append({'start': Xu,
                          'end': Xv,
                          'color': color,
                          'width': thk
                          }
                         )

        return compas_rhino.draw_lines(lines, layer=layer, clear=False, redraw=True)

    def draw_cracks(self, color_intrados=(0, 0, 200), color_extrados=(0, 200, 0), layer=None, tol=10e-5):
        """Draw the intersection of the thrust network with intrados and extrados.

        Parameters
        ----------
        vertices : list, optional
            A selection of vertices to draw.
            The default is ``None``, in which case all vertices are drawn.
        displacement : list, optional
            A displacement to add mesh to the scene.
        tol : float, optional
            Acceptable error to find an intersection.

        Returns
        -------
        tuple with 2 lists
            List of the GUIDs of the created Rhino object, and the keys associated to it in the form diagram.

        """
        layer = layer or self.layer
        form = self.diagram
        intra_vertices = []
        extra_vertices = []
        intra_keys = []
        extra_keys = []
        for key in form.vertices():
            x, y, z = form.vertex_coordinates(key)
            lb = form.vertex_attribute(key, 'lb')
            ub = form.vertex_attribute(key, 'ub')
            if lb:
                if abs(z - lb) < tol:
                    intra_keys.append(key)
                    intra_vertices.append({
                        'pos': [x, y, z],
                        'name': "{}.vertex.{}".format("Intrados", key),
                        'color': self.color_vertex_intrados})
            if ub:
                if abs(z - ub) < tol:
                    extra_keys.append(key)
                    extra_vertices.append({
                        'pos': [x, y, z],
                        'name': "{}.vertex.{}".format("Extrados", key),
                        'color': self.color_vertex_extrados})

        guids_intra = compas_rhino.draw_points(intra_vertices, layer=layer + '::Intrados', color=color_intrados, clear=False, redraw=False)
        guids_extra = compas_rhino.draw_points(extra_vertices, layer=layer + '::Extrados', color=color_extrados, clear=False, redraw=False)
        return (guids_intra + guids_extra, intra_keys + extra_keys)

    # update this function
    def draw_from_attributes(self, attribute='target', name='mesh', displacement=None, color=None, join_faces=True, layer=None):
        """Draw a copy mesh of the form diagram whose height is defined based on the a given attribure.

        Parameters
        ----------
        attribute : str
            Attribute to base the heights of the network on.
        scale : float, optional
            The scaling factor for the load force vectors.
            The default value is ``0.01``
        layer : str, optional
            The layer to draw the output.

        Returns
        -------
        list
            A list with the guid of the corresponding load force vectors in Rhino.

        Notes
        -----
        The magnitude of the externally applied load at a vetex the attribute  `pz`.

        """

        faces = list(self.diagram.faces())
        layer = layer or self.layer
        vertex_xyz = self.vertex_xyz
        for key in vertex_xyz:
            z = self.diagram.vertex_attribute(key, attribute)  # Check my forms and remove this
            if isinstance(z, list):
                z = z[0]
            vertex_xyz[key][2] = z
            if displacement:
                vertex_xyz[key][0] += displacement[0]
                vertex_xyz[key][1] += displacement[1]
        face_color = colordict(color, faces, default=self.color_faces)
        facets = []
        for face in faces:
            facets.append({
                'points': [vertex_xyz[vertex] for vertex in self.diagram.face_vertices(face)],
                'name': "{}.face.{}".format(name, face),
                'color': face_color[face]})
        guids = compas_rhino.draw_faces(facets, layer=layer, clear=False, redraw=False)
        if not join_faces:
            return guids
        guid = compas_rhino.rs.JoinMeshes(guids, delete_input=True)
        compas_rhino.rs.ObjectLayer(guid, layer)
        compas_rhino.rs.ObjectName(guid, '{}'.format(name))
        if color:
            compas_rhino.rs.ObjectColor(guid, color)
        return

    def draw_loads(self, color=(255, 0, 0), scale=0.01, layer="FormDiagram::Loads", tol=1e-3):
        """Draw the externally applied loads at all vertices of the diagram.

        Parameters
        ----------
        color : list or tuple, optional
            The RGB color specification for load forces.
            The specification must be in integer format, with each component between 0 and 255.
        scale : float, optional
            The scaling factor for the load force vectors.
            The default value is ``0.01``
        layer : str, optional
            The layer to draw the output.

        Returns
        -------
        list
            A list with the guid of the corresponding load force vectors in Rhino.

        Notes
        -----
        The magnitude of the externally applied load at a vetex the attribute  `pz`.

        """
        vertex_xyz = self.vertex_xyz
        layer = layer or self.layer
        lines = []

        for vertex in self.diagram.vertices():
            a = vertex_xyz[vertex]
            pz = -1 * self.diagram.vertex_attribute(vertex, 'pz')
            load = scale_vector((0, 0, 1), scale * pz)
            b = add_vectors(a, load)
            lines.append({'start': a, 'end': b, 'color': color, 'arrow': "start"})

        return compas_rhino.draw_lines(lines, layer=layer, clear=False, redraw=False)

    def draw_reactions(self, color=None, scale=0.03, draw_as_pipes=False, layer="FormDiagram::Reactions", tol=1e-3):
        """Draw the reaction forces.

        Parameters
        ----------
        color : list or tuple
            The RGB color specification for load forces.
            The specification must be in integer format, with each component between 0 and 255.
        scale : float, optional
            The scaling factor for the load force vectors.
            The default value is ``0.01``
        layer : str, optional
            The layer to draw the output.

        Returns
        -------
        list
            A list with the guid of the corresponding load force vectors in Rhino.

        Notes
        -----
        The magnitude of the externally applied load at a vetex the attribute  `pz`.

        """
        layer = layer or self.layer
        color = color or self.color_reaction
        vertex_xyz = self.vertex_xyz
        lines = []
        cylinders = []
        labels = []
        for key in self.diagram.vertices_where({'is_fixed': True}):
            a = vertex_xyz[key]
            r = self.diagram.vertex_attributes(key, ['_rx', '_ry', '_rz'])
            res2 = r[0]**2 + r[1]**2 + r[2]**2
            res = 0
            if res2 > 0:
                res = round(math.sqrt(res2), 1)
            if not any(r):  # If receives null vector or None
                continue
            r = scale_vector(r, -scale)
            if length_vector(r) < tol:
                continue

            b = add_vectors(a, r)
            lines.append({'start': a,
                          'end': b,
                          'color': color,
                          'arrow': "start"  # 'arrow': "end"
                          }
                         )
            labels.append({'pos': [(a[0] + b[0])/2, (a[1] + b[1])/2, (a[2] + b[2])/2], 'text': str(res), 'color': (255, 0, 0)})
            if draw_as_pipes:
                force = norm_vector(self.diagram.vertex_attributes(key, ['_rx', '_ry', '_rz']))
                print(force)
                cylinders.append({
                    'start': a,
                    'end': b,
                    'radius': self.pipes_scale * sqrt(abs(force)/pi),
                    'color': color
                })

        compas_rhino.draw_labels(labels, layer="FormDiagram::ForceLabels")

        if draw_as_pipes:
            return compas_rhino.draw_cylinders(cylinders, layer=layer, clear=False, redraw=False)
        else:
            return compas_rhino.draw_lines(lines, layer=layer, clear=False, redraw=False)

    def draw_forcepipes(self, color_compression=(255, 0, 0), color_tension=(0, 0, 255),
                        compression_negative=True, pipes_scale=None, tol=1e-3, layer="FormDiagram::ForcePipes"):
        """Draw the forces in the internal edges as pipes with color and thickness matching the force value.

        Parameters
        ----------
        color_compression
        color_tension
        scale
        tol

        Returns
        -------
        list
            The GUIDs of the created Rhino objects.
        """
        layer = layer or self.layer
        vertex_xyz = self.vertex_xyz
        if not pipes_scale:
            pipes_scale = self.pipes_scale
        cylinders = []
        for edge in self.diagram.edges_where({'_is_edge': True}):
            u, v = edge
            start = vertex_xyz[u]
            end = vertex_xyz[v]
            length = self.diagram.edge_length(*edge)
            q = self.diagram.edge_attribute(edge, 'q')
            force = q * length
            if abs(force) < tol:
                continue
            radius = pipes_scale * sqrt(abs(force)/pi)
            if compression_negative:
                pipe_color = color_compression if force < 0 else color_tension
            else:
                pipe_color = color_tension if force < 0 else color_compression
            cylinders.append({
                'start': start,
                'end': end,
                'radius': radius,
                'color': pipe_color
            })
        return compas_rhino.draw_cylinders(cylinders, layer=layer, clear=False, redraw=False)

    def draw_forcelabels(self, layer="FormDiagram::ForceLabels"):
        """Draw the magnitue of the forces in the diagram.

        """
        forces = []
        for edge in self.diagram.edges_where({'_is_edge': True}):
            u, v = edge
            pt = self.diagram.edge_midpoint(u, v)
            f = self.diagram.edge_attribute((u, v), 'f')
            if abs(f) > 1e-4:
                forces.append({'text': str(round(f, 1)), 'pos': pt})

        return compas_rhino.draw_labels(forces, layer=layer)
