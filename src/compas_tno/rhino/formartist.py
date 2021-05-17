from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino

from functools import partial
from math import fabs
from compas_tno.rhino.diagramartist import DiagramArtist
from compas.utilities import color_to_colordict

import rhinoscriptsyntax as rs

colordict = partial(color_to_colordict, colorformat='rgb', normalize=False)


__all__ = ['FormArtist']


class FormArtist(DiagramArtist):
    """Artist for form diagram in AGS.

    Parameters
    ----------
    form: compas_ags.diagrams.FormDiagram
        The form diagram to draw.

    Attributes
    ----------
    color_compression : 3-tuple
        Default color for compression.
    color_tension : 3-tuple
        Default color for tension.
    scale_forces : float
        Scale factor for the force pipes.
    tol_forces : float
        Tolerance for force magnitudes.
    """

    def __init__(self, form, layer=None):
        super(FormArtist, self).__init__(form, layer=layer)
        self.color_compression = (255, 0, 0)
        self.color_tension = (0, 255, 0)
        self.color_mesh_thrust = (0, 0, 0)
        self.color_mesh_intrados = (0, 0, 0)
        self.color_mesh_extrados = (0, 0, 0)
        self.color_mesh_middle = (0, 0, 0)
        self.color_vertex_extrados = (0, 255, 0)
        self.color_vertex_intrados = (0, 0, 255)
        self.color_faces = (0, 0, 0)
        self.scale_forces = 0.001
        self.tol_forces = 0.001
        self.radius_sphere = 0.15

    def draw_edges(self, edges=None, color=None, displacement=None, layer='Thrust'):
        """Draw a selection of edges.

        Parameters
        ----------
        edges : list, optional
            A selection of edges to draw.
            The default is ``None``, in which case all edges are drawn.
        color : tuple or dict of tuple, optional
            The color specififcation for the edges.
            The default color is black, ``(0, 0, 0)``.

        Returns
        -------
        list
            The GUIDs of the created Rhino objects.

        """
        edges = edges or list(self.diagram.edges())
        vertex_xyz = self.vertex_xyz
        if displacement:
            for key in vertex_xyz:
                vertex_xyz[key][0] += displacement[0]
                vertex_xyz[key][1] += displacement[1]
        edge_color = colordict(color, edges, default=self.color_edges)
        lines = []
        for edge in edges:
            lines.append({
                'start': vertex_xyz[edge[0]],
                'end': vertex_xyz[edge[1]],
                'color': edge_color[edge],
                'name': "{}.edge.{}-{}".format(self.diagram.name, *edge)})
        return compas_rhino.draw_lines(lines, layer=layer, clear=False, redraw=False)

    def draw_forcepipes(self, color_compression=None, color_tension=None, scale=None, tol=None, displacement=None, layer='Pipes'):
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
        color_compression = color_compression or self.color_compression
        color_tension = color_tension or self.color_tension
        scale = scale or self.scale_forces
        tol = tol or self.tol_forces
        vertex_xyz = self.vertex_xyz
        if displacement:
            for key in vertex_xyz:
                vertex_xyz[key][0] += displacement[0]
                vertex_xyz[key][1] += displacement[1]
        edges = []
        pipes = []
        for edge in self.diagram.edges():
            length = self.diagram.edge_length(*edge)
            q = self.diagram.edge_attribute(edge, 'q')
            force = q * length
            if force > tol:
                radius = fabs(scale * force)
            else:
                continue
            edges.append(edge)
            color = color_compression if force > 0 else color_tension
            pipes.append({'points': [vertex_xyz[edge[0]], vertex_xyz[edge[1]]],
                          'color': color,
                          'radius': radius,
                          'name': "{}.force.{}-{}".format(self.diagram.name, *edge)})
        return compas_rhino.draw_pipes(pipes, layer=layer, clear=False, redraw=False)

    def draw_thrust(self, displacement=None):
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

    def draw_cracks(self, vertices=None, displacement=None, tol_cracks=10e-5, layer='Cracks', spheres=False):
        """Draw the intersection of the thrust network with intrados and extrados.

        Parameters
        ----------
        vertices : list, optional
            A selection of vertices to draw.
            The default is ``None``, in which case all vertices are drawn.
        displacement : list, optional
            A displacement to add mesh to the scene.
        tol_cracks : float, optional
            Acceptable error to find an intersection.

        Returns
        -------
        list
            The GUIDs of the created Rhino object.

        """
        layer = layer or self.layer
        form = self.diagram
        intra_vertices = []
        extra_vertices = []
        for key in form.vertices():
            x, y, z = form.vertex_coordinates(key)
            if displacement:
                x += displacement[0]
                y += displacement[1]
            lb = form.vertex_attribute(key, 'lb')
            ub = form.vertex_attribute(key, 'ub')
            if isinstance(lb, list):
                lb = lb[0]
                ub = ub[0]
            if lb:
                if abs(z - lb) < tol_cracks:
                    if spheres:
                        sphere = compas_rhino.rs.AddSphere([x, y, z], self.radius_sphere)
                        rs.ObjectColor(sphere, self.color_vertex_intrados)
                        rs.ObjectLayer(sphere, layer)
                    else:
                        intra_vertices.append({
                            'pos': [x, y, z],
                            'name': "{}.vertex.{}".format("Intrados", key),
                            'color': self.color_vertex_intrados})
            if ub:
                if abs(z - ub) < tol_cracks:
                    if spheres:
                        sphere = compas_rhino.rs.AddSphere([x, y, z], self.radius_sphere)
                        rs.ObjectColor(sphere, self.color_vertex_extrados)
                        rs.ObjectLayer(sphere, layer)
                    else:
                        extra_vertices.append({
                            'pos': [x, y, z],
                            'name': "{}.vertex.{}".format("Extrados", key),
                            'color': self.color_vertex_extrados})

        compas_rhino.draw_points(intra_vertices, layer=layer, clear=False, redraw=False)
        compas_rhino.draw_points(extra_vertices, layer=layer, clear=False, redraw=False)
        return

    def draw_from_attributes(self, attribute='target', name='mesh', displacement=None, color=None, join_faces=True, layer=None):

        faces = list(self.diagram.faces())
        layer = layer or self.layer
        vertex_xyz = self.vertex_xyz
        for key in vertex_xyz:
            z = self.diagram.vertex_attribute(key, attribute)  # Check my forms and remove this
            if type(z) == list:
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

    def draw_intrados(self, displacement=None, layer='Intrados'):

        self.draw_from_attributes(attribute='lb', name='Intrados', color=self.color_mesh_intrados, displacement=displacement, layer=layer)

        return

    def draw_extrados(self, displacement=None, layer='Extrados'):

        self.draw_from_attributes(attribute='ub', name='Extrados', color=self.color_mesh_extrados, displacement=displacement, layer=layer)

        return

    def draw_middle(self, displacement=None, layer='Middle'):

        self.draw_from_attributes(attribute='target', name='Middle', color=self.color_mesh_middle, displacement=displacement, layer=layer)

        return

    def draw_reactions(self, layer='Reactions', displacement=None, TextDot=False):

        layer = layer or self.layer
        # rs.CurrentLayer(layer)
        for key in self.diagram.vertices_where({'is_fixed': True}):
            xb, yb, zb = self.diagram.vertex_coordinates(key)
            if displacement:
                xb += displacement[0]
                yb += displacement[1]
            rx = self.diagram.vertex_attribute(key, '_rx')
            ry = self.diagram.vertex_attribute(key, '_ry')
            rz = self.diagram.vertex_attribute(key, '_rz')
            norm = (rx ** 2 + ry ** 2 + rz ** 2) ** (1/2)
            if rz < 0.0 and norm > 0.0:
                sp = [xb, yb, zb]
                dz = rz/norm
                mult = zb/dz
                dz *= mult
                dx_ = mult * rx/norm
                dy_ = mult * ry/norm
                ep = [xb - dx_, yb - dy_, zb - dz]
                id = compas_rhino.rs.AddLine(sp, ep)
                compas_rhino.rs.ObjectName(id, str(norm))
                compas_rhino.rs.ObjectLayer(id, layer)
                # rs.CurrentLayer(reac_val)
                if TextDot:
                    rs.AddTextDot('({0:.1f};{1:.1f})'.format(rx, ry), sp, layer=layer)
                # rs.CurrentLayer(reac_layer)

        return

    def redraw(self):

        rs.EnableRedraw(True)

        return
