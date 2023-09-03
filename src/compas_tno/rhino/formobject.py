from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino
from compas_tno.rhino.diagramobject import DiagramObject
from compas.geometry import add_vectors
from compas.geometry import scale_vector
from compas.geometry import norm_vector
from compas.geometry import centroid_points


__all__ = ['FormObject']


class FormObject(DiagramObject):
    """A form object represents a form diagram in the Rhino model space.
    """

    SETTINGS = {
        'show.vertices': False,
        'show.supports': True,
        'show.edges': True,
        'show.vertexlabels': False,
        'show.vertexloads': False,
        'show.edgelabels': False,
        'show.reactionlabels': False,
        'show.cracks': False,
        'show.vertices_bound': False,
        'show.forcecolors': False,
        'show.forcelabels': False,
        'show.forcepipes': False,
        'show.loadvectors': False,
        'show.reactionvectors': False,

        'is.compression_negative': True,

        'color.vertices': (0, 0, 0),
        'color.vertexlabels': (255, 255, 255),
        'color.vertexloads': (255, 255, 255),
        'color.vertices:is_fixed': (255, 0, 0),
        'color.vertices:upper_bound': (0, 200, 0),
        'color.vertices:lower_bound': (0, 0, 200),
        'color.edges': (0, 0, 0),
        'color.edges:is_ind': (0, 0, 0),
        'color.edges:is_external': (0, 255, 0),
        'color.edges:is_reaction': (0, 255, 0),
        'color.edges:is_load': (0, 255, 0),
        'color.reactionvectors': (255, 0, 0),
        'color.loadvectors': (0, 200, 200),
        'color.faces': (210, 210, 210),
        'color.compression': (255, 0, 0),
        'color.tension': (0, 0, 255),

        'scale.forcepipes': 0.001,
        'scale.vectors': 0.01,

        'tol.edges': 1e-3,
        'tol.forces': 1e-3,
    }

    def __init__(self, diagram, *args, **kwargs):
        super(FormObject, self).__init__(diagram, *args, **kwargs)
        self.settings.update(FormObject.SETTINGS)
        settings = kwargs.get('settings') or {}
        if settings:
            self.settings.update(settings)
        self._guid_pipes = {}
        self._guid_cracks = {}
        self._guid_vectors = {}
        self._guid_reactions = {}
        self._guids = {}

    @property
    def guids(self):
        guids = super(FormObject, self).guids
        guids += list(self.guid_pipes.keys())
        guids += list(self.guid_cracks.keys())
        guids += list(self.guid_vectors.keys())
        guids += list(self.guid_reactions.keys())
        return guids

    @guids.setter
    def guids(self, values):
        self._guids = dict(values)

    @property
    def guid_pipes(self):
        """Map between Rhino object GUIDs and form diagram edge force identifiers."""
        return self._guid_pipes

    @guid_pipes.setter
    def guid_pipes(self, values):
        self._guid_pipes = dict(values)

    @property
    def guid_cracks(self):
        """Map between Rhino object GUIDs and thrust network cracks."""
        return self._guid_cracks

    @guid_cracks.setter
    def guid_cracks(self, values):
        self._guid_cracks = dict(values)

    @property
    def guid_vectors(self):
        """Map between Rhino object GUIDs and thrust network vectors to represent loads."""
        return self._guid_vectors

    @guid_vectors.setter
    def guid_vectors(self, values):
        self._guid_vectors = dict(values)

    @property
    def guid_reactions(self):
        """Map between Rhino object GUIDs and thrust network vectors to represent reactions."""
        return self._guid_reactions

    @guid_reactions.setter
    def guid_reactions(self, values):
        self._guid_reactions = dict(values)

    def clear(self):
        super(FormObject, self).clear()
        compas_rhino.delete_objects(self.guids, purge=True)
        self._guid_pipes = {}
        self._guid_cracks = {}
        self._guid_vectors = {}
        self._guid_reactions = {}
        self._guids = {}

    def draw(self):
        """Draw the form diagram.

        The visible components, display properties and visual style of the form diagram
        drawn by this method can be fully customised using the configuration items
        in the settings dict: ``FormArtist.settings``.

        The method will clear the scene of any objects it has previously drawn
        and keep track of any newly created objects using their GUID.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        self.clear()
        if not self.visible:
            return

        self.artist.vertex_xyz = self.vertex_xyz

        guids = []
        vertices = []

        # vertices (including supports)
        if self.settings['show.vertices']:
            vertices = list(self.diagram.vertices())
            color = {}
            color.update({vertex: self.settings['color.vertices'] for vertex in vertices})
            color.update({vertex: self.settings['color.vertices:is_fixed'] for vertex in self.diagram.vertices_where({'is_fixed': True})})
            guids = self.artist.draw_vertices(color=color)
            self.guid_vertex = zip(guids, vertices)

            # vertex labels
            if self.settings['show.vertexlabels']:
                text = {vertex: index for index, vertex in enumerate(vertices)}
                color = {}
                color.update({vertex: self.settings['color.vertexlabels'] for vertex in vertices})
                color.update({vertex: self.settings['color.vertices:is_fixed'] for vertex in self.diagram.vertices_where({'is_fixed': True})})
                guids = self.artist.draw_vertexlabels(text=text, color=color)
                self.guid_vertexlabel = zip(guids, vertices)

        # supports only
        elif self.settings['show.supports']:
            vertices = list(self.diagram.vertices_where({'is_fixed': True}))
            color = {}
            color.update({vertex: self.settings['color.vertices:is_fixed'] for vertex in self.diagram.vertices_where({'is_fixed': True})})
            guids = self.artist.draw_vertices(vertices=vertices, color=color)
            self.guid_vertex = zip(guids, vertices)

        # upper and lower bounds
        if self.settings['show.vertices_bound']:
            points = []
            vertices = []
            for key in self.diagram.vertices():
                lb = self.diagram.vertex_attribute(key, 'lb')
                ub = self.diagram.vertex_attribute(key, 'ub')
                x, y, z = self.diagram.vertex_coordinates(key)
                if lb:
                    vertices.append(key)
                    points.append({
                        'pos': [x, y, lb],
                        'color': self.settings['color.vertices:lower_bound']
                    })
                if ub:
                    vertices.append(key)
                    points.append({
                        'pos': [x, y, ub],
                        'color': self.settings['color.vertices:upper_bound']
                    })

            guids = compas_rhino.draw_points(points, layer="TNO::Shape::Bounds", clear=False, redraw=False)
            self.guid_cracks = zip(guids, vertices)

        if self.settings['show.cracks']:
            guids, vertices = self.artist.draw_cracks(
                color_intrados=self.settings['color.vertices:lower_bound'],
                color_extrados=self.settings['color.vertices:upper_bound'],
                layer='TNO::Shape'
            )
            self.guid_cracks = zip(guids, vertices)

        # vertex loads
        if self.settings['show.vertexloads']:
            vertices = list(self.diagram.vertices())
            text = {vertex: '{0:.1f}'.format(self.diagram.vertex_attribute(vertex, 'pz')) for vertex in self.diagram.vertices()}
            color = {}
            color.update({vertex: self.settings['color.vertexloads'] for vertex in vertices})
            color.update({vertex: self.settings['color.vertices:is_fixed'] for vertex in self.diagram.vertices_where({'is_fixed': True})})
            guids = self.artist.draw_vertexlabels(text=text, color=color)
            self.guid_vertexlabel = zip(guids, vertices)

        # edges
        if self.settings['show.edges']:
            edges = list(self.diagram.edges_where({'_is_edge': True}))
            color = {}
            color.update({edge: self.settings['color.edges'] for edge in edges})
            color.update({edge: self.settings['color.edges:is_external'] for edge in self.diagram.edges_where({'is_external': True})})
            color.update({edge: self.settings['color.edges:is_load'] for edge in self.diagram.edges_where({'is_load': True})})
            color.update({edge: self.settings['color.edges:is_reaction'] for edge in self.diagram.edges_where({'is_reaction': True})})
            color.update({edge: self.settings['color.edges:is_ind'] for edge in self.diagram.edges_where({'is_ind': True})})

            # force colors
            if self.settings['show.forcecolors']:
                tol = self.settings['tol.forces']
                for edge in self.diagram.edges_where({'is_external': False}):
                    if self.diagram.edge_attribute(edge, 'f') > + tol:
                        color[edge] = self.settings['color.tension']
                    elif self.diagram.edge_attribute(edge, 'f') < - tol:
                        color[edge] = self.settings['color.compression']

            print('draw_edges run')
            guids = self.artist.draw_edges(edges=edges, color=color)
            self.guid_edge = zip(guids, edges)

            guid_edgelabel = []

            # edge labels
            if self.settings['show.edgelabels']:
                text = {edge: index for index, edge in enumerate(edges)}
                color = {}
                color.update({edge: self.settings['color.edges'] for edge in edges})
                color.update({edge: self.settings['color.edges:is_external'] for edge in self.diagram.edges_where({'is_external': True})})
                color.update({edge: self.settings['color.edges:is_load'] for edge in self.diagram.edges_where({'is_load': True})})
                color.update({edge: self.settings['color.edges:is_reaction'] for edge in self.diagram.edges_where({'is_reaction': True})})
                color.update({edge: self.settings['color.edges:is_ind'] for edge in self.diagram.edges_where({'is_ind': True})})

                # force colors
                if self.settings['show.forcecolors']:
                    tol = self.settings['tol.forces']
                    for edge in self.diagram.edges_where({'is_external': False}):
                        if self.diagram.edge_attribute(edge, 'f') > + tol:
                            color[edge] = self.settings['color.tension']
                        elif self.diagram.edge_attribute(edge, 'f') < - tol:
                            color[edge] = self.settings['color.compression']

                guids = self.artist.draw_edgelabels(text=text, color=color)
                guid_edgelabel += zip(guids, edges)

            # force labels
            if self.settings['show.forcelabels']:
                text = {}
                for index, edge in enumerate(self.diagram.edges_where({'is_external': True})):
                    f = self.diagram.edge_attribute(edge, 'f')
                    if f != 0.0:
                        text[edge] = "{:.4g}kN".format(abs(f))
                color = {}
                # color.update({edge: self.settings['color.edges'] for edge in self.diagram.edges()})
                color.update({edge: self.settings['color.edges:is_external'] for edge in self.diagram.edges_where({'is_external': True})})
                color.update({edge: self.settings['color.edges:is_load'] for edge in self.diagram.edges_where({'is_load': True})})
                color.update({edge: self.settings['color.edges:is_reaction'] for edge in self.diagram.edges_where({'is_reaction': True})})
                color.update({edge: self.settings['color.edges:is_ind'] for edge in self.diagram.edges_where({'is_ind': True})})

                guids = self.artist.draw_edgelabels(text=text, color=color)
                guid_edgelabel += zip(guids, edges)

            # reaction labels
            if self.settings['show.reactionlabels']:
                vertices = list(self.diagram.vertices_where({'is_fixed': True}))
                labels = []
                for vertex in vertices:
                    a = self.vertex_xyz[vertex]
                    r = self.diagram.vertex_attributes(vertex, ['_rx', '_ry', '_rz'])
                    if not any(r):  # If receives null vector or None
                        text = "0.0 kN"
                        b = a
                    else:
                        text = "{:.4g}kN".format(norm_vector(r))
                        r = scale_vector(r, -1 * self.settings['scale.vectors'])
                        b = add_vectors(a, r)
                    labels.append({
                        'pos': centroid_points([a, b]),
                        'name': "{}.reactionlabel.{}".format(self.diagram.name, vertex),
                        'color': self.settings['color.reactionvectors'],
                        'text': text})

                guids = compas_rhino.draw_labels(labels, layer=self.layer, clear=False, redraw=False)
                guid_edgelabel += zip(guids, vertices)

            self.guid_edgelabel = guid_edgelabel

        # force pipes
        if self.settings['show.forcepipes']:
            edges = list(self.diagram.edges_where({'_is_edge': True}))
            guids = self.artist.draw_forcepipes(
                color_compression=self.settings['color.compression'],
                color_tension=self.settings['color.tension'],
                compression_negative=self.settings['is.compression_negative'],
                scale=self.settings['scale.forcepipes'],
                tol=self.settings['tol.forces'])

            self.guid_pipes = zip(guids, edges)

        # load vectors
        if self.settings['show.loadvectors']:
            vertices = list(self.diagram.vertices())
            guids = self.artist.draw_loads(
                color=self.settings['color.loadvectors'],
                scale=self.settings['scale.vectors'],
                layer="TNO::Loads",
                tol=self.settings['tol.forces']
            )

            self.guid_vectors = zip(guids, edges)

        # reaction vectors
        if self.settings['show.reactionvectors']:
            vertices = list(self.diagram.vertices_where({'is_fixed': True}))
            guids = self.artist.draw_reactions(
                color=self.settings['color.reactionvectors'],
                scale=self.settings['scale.vectors'],
                layer="TNO::Reactions",
                tol=self.settings['tol.forces']
            )

            self.guid_reactions = zip(guids, edges)

        self.guid = self.guid_edge

        self.redraw()
