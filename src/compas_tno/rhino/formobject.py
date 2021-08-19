from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino
from compas_tno.rhino.diagramobject import DiagramObject


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
        'show.cracks': False,
        'show.forcecolors': False,
        'show.forcelabels': False,
        'show.forcepipes': False,
        'show.vertices_bound': False,

        'color.vertices': (0, 0, 0),
        'color.vertexlabels': (255, 255, 255),
        'color.vertexloads': (255, 255, 255),
        'color.vertices:is_fixed': (255, 0, 0),
        'color.vertices:upper_bound': (0, 255, 0),
        'color.vertices:lower_bound': (0, 0, 255),
        'color.edges': (0, 0, 0),
        'color.edges:is_ind': (0, 0, 0),
        'color.edges:is_external': (0, 255, 0),
        'color.edges:is_reaction': (0, 255, 0),
        'color.edges:is_load': (0, 255, 0),
        'color.faces': (210, 210, 210),
        'color.compression': (255, 0, 0),
        'color.tension': (0, 0, 255),
        'scale.forces': None,

        'tol.edges': 1e-3,
        'tol.forces': 1e-3,
    }

    def __init__(self, diagram, *args, **kwargs):
        super(FormObject, self).__init__(diagram, *args, **kwargs)
        self.settings.update(FormObject.SETTINGS)
        settings = kwargs.get('settings') or {}
        if settings:
            self.settings.update(settings)
        self._guid_force = {}

    @property
    def guids(self):
        guids = super(FormObject, self).guids
        guids += list(self.guid_force.keys())
        return guids

    @property
    def guid_force(self):
        """Map between Rhino object GUIDs and form diagram edge force identifiers."""
        return self._guid_force

    @guid_force.setter
    def guid_force(self, values):
        self._guid_force = dict(values)

    def clear(self):
        super(FormObject, self).clear()
        compas_rhino.delete_objects(self.guids, purge=True)
        self._guid_force = {}

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
            points_bounds = []
            vertices_bounds = []
            for key in self.diagram.vertices():
                lb = self.diagram.vertex_attribute(key, 'lb')
                ub = self.diagram.vertex_attribute(key, 'ub')
                x, y, z = self.diagram.vertex_coordinates(key)
                if lb:
                    vertices_bounds.append(key)
                    points_bounds.append({
                        'pos': [x, y, lb],
                        'color': self.settings['color.vertices:lower_bound']
                    })
                if ub:
                    vertices_bounds.append(key)
                    points_bounds.append({
                        'pos': [x, y, ub],
                        'color': self.settings['color.vertices:upper_bound']
                    })

            if points_bounds:
                guids_bounds = compas_rhino.draw_points(points_bounds, layer="TNO::Shape::Bounds", clear=False, redraw=False)
                guids = guids + guids_bounds
                vertices = vertices + vertices_bounds
                self.guid_vertex = zip(guids, vertices)

        if self.settings['show.cracks']:
            points_cracks = []
            vertices_cracks = []
            for key in self.diagram.vertices():
                lb = self.diagram.vertex_attribute(key, 'lb')
                ub = self.diagram.vertex_attribute(key, 'ub')
                x, y, z = self.diagram.vertex_coordinates(key)
                if ub is not None and lb is not None:
                    if abs(ub - z) < 10e-4:
                        vertices_cracks.append(key)
                        points_cracks.append({
                            'pos': [x, y, ub],
                            'color': self.settings['color.vertices:upper_bound']
                        })
                    elif abs(lb - z) < 10e-4:
                        vertices_cracks.append(key)
                        points_cracks.append({
                            'pos': [x, y, lb],
                            'color': self.settings['color.vertices:lower_bound']
                        })
                    else:
                        pass

            if points_cracks:
                guids_cracks = compas_rhino.draw_points(points_cracks, layer="TNO::Shape::Cracks", clear=False, redraw=False)
                guids = guids + guids_cracks
                vertices = vertices + vertices_cracks
                self.guid_vertex = zip(guids, vertices)

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

            self.guid_edgelabel = guid_edgelabel

        # force pipes
        if self.settings['show.forcepipes']:
            guids = self.artist.draw_forcepipes(
                color_compression=self.settings['color.compression'],
                color_tension=self.settings['color.tension'],
                scale=self.settings['scale.forces'],
                tol=self.settings['tol.forces'])

            self.guid_force = zip(guids, edges)

        self.redraw()
