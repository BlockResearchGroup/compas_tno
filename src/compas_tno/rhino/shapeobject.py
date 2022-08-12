from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino
from compas_ui.objects import Object
from compas.geometry import Point


__all__ = ['ShapeObject']


class ShapeObject(Object):
    """A form object represents a Shape in the Rhino model space.
    """

    SETTINGS = {
        'show.intrados': True,
        'show.extrados': True,
        'show.middle': False,

        'color.intrados': (0, 0, 0),
        'color.extrados': (0, 0, 0),
        'color.middle': (0, 0, 0),
    }

    def __init__(self, shape, scene=None, name=None, layer=None, visible=True, settings=None):
        super(ShapeObject, self).__init__(shape, scene, name, layer, visible)
        self.settings.update(ShapeObject.SETTINGS)
        self._shape = shape
        self._guids = []
        self._guid_intrados = {}
        self._guid_extrados = {}
        self._guid_middle = {}
        self._location = None
        self._scale = None
        self._rotation = None
        self._anchor = None
        self.layer = layer
        # TODO: add more guid labeles that are relevant to the shape object
        self.settings.update(type(self).SETTINGS)
        if settings:
            self.settings.update(settings)

    @property
    def shape(self):
        """The shape associated with the object."""
        return self._shape

    @shape.setter
    def shape(self, shape):
        self._shape = shape

    @property
    def location(self):
        """:class:`compas.geometry.Point`:
        The location of the object.
        Default is the origin of the world coordinate system.
        The object transformation is applied relative to this location.
        Setting this location will make a copy of the provided point object.
        Moving the original point will thus not affect the object's location.
        """
        if not self._location:
            self._location = Point(0, 0, 0)
        return self._location

    @location.setter
    def location(self, location):
        self._location = Point(*location)

    @property
    def scale(self):
        """float:
        A uniform scaling factor for the object in the scene.
        The scale is applied relative to the location of the object in the scene.
        """
        if not self._scale:
            self._scale = 1.0
        return self._scale

    @scale.setter
    def scale(self, scale):
        self._scale = scale

    @property
    def rotation(self):
        """list of float:
        The rotation angles around the 3 axis of the coordinate system
        with the origin placed at the location of the object in the scene.
        """
        if not self._rotation:
            self._rotation = [0, 0, 0]
        return self._rotation

    @rotation.setter
    def rotation(self, rotation):
        self._rotation = rotation

    @property
    def anchor(self):
        """The vertex of the mesh that is anchored to the location of the object."""
        return self._anchor

    @anchor.setter
    def anchor(self, vertex):
        if self.shape.middle.has_vertex(vertex):
            self._anchor = vertex

    # @property  # See if this property makes sense or it needs to be vertex_xyz_intrados ...
    # def vertex_xyz(self):
    #     """dict : The view coordinates of the mesh object."""
    #     origin = Point(0, 0, 0)
    #     if self.anchor is not None:
    #         xyz = self.mesh.vertex_attributes(self.anchor, 'xyz')
    #         point = Point(* xyz)
    #         T1 = Translation.from_vector(origin - point)
    #         S = Scale.from_factors([self.scale] * 3)
    #         R = Rotation.from_euler_angles(self.rotation)
    #         T2 = Translation.from_vector(self.location)
    #         X = T2 * R * S * T1
    #     else:
    #         S = Scale.from_factors([self.scale] * 3)
    #         R = Rotation.from_euler_angles(self.rotation)
    #         T = Translation.from_vector(self.location)
    #         X = T * R * S
    #     intrados_transformed = self.intrados.transformed(X)
    #     extrados_transformed = self.extrados.transformed(X)
    #     middle_transformed = self.middle.transformed(X)
    #     vertex_xyz = {vertex: mesh.vertex_attributes(vertex, 'xyz') for vertex in mesh.vertices()}
    #     return vertex_xyz

    @property
    def guid_intrados(self):
        """dict: Map between Rhino object GUIDs and intrados identifiers."""
        return self._guid_intrados

    @guid_intrados.setter
    def guid_intrados(self, values):
        self._guid_intrados = dict(values)

    @property
    def guid_extrados(self):
        """dict: Map between Rhino object GUIDs and intrados identifiers."""
        return self._guid_extrados

    @guid_extrados.setter
    def guid_extrados(self, values):
        self._guid_extrados = dict(values)

    @property
    def guid_middle(self):
        """dict: Map between Rhino object GUIDs and intrados identifiers."""
        return self._guid_middle

    @guid_middle.setter
    def guid_middle(self, values):
        self._guid_middle = dict(values)

    @property
    def guids(self):
        """list: The GUIDs of all Rhino objects created by this artist."""
        guids = self._guids
        guids += list(self.guid_intrados.keys())
        guids += list(self.guid_extrados.keys())
        guids += list(self.guid_middle.keys())
        return guids

    def clear(self):
        """Clear all Rhino objects associated with this object.
        """
        compas_rhino.delete_objects(self.guids, purge=True)
        self._guids = []
        self._guid_intrados = {}
        self._guid_extrados = {}
        self._guid_middle = {}

    def draw(self):
        """Draw the shape object.

        The visible components, display properties and visual style of the shape
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

        self.artist.shape = self.shape

        if not self.visible:
            return

        # intrados
        if self.settings['show.intrados']:
            guid = self.artist.draw_intrados(color=self.settings['color.intrados'], layer=self.layer + '::Intrados')
            self.guid_intrados = zip(guid, 'intrados')

        # extrados
        if self.settings['show.extrados']:
            guid = self.artist.draw_extrados(color=self.settings['color.extrados'], layer=self.layer + '::Extrados')
            self.guid_extrados = zip(guid, 'extrados')

        # middle
        if self.settings['show.middle']:
            guid = self.artist.draw_middle(color=self.settings['color.middle'], layer=self.layer + '::Middle')
            self.guid_middle = zip(guid, 'middle')

        self.redraw()
