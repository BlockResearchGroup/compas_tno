from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import compas_rhino
from compas_ui.objects import Object
from compas.geometry import Point


__all__ = ['OptimiserObject']


class OptimiserObject(Object):
    """A form object represents a Shape in the Rhino model space.
    """

    SETTINGS = {
        'SolverSelection.library': 'SLSQP',
        'SolverSelection.solver': 'SLSQP',
        'optimisation.constraints': ['funicular', 'envelope'],
        'optimisation.variables': ['q', 'zb'],
        'optimisation.features': ['fixed'],
        'optimisation.objective': 'min',
        'optimisation.gradient': True,
        'optimisation.jacobian': True,
        'optimisation.starting_point': 'current',
        'parameters.qmin': -1e+4,
        'parameters.qmax': + 1e-8,
        'parameters.max_iter': 500,
    }

    def __init__(self, optimiser, scene=None, name=None, layer=None, visible=True, settings=None):
        super(OptimiserObject, self).__init__(optimiser, scene, name, layer, visible)
        self.settings.update(OptimiserObject.SETTINGS)
        self.optimiser = optimiser
        self._guids = []
        self._location = None
        self._scale = None
        self._rotation = None
        self._anchor = None
        self.layer = layer
        self.settings.update(type(self).SETTINGS)
        if settings:
            self.settings.update(settings)

    def update_object_from_optimiser(self):
        """Update the settings of the object based on the optimiser settings.
        """

        if not self.optimiser:
            return

        keys = self.settings.keys()

        for key in keys:
            key_end = key.split('.')[1]
            self.settings[key] = self.optimiser.settings[key_end]

    def update_optimiser_from_object(self):
        """Update the settings of the optimiser based on the object settings.
        """

        if not self.optimiser:
            return

        keys = self.settings.keys()

        for key in keys:
            key_end = key.split('.')[1]
            self.optimiser.settings[key_end] = self.settings[key]

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
        if self.middle.has_vertex(vertex):
            self._anchor = vertex

    @property
    def guids(self):
        """list: The GUIDs of all Rhino objects created by this artist."""
        guids = self._guids
        return guids

    def clear(self):
        """Clear all Rhino objects associated with this object.
        """
        compas_rhino.delete_objects(self.guids, purge=True)
        self._guids = []

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
        self.update_optimiser_from_object()

        self.clear()

        if not self.visible:
            return

        self.redraw()
