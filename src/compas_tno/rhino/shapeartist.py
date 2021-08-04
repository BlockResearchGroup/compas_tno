from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas_rhino.artists import BaseArtist
from compas_rhino.artists import MeshArtist


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
        self.artist_intrados = MeshArtist(shape.intrados)
        self.artist_extrados = MeshArtist(shape.extrados)
        self.artist_middle = MeshArtist(shape.middle)

    def draw_intrados(self, color=(0, 0, 0), displacement=None, layer=None):

        self.artist_intrados.layer = layer

        return self.artist_intrados.draw_mesh(color=color)

    def draw_extrados(self, color=(0, 0, 0), displacement=None, layer=None):

        self.artist_extrados.layer = layer

        return self.artist_extrados.draw_mesh(color=color)

    def draw_middle(self, color=(0, 0, 0), displacement=None, layer=None):

        self.artist_middle.layer = layer

        return self.artist_middle.draw_mesh(color=color)
