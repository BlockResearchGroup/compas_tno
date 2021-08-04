from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from functools import partial
from compas.utilities import color_to_colordict
from compas_rhino.artists import MeshArtist


colordict = partial(color_to_colordict, colorformat='rgb', normalize=False)


__all__ = ['DiagramArtist']


class DiagramArtist(MeshArtist):
    """Base artist for diagrams in TNO.

    Attributes
    ----------
    diagram : :class:`compas_tno.diagrams.Diagram`
        The diagram associated with the artist.

    """

    @property
    def diagram(self):
        """The diagram assigned to the artist."""
        return self.mesh

    @diagram.setter
    def diagram(self, diagram):
        self.mesh = diagram
