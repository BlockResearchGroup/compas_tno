from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas_tno.shapes.crossvault import cross_vault_highfields
from compas_tno.shapes.pavillionvault import pavillion_vault_highfields
from compas_tno.shapes.dome import set_dome_heighfield

from numpy import array

from scipy import interpolate

__all__ = ['Shape']

class Shape(object):

    """The ``Shape`` class deals with the geometry of a masonry vault, where the definition of Extrados, Intrados and Middle surfaces are of interest.

    Notes
    -----
    A ``Shape`` has the following constructor functions

    *   ``from_library`` : Construct the shape from a dictionary with instructions, the library supports the creation of parametric arches, domes, and vaults.
    *   ``from_rhinomesh`` : Construct Extrados, Intrados and Middle surfaces from RhinoMeshes.
    *   ``from_rhinosurface`` : Construct Extrados, Intrados and Middle surfaces using the U and V isolines.

    A ``Shape`` has the following elements:

    *   ``data``

        *   ``type``  : The type of the vault to be constructed.
        *   ``xy_span``  : Planar range of the structure.
        *   ``discretisation``  : Density of the highfield.
        *   ``thickness`` : Mean thickness of the vault.
        *    ``t`` : The numerical for the negative bounds on the height of the fixed nodes.

    *   ``discretisation``

        *   ``array``    : (2 x n) array with the coordinate of the n points with intrados, extrados and middle evaluated.

    *   ``intrados``

        *   ``Mesh``    : Mesh representing intrados.

    *   ``extrados``

        *   ``Mesh``    : Mesh representing extrados.

    *   ``middle``

        *   ``Mesh``    : Mesh representing middle surface

    """

    __module__ = 'compas_tna.shapes'

    def __init__(self):
        super(Shape, self).__init__()
        self.data = {
            'type'      : None,
            'thk'       : 1.0,
            'discretisation'   : [100, 100],
            'xy_span'   : [0.0, 10.0],
            't'         : 10.0,
        }
        self.intrados = None
        self.extrados = None
        self.middle = None
        self.total_selfweight = 0.0
        self.area = 0.0
        self.ro = 20.0

    @classmethod
    def from_library(cls, data):
        """Construct a FormDiagram from a library.

        Parameters
        ----------
        data : dictionary
            Dictionary with the options to create the vault.

        *   ``type``  : The type of the vault to be constructed. Supported: arch, dome, crossvault, pointedvault, pavillionvault.
        *   ``xy_span``  : Planar range of the structure.
        *   ``discretisation``  : Density of the highfield.
        *   ``thickness`` : Mean thickness of the vault.

        Returns
        -------
        Shape
            A Shape object.

        Examples
        --------
        .. code-block:: python

            import compas
            from compas.files import OBJ
            from compas_tna.diagrams import FormDiagram

            # W.I.P.

        """
        shape = cls()

        typevault = data['type']
        thk = data['thk']
        discretisation = data['discretisation']
        t = data['t']

        if typevault == 'crossvault':
            xy_span = data['xy_span']
            intrados, extrados, middle = cross_vault_highfields(xy_span, thk=thk, discretisation=discretisation, t=t)
        elif typevault == 'pavillionvault':
            xy_span = data['xy_span']
            intrados, extrados, middle = pavillion_vault_highfields(xy_span, thk=thk, discretisation=discretisation, t=t)
        elif typevault == 'dome':
            center = data['center']
            radius = data['radius']
            intrados, extrados, middle = set_dome_heighfield(center, radius=radius, thk=thk, discretisation=discretisation, t=t)

        shape.data = data
        shape.intrados = intrados
        shape.extrados = extrados
        shape.middle = middle

        shape.area = middle.area()
        shape.volume = shape.compute_volume()
        shape.total_selfweight = shape.compute_selfweight()

        return shape


    @classmethod
    def from_rhinosurface(cls):
        ''' Work in progress'''
        shape = cls()

        return shape

    @classmethod
    def from_rhinomesh(cls, data):
        ''' Work in progress'''
        shape = cls()

        return shape

    def get_ub(self, x, y):
        """Get the height of the extrados in the point."""

        vertices = array(self.extrados.get_vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_lb(self, x, y):
        """Get the height of the intrados in the point."""

        vertices = array(self.intrados.get_vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_middle(self, x, y):
        """Get the height of the target/middle surface in the point."""

        vertices = array(self.middle.get_vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_target(self, x, y):
        """Get the height of the target/middle surface in the point."""

        z = get_middle(x,y)

        return z

    def compute_selfweight(self):
        """Compute and returns the selfweight of the structure based on the area and thickness in the data."""

        middle = self.middle
        thk = self.data['thk']
        area = middle.area()
        ro = self.ro
        total_selfweight = thk * area * ro

        return total_selfweight

    def compute_volume(self):
        """Compute and returns the vollume of the structure based on the area and thickness in the data."""

        middle = self.middle
        thk = self.data['thk']
        area = middle.area()
        volume = thk * area

        return volume

