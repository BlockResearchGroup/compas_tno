from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from compas_tno.shapes.crossvault import cross_vault_highfields
from compas_tno.shapes.pavillionvault import pavillion_vault_highfields
from compas_tno.shapes.dome import set_dome_heighfield
from compas_tno.shapes.dome import set_dome_with_spr
from compas_tno.shapes.dome import set_dome_polar_coord
from compas_tno.shapes.circular_arch import arch_shape
from compas_tno.shapes.pointed_crossvault import pointed_vault_heightfields

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
        self.fill = False
        self.fill_ro = 20.0

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
        elif typevault == 'dome_polar':
            center = data['center']
            radius = data['radius']
            intrados, extrados, middle = set_dome_polar_coord(center, radius=radius, thk=thk, discretisation=discretisation, t=t)
        elif typevault == 'dome_spr':
            center = data['center']
            radius = data['radius']
            theta = data['theta']
            intrados, extrados, middle = set_dome_with_spr(center, radius=radius, thk=thk, theta=theta, discretisation=discretisation, t=t)
        elif typevault == 'arch':
            H = data['H']
            L = data['L']
            b = data['b']
            x0 = data['x0']
            intrados, extrados, middle = arch_shape(H=H, L=L, x0=x0, thk=thk, total_nodes=discretisation, b=b, t=t)
        elif typevault == 'pointed_crossvault':
            xy_span = data['xy_span']
            hc = data['hc']
            hm = data['hm']
            he = data['he']
            intrados, extrados, middle = pointed_vault_heightfields(xy_span=xy_span, discretisation=discretisation, hc=hc, he=he, hm=hm, thk=thk)
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
        NotImplementedError

    @classmethod
    def from_rhinomesh(cls, data):
        ''' Work in progress'''
        NotImplementedError

    @classmethod
    def from_assembly():
        NotImplementedError

    def get_ub(self, x, y):
        """Get the height of the extrados in the point."""

        vertices = array(self.extrados.vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_ub_fill(self, x, y):
        """Get the height of the fill in the point."""

        vertices = array(self.extrados_fill.vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_lb(self, x, y):
        """Get the height of the intrados in the point."""

        vertices = array(self.intrados.vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_middle(self, x, y):
        """Get the height of the target/middle surface in the point."""

        vertices = array(self.middle.vertices_attributes('xyz'))
        z = float(interpolate.griddata(vertices[:,:2], vertices[:,2], [x,y]))

        return z

    def get_target(self, x, y):
        """Get the height of the target/middle surface in the point."""

        z = self.get_middle(x,y)

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

    def add_fill_with_height(self, height, fill_ro=20.0):

        from copy import deepcopy
        from compas.geometry import subtract_vectors
        from compas.geometry import length_vector
        from compas.geometry import cross_vectors

        self.extrados_fill = deepcopy(self.extrados)
        volume = 0.
        proj_area_total = 0.

        for key in self.extrados_fill.vertices():
            zi = self.extrados_fill.vertex_attribute(key, 'z')
            if zi < height:
                self.extrados_fill.vertex_attribute(key, 'z', height)
                ##### this should become a method
                area = 0.
                p0 = self.extrados_fill.vertex_coordinates(key)
                p0[2] = 0
                for nbr in self.extrados_fill.halfedge[key]:
                    p1 = self.extrados_fill.vertex_coordinates(nbr)
                    p1[2] = 0
                    v1 = subtract_vectors(p1, p0)
                    fkey = self.extrados_fill.halfedge[key][nbr]
                    if fkey is not None:
                        p2 = self.extrados_fill.face_centroid(fkey)
                        p2[2] = 0
                        v2 = subtract_vectors(p2, p0)
                        area += length_vector(cross_vectors(v1, v2))
                    fkey = self.extrados_fill.halfedge[nbr][key]
                    if fkey is not None:
                        p3 = self.extrados_fill.face_centroid(fkey)
                        p3[2] = 0
                        v3 = subtract_vectors(p3, p0)
                        area += length_vector(cross_vectors(v1, v3))
                proj_area = area * 0.25
                proj_area_total += proj_area
                ##### this should become a method
                vol_i = proj_area*(height - zi)
                volume += vol_i

        self.fill = True
        self.fill_volume = volume
        print('Proj area total of shape', proj_area_total)
        self.fill_ro = fill_ro

        return


    def compute_fill_weight(self):

        return self.fill_volume * self.fill_ro

    def add_fill_with_angle(self, angle):
        pass
