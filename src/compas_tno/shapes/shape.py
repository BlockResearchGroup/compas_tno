from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from copy import deepcopy
from compas.geometry import subtract_vectors
from compas.geometry import length_vector
from compas.geometry import cross_vectors
from compas.geometry import normalize_vector

from compas_tno.shapes import MeshDos

from compas.datastructures import Datastructure


__all__ = ['Shape']


class Shape(Datastructure):

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

    # __module__ = 'compas_tno.shapes'

    def __init__(self):
        super(Shape, self).__init__()
        self.datashape = {
            'type': None,
            'thk': 1.0,
            'discretisation': [20, 20],
            'xy_span': [0.0, 10.0],
            't': 10.0,
        }
        self.intrados = None
        self.extrados = None
        self.middle = None
        self.name = "shape"
        self.total_selfweight = 0.0
        self.area = 0.0
        self.ro = 20.0
        self.fill = False
        self.fill_ro = 20.0

    @property
    def data(self):
        """dict : A data dict representing the shape data structure for serialization.
        """
        dataintrados = None
        dataextrados = None
        datamiddle = None

        if self.intrados:
            dataintrados = self.intrados.to_data()
        if self.extrados:
            dataextrados = self.extrados.to_data()
        if self.middle:
            datamiddle = self.middle.to_data()

        data = {
            'datashape': self.datashape,
            'intrados': dataintrados,
            'extrados': dataextrados,
            'middle': datamiddle,
            'name': self.name,
            'total_selfweight': self.total_selfweight,
            'area': self.area,
            'ro': self.ro,
            'fill': self.fill,
            'fill_ro': self.fill_ro
        }
        return data

    @data.setter
    def data(self, data):
        if 'data' in data:
            data = data['data']
        self.datashape = data.get('datashape') or {}

        self.name = data.get('name') or "shape"
        self.total_selfweight = data.get('total_selfweight') or 0.0
        self.area = data.get('area') or 0.0
        self.ro = data.get('ro') or 20.0
        self.fill = data.get('fill') or False
        self.fill_ro = data.get('ro') or 20.0

        dataintrados = data.get('intrados')
        dataextrados = data.get('extrados')
        datamiddle = data.get('middle')

        self.intrados = MeshDos.from_data(dataintrados) if dataintrados else None
        self.extrados = MeshDos.from_data(dataextrados) if dataintrados else None
        self.middle = MeshDos.from_data(datamiddle) if dataintrados else None

    @classmethod
    def from_library(cls, data):
        """Construct a Shape from a library.

        Parameters
        ----------
        data : dictionary
            Dictionary with the settings to create the vault.

        Returns
        -------
        Shape
            A Shape object.

        """

        from compas_tno.shapes import cross_vault_highfields
        from compas_tno.shapes import pavillion_vault_highfields
        from compas_tno.shapes import set_dome_heighfield
        from compas_tno.shapes import set_dome_with_spr
        from compas_tno.shapes import set_dome_polar_coord
        from compas_tno.shapes import arch_shape
        from compas_tno.shapes import pointed_arch_shape
        from compas_tno.shapes import pointed_vault_heightfields
        from compas_tno.shapes import domical_vault
        from compas_tno.shapes import parabolic_shell_highfields

        shape = cls()

        typevault = data['type']
        thk = data['thk']
        discretisation = data['discretisation']
        t = data.get('t', 0.0)
        expanded = data.get('expanded', False)

        if typevault == 'crossvault':
            xy_span = data['xy_span']
            intrados, extrados, middle = cross_vault_highfields(xy_span, thk=thk, discretisation=discretisation, t=t, expanded=expanded)
        elif typevault == 'pavillionvault':
            xy_span = data['xy_span']
            intrados, extrados, middle = pavillion_vault_highfields(xy_span, thk=thk, discretisation=discretisation, t=t)
        elif typevault == 'dome':
            center = data['center']
            radius = data['radius']
            intrados, extrados, middle = set_dome_heighfield(center, radius=radius, thk=thk, discretisation=discretisation, t=t, expanded=expanded)
        elif typevault == 'dome_polar':
            center = data['center']
            radius = data['radius']
            intrados, extrados, middle = set_dome_polar_coord(center, radius=radius, thk=thk, discretisation=discretisation, t=t)
        elif typevault == 'dome_spr':
            center = data['center']
            radius = data['radius']
            theta = data['theta']
            intrados, extrados, middle = set_dome_with_spr(center, radius=radius, thk=thk, theta=theta, discretisation=discretisation, t=t)
        elif typevault == 'parabolic_shell':
            xy_span = data['xy_span']
            intrados, extrados, middle = parabolic_shell_highfields(xy_span, thk=thk, discretisation=discretisation, t=t)
        elif typevault == 'arch':
            H = data['H']
            L = data['L']
            b = data['b']
            x0 = data['x0']
            intrados, extrados, middle = arch_shape(H=H, L=L, x0=x0, thk=thk, total_nodes=discretisation, b=b, t=t)
        elif typevault == 'pointed_arch':
            hc = data['hc']
            L = data['L']
            b = data['b']
            x0 = data['x0']
            intrados, extrados, middle = pointed_arch_shape(hc=hc, L=L, x0=x0, thk=thk, total_nodes=discretisation, b=b, t=t)
        elif typevault == 'pointed_crossvault':
            xy_span = data['xy_span']
            hc = data['hc']
            hm = data.get('hm', None)
            he = data.get('he', None)
            intrados, extrados, middle = pointed_vault_heightfields(xy_span=xy_span, discretisation=discretisation, t=t, hc=hc, he=he, hm=hm, thk=thk)
        elif typevault == 'domicalvault':
            xy_span = data['xy_span']
            center = [sum(xy_span[0])/2, sum(xy_span[1])/2]
            radius = data.get('radius', None)
            intrados, extrados, middle = domical_vault(xy_span=xy_span, thk=thk, radius=radius, center=center, t=t, discretisation=discretisation)

        shape.datashape = data
        shape.intrados = intrados
        shape.extrados = extrados
        shape.middle = middle

        shape.area = middle.area()
        shape.volume = shape.compute_volume()
        shape.total_selfweight = shape.compute_selfweight()

        return shape

    @classmethod
    def from_library_proxy(cls, data):
        """Construct a Shape from a library using a proxy implementation (compatible with Rhinoceros).

        Parameters
        ----------
        data : dictionary
            Dictionary with the setttings to create the vault.

        Returns
        -------
        Shape
            A Shape object.

        """

        from compas.rpc import Proxy

        proxy = Proxy()
        proxy.package = 'compas_tno.shapes'

        typevault = data['type']
        thk = data['thk']
        discretisation = data['discretisation']
        t = data.get('t', 0.0)
        expanded = data.get('expanded', False)

        if typevault == 'crossvault':
            xy_span = data['xy_span']
            intra_data, extra_data, middle_data = proxy.cross_vault_highfields_proxy(xy_span, thk=thk, discretisation=discretisation, t=t, expanded=expanded)
        elif typevault == 'pavillionvault':
            xy_span = data['xy_span']
            intra_data, extra_data, middle_data = proxy.pavillion_vault_highfields_proxy(xy_span, thk=thk, discretisation=discretisation, t=t, expanded=expanded)
        elif typevault == 'pointed_crossvault':
            xy_span = data['xy_span']
            hc = data['hc']
            hm = data.get('hm', None)
            he = data.get('he', None)
            intra_data, extra_data, middle_data = proxy.pointed_vault_heightfields_proxy(xy_span=xy_span, thk=thk, discretisation=discretisation, t=t, hc=hc, he=he, hm=hm)
        elif typevault == 'dome':
            center = data['center']
            radius = data['radius']
            intrados, extrados, middle = proxy.set_dome_heighfield_proxy(center, radius=radius, thk=thk, discretisation=discretisation, t=t, expanded=expanded)

        intrados = MeshDos.from_data(intra_data)
        extrados = MeshDos.from_data(extra_data)
        middle = MeshDos.from_data(middle_data)

        shape = cls()
        shape.datashape = data
        shape.intrados = intrados
        shape.extrados = extrados
        shape.middle = middle

        shape.area = middle.area()
        shape.volume = shape.compute_volume()
        shape.total_selfweight = shape.compute_selfweight()

        return shape

    @classmethod
    def from_meshes(cls, intrados, extrados, middle=None, treat_creases=False, data=None):
        """Construct a Shape from meshes for intrados, extrados, and middle.

        Parameters
        ----------
        intrados : mesh
            Mesh for intrados.
        extrados : mesh
            Mesh for extrados.
        middle : mesh, optional
            Mesh for middle.
            The default value is ``None``, in a case in which no middle surface is assigned
        data : dict, optional
            Dictionary with the data about the structure.
            The default value is ``None``, no data is assigned.

        Returns
        -------
        Shape
            A Shape object.

        """
        shape = cls()

        shape.intrados = MeshDos.from_mesh(intrados)
        shape.extrados = MeshDos.from_mesh(extrados)
        if middle:
            shape.middle = MeshDos.from_mesh(middle)
        else:
            shape.interpolate_middle_from_ub_lb()
            shape.middle = MeshDos.from_mesh(shape.intrados.copy())
            for key in shape.intrados.vertices():
                ub = shape.extrados.vertex_attribute(key, 'z')
                lb = shape.intrados.vertex_attribute(key, 'z')
                shape.middle.vertex_attribute(key, 'z', (ub + lb)/2)

        if data:
            shape.datashape = data

        return shape

    @classmethod
    def from_pointcloud(cls, intrados_pts, extrados_pts, middle=None, data={'type': 'general', 't': 0.0}):
        """Construct a Shape from a pointcloud.

        Parameters
        ----------
        intrados_pts : list
            List of points collected in the intrados.
        extrados : mesh
            List of points collected in the extrados.
        middle : mesh (None)
            Mesh for middle.
        data : dict (None)
            Dictionary with the data in required.

        Returns
        -------
        Shape
            A Shape object.

        """

        intrados_mesh = MeshDos.from_points_delaunay(intrados_pts)
        extrados_mesh = MeshDos.from_points_delaunay(extrados_pts)
        shape = cls().from_meshes(intrados_mesh, extrados_mesh, middle=middle, data=data)

        return shape

    @classmethod
    def from_middle(cls, middle, thk=0.50, treat_creases=False, printout=False, data={'type': 'general', 't': 0.0, 'xy_span': [[0.0, 10.0], [0.0, 10.0]]}):
        """Construct a Shape from a middle mesh and a thickness (offset happening to both sides).

        Parameters
        ----------
        middle : Mesh
            Mesh with middle. Topology will be considered.
        data : dict (None)
            Dictionary with the data in required.

        Returns
        -------
        Shape
            A Shape object.

        """
        middle = MeshDos.from_mesh(middle)
        if treat_creases:
            middle.identify_creases_at_diagonals(xy_span=data['xy_span'])
            middle.store_normals(correct_creases=True)
        else:
            middle.store_normals()
        extrados_mesh = middle.offset_mesh(thk/2, direction='up')
        intrados_mesh = middle.offset_mesh(thk/2, direction='down')
        data['thk'] = thk
        shape = cls().from_meshes(intrados_mesh, extrados_mesh, middle=middle, data=data)

        return shape

    # @classmethod
    # def from_middle_pointcloud(cls, middle_pts, topology=None, thk=0.50, treat_creases=False, printout=False, data={'type': 'general', 't': 0.0, 'xy_span': [[0.0, 10.0], [0.0, 10.0]]}):
    #     """Construct a Shape from a pointcloud.

    #     Parameters
    #     ----------
    #     middle_pts : list
    #         List of points collected in the intrados.
    #     topology : mesh (None)
    #         If a mesh is given the middle, intrados and extrados meshes will repeat the topology.
    #     data : dict (None)
    #         Dictionary with the data in required.

    #     Returns
    #     -------
    #     Shape
    #         A Shape object.

    #     """

    #     if not topology:
    #         middle = MeshDos.from_points_delaunay(middle_pts)
    #     else:
    #         middle = MeshDos.from_topology_and_pointcloud(topology, array(middle_pts))
    #     if treat_creases:
    #         middle.identify_creases_at_diagonals(xy_span=data['xy_span'])
    #         middle.store_normals(correct_creases=True)
    #     else:
    #         middle.store_normals()
    #     extrados_mesh = middle.offset_mesh(thk/2, direction='up')
    #     intrados_mesh = middle.offset_mesh(thk/2, direction='down')
    #     data['thk'] = thk
    #     shape = cls().from_meshes(intrados_mesh, extrados_mesh, middle=middle, data=data)

    #     return shape

    # @classmethod
    # def from_pointcloud_and_formdiagram(cls, form, intrados_pts, extrados_pts, middle=None, data={'type': 'general', 't': 0.0}):
    #     """Construct a Shape from a pointcloud and a formdiagram that will have its topology copied.

    #     Parameters
    #     ----------
    #     form: FormDiagram
    #         Form Diagram with the topology to be used.
    #     intrados_pts : list
    #         List of points collected in the intrados.
    #     extrados : mesh
    #         List of points collected in the extrados.
    #     middle : mesh (None)
    #         Mesh for middle.
    #     data : dict (None)
    #         Dictionary with the data in required.

    #     Returns
    #     -------
    #     Shape
    #         A Shape object.

    #     """

    #     intrados_mesh = MeshDos.from_topology_and_pointcloud(form, array(intrados_pts))
    #     extrados_mesh = MeshDos.from_topology_and_pointcloud(form, array(extrados_pts))
    #     shape = cls().from_meshes(intrados_mesh, extrados_mesh, middle=middle, data=data)

    #     return shape

    @classmethod
    def from_meshes_and_formdiagram(cls, form, intrados, extrados, middle=None, data={'type': 'general', 't': 0.0}):
        """Construct a Shape from a pair of meshes and a formdiagram that will have its topology copied and normals copied.

        Parameters
        ----------
        form: FormDiagram
            Form Diagram with the topology to be used.
        intrados : mesh
            Mesh for intrados.
        extrados : mesh
            Mesh for extrados.
        middle : mesh (None)
            Mesh for middle.
        data : dict (None)
            Dictionary with the data in required.

        Returns
        -------
        Shape
            A Shape object.

        """

        intrados_mesh = MeshDos.from_topology_and_mesh(form, intrados, keep_normals=True)
        extrados_mesh = MeshDos.from_topology_and_mesh(form, extrados, keep_normals=True)
        shape = cls().from_meshes(intrados_mesh, extrados_mesh, middle=middle, data=data)

        return shape

    @classmethod
    def from_formdiagram_and_attributes(cls, form, data=None):
        """Construct a Shape from the form diagram and its attributes 'ub' and 'lb'.

        Parameters
        ----------
        form: FormDiagram
            Form Diagram with the topology to be used.
        intrados : mesh
            Mesh for intrados.
        extrados : mesh
            Mesh for extrados.
        middle : mesh (None)
            Mesh for middle.
        data : dict (None)
            Dictionary with the data if required.

        Returns
        -------
        Shape
            A Shape object.

        """

        vertices, faces = form.to_vertices_and_faces()
        intra = MeshDos.from_vertices_and_faces(vertices, faces)
        extra = MeshDos.from_vertices_and_faces(vertices, faces)
        middle = MeshDos.from_vertices_and_faces(vertices, faces)

        for key in form.vertices():
            intra.vertex_attribute(key, 'z', form.vertex_attribute(key, 'lb'))
            extra.vertex_attribute(key, 'z', form.vertex_attribute(key, 'ub'))
            middle.vertex_attribute(key, 'z', form.vertex_attribute(key, 'target'))

        shape = cls().from_meshes(intra, extra, middle=middle, data=data)

        return shape

    @classmethod
    def from_assembly(self):
        NotImplementedError

    def add_damage_from_meshes(self, intrados, extrados):

        self.intrados_damage = intrados
        self.extrados_damage = extrados

        return

    def interpolate_middle_from_ub_lb(self, intrados=None, extrados=None):

        if not intrados:
            intrados = self.intrados

        if not extrados:
            extrados = self.extrados

        self.middle = MeshDos.from_mesh(self.intrados.copy())

        for key in self.intrados.vertices():
            ub = self.extrados.vertex_attribute(key, 'z')
            lb = self.intrados.vertex_attribute(key, 'z')
            self.middle.vertex_attribute(key, 'z', (ub + lb)/2)

        return

    def store_normals(self, mark_fixed_LB=True, plot=False):

        intrados = self.intrados
        extrados = self.extrados

        intrados.store_normals()
        extrados.store_normals()

        intrados.vertices_attribute('is_outside', False)
        extrados.vertices_attribute('is_outside', False)

        if mark_fixed_LB:  # Look for the vertices with LB == 't'
            try:
                t = self.datashape['t']
            except BaseException:
                print('No t is assigned to vertices in intrados.')
                t = None
            if t is not None:
                for key in intrados.vertices():
                    if abs(intrados.vertex_attribute(key, 'z') - t) < 10e-3:
                        intrados.vertex_attribute(key, 'is_outside', True)
                    if abs(extrados.vertex_attribute(key, 'z') - t) < 10e-3:
                        extrados.vertex_attribute(key, 'is_outside', True)

        if plot:
            intrados.plot_normals()
            extrados.plot_normals()

        return

    def analytical_normals(self, assume_shape=None, mark_fixed_LB=True):

        if assume_shape:
            data = assume_shape
        else:
            data = self.datashape

        if data['type'] == 'dome':
            xc, yc = data['center']
        else:
            return Exception

        intrados = self.intrados
        extrados = self.extrados

        for key in intrados.vertices():
            x, y, zlb = intrados.vertex_coordinates(key)
            _, _, zub = extrados.vertex_coordinates(key)
            nlb = normalize_vector([(x - xc), (y - yc), zlb])
            nub = normalize_vector([(x - xc), (y - yc), zub])
            intrados.vertex_attribute(key, 'n', nlb)
            extrados.vertex_attribute(key, 'n', nub)

        intrados.vertices_attribute('is_outside', False)
        extrados.vertices_attribute('is_outside', False)

        if mark_fixed_LB:  # Look for the vertices with LB == 't'
            try:
                t = self.datashape['t']
            except BaseException:
                print('No t is assigned to vertices in intrados.')
                t = None
            if t is not None:
                for key in intrados.vertices():
                    if abs(intrados.vertex_attribute(key, 'z') - t) < 10e-3:
                        intrados.vertex_attribute(key, 'is_outside', True)
                    if abs(extrados.vertex_attribute(key, 'z') - t) < 10e-3:
                        extrados.vertex_attribute(key, 'is_outside', True)

        return

    def compute_selfweight(self):
        """Compute and returns the total selfweight of the structure based on the area and thickness in the data.

        Returns
        -------
        swt : float
            The selfweight.
        """

        middle = self.middle
        ro = self.ro
        try:
            thk = self.datashape['thk']
            area = middle.area()
            total_selfweight = thk * area * ro
        except BaseException:
            intrados = self.intrados
            extrados = self.extrados
            total_selfweight = 0
            for key in intrados.vertices():
                h = extrados.vertex_attribute(key, 'z') - intrados.vertex_attribute(key, 'z')
                vol = h * middle.vertex_projected_area(key)
                total_selfweight += vol * ro

        return total_selfweight

    def compute_volume(self):
        """Compute and returns the vollume of the structure based on the area and thickness in the data.

        Returns
        -------
        vol : float
            The volume.
        """

        middle = self.middle
        thk = self.datashape['thk']
        area = middle.area()
        volume = thk * area

        return volume

    def add_fill_with_height(self, height, fill_ro=20.0):
        """Compute the weight of the fill applied based on the height.

        Returns
        -------
        height : float
            The height until which infill should be observed.
        fill_ro : float (20)
            The density of the infill.
        """

        self.extrados_fill = deepcopy(self.extrados)
        volume = 0.
        proj_area_total = 0.

        for key in self.extrados_fill.vertices():
            zi = self.extrados_fill.vertex_attribute(key, 'z')
            if zi < height:
                self.extrados_fill.vertex_attribute(key, 'z', height)
                # this should become a method
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
                # this should become a method
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
