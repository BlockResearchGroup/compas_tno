import math  # noqa: I001
from copy import deepcopy
from typing import Optional
from typing import TYPE_CHECKING

from compas.data import Data
from compas.datastructures import Mesh
from compas.geometry import cross_vectors
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import subtract_vectors

from compas_tno.shapes import MeshDos

from compas_tno.shapes import arch_shape
from compas_tno.shapes import arch_shape_polar
from compas_tno.shapes import cross_vault_highfields
from compas_tno.shapes import domical_vault
from compas_tno.shapes import parabolic_shell_highfields
from compas_tno.shapes import pavillion_vault_highfields
from compas_tno.shapes import pointed_arch_shape
from compas_tno.shapes import pointed_vault_heightfields
from compas_tno.shapes import set_dome_heighfield
from compas_tno.shapes import set_dome_polar_coord
from compas_tno.shapes import set_dome_with_spr

from compas_tno.utilities import create_mesh_from_topology_and_basemesh
from compas_tno.utilities import create_mesh_from_topology_and_pointcloud
from compas_tno.utilities import mesh_from_pointcloud
from compas_tno.utilities import interpolate_from_pointcloud

if TYPE_CHECKING:
    from compas_tno.diagrams import FormDiagram


class Shape(Data):
    """The ``Shape`` class deals with the geometry of a masonry vault, where the definition of Extrados, Intrados and Middle surfaces are of interest.

    The class imports the attributes and set ups from ``compas_tna.diagrams.FormDiagram`` and include some
    functionalities useful for the assessment of masonry structures. Most relevant edge and vertex attributes
    are listed below:

    Parameters
    ----------
    intrados : MeshDos
        Mesh representing the intrados of the masonry
    extrados : MeshDos
        Mesh representing thhe extrados of the masonry
    middle : MeshDos
        Mesh representing the midde surface of the masonry
    name : str
        The name of the shape
    total_selfweight : float
        Total weight of the masonry in kN.
    area : float
        Area of the middle surface of the masonry
    ro : float
        The density of the masonry
    fill : bool
        Whether or not a fill should be considered up to a given elevation
    fill_ro : float
        The density of the masonry in the fill

    Attributes
    ----------
    type : str
        The type of the masory. Types include:
        ['crossvault', 'pavillionvault', 'dome', 'dome_polar', 'dome_spr', 'parabolic_shell', 'arch', 'pointed_arch', 'pointed_crossvault']
    thk : float
        The thickness of the masonry
    discretisation : int
        The discretisation level
    xy_span : list
        The limits of the xy box (for rectangular form diagram)
    t : float
        The parameter assumed when no intrados projection is found

    Examples
    --------
    >>> from compas_tno.shapes import Shape
    >>> data = {"type": "crossvault", "thk": 0.5, "discretisation": 10, "xy_span": [[0.0, 10.0], [0.0, 10.0]], "t": 0.0}
    >>> shape = Shape.from_library(data)

    """

    intrados: MeshDos
    extrados: MeshDos
    middle: MeshDos

    def __init__(self, middle: MeshDos, intrados: MeshDos, extrados: MeshDos, fill: Optional[Mesh] = None, **kwargs) -> None:
        super().__init__(**kwargs)

        self.parameters = {}

        self.middle = middle
        self.intrados = intrados
        self.extrados = extrados
        self.fill = fill

        self.ro = 20.0
        self.fill_ro = 14.0

        self._area = 0.0
        self._volume = 0.0
        self._total_selfweight = 0.0

    @property
    def __data__(self) -> dict:
        raise NotImplementedError

    @classmethod
    def __from_data__(cls, data: dict) -> Data:
        return super().__from_data__(data)

    # =============================================================================
    # Properties
    # =============================================================================

    @property
    def area(self):
        if not self._area:
            self._area = self.middle.area()
        return self._area

    @property
    def volume(self):
        if not self._volume:
            self._volume = self.compute_volume()
        return self._volume

    @property
    def total_selfweight(self):
        if not self._total_selfweight:
            self._total_selfweight = self.compute_selfweight()
        return self._total_selfweight

    # =============================================================================
    # Constructors
    # =============================================================================

    @classmethod
    def from_library(cls, data: dict) -> "Shape":
        """Construct a Shape from a library.

        Parameters
        ----------
        data : dict
            Dictionary with the settings to create the vault.

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        typevault = data.pop("type")

        if typevault == "crossvault":
            return cls.create_crossvault(**data)

        if typevault == "pavillionvault":
            return cls.create_pavillionvault(**data)

        if typevault == "dome":
            return cls.create_dome(**data)

        if typevault == "dome_polar":
            return cls.create_dome_polar(**data)

        if typevault == "arch":
            return cls.create_arch(**data)

        if typevault == "arch_polar":
            return cls.create_arch_polar(**data)

        if typevault == "pointed_crossvault":
            return cls.create_pointedcrossvault(**data)

        # if typevault == "pointed_arch":
        #     hc = data["hc"]
        #     L = data["L"]
        #     b = data["b"]
        #     x0 = data["x0"]
        #     intrados, extrados, middle = pointed_arch_shape(hc=hc, L=L, x0=x0, thk=thk, total_nodes=discretisation, b=b, t=t)

        # if typevault == "dome_spr":
        #     center = data["center"]
        #     radius = data["radius"]
        #     theta = data["theta"]
        #     intrados, extrados, middle = set_dome_with_spr(center, radius=radius, thk=thk, theta=theta, discretisation=discretisation, t=t)

        # if typevault == "parabolic_shell":
        #     xy_span = data["xy_span"]
        #     intrados, extrados, middle = parabolic_shell_highfields(xy_span, thk=thk, discretisation=discretisation, t=t)

        # if typevault == "domicalvault":
        #     xy_span = data["xy_span"]
        #     center = [sum(xy_span[0]) / 2, sum(xy_span[1]) / 2]
        #     radius = data.get("radius", None)
        #     intrados, extrados, middle = domical_vault(xy_span=xy_span, thk=thk, radius=radius, center=center, t=t, discretisation=discretisation)

        raise NotImplementedError

    @classmethod
    def from_meshes(cls, intrados: Mesh, extrados: Mesh, middle: Optional[Mesh] = None, fill: Optional[Mesh] = None) -> "Shape":
        """Construct a Shape from meshes for intrados, extrados, and middle.

        Parameters
        ----------
        intrados : :class:`~compas.datastructures.Mesh`
            Mesh for intrados.
        extrados : :class:`~compas.datastructures.Mesh`
            Mesh for extrados.
        middle : :class:`~compas.datastructures.Mesh`, optional
            Mesh for middle.
        fill : :class:`~compas.datastructures.Mesh`, optional
            Mesh for fill.

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        intrados = MeshDos.from_mesh(intrados)
        extrados = MeshDos.from_mesh(extrados)
        if middle:
            middle = MeshDos.from_mesh(middle)
        else:
            middle: MeshDos = intrados.copy()
            for vertex in intrados.vertices():
                ub = extrados.vertex_attribute(vertex, "z")
                lb = intrados.vertex_attribute(vertex, "z")
                middle.vertex_attribute(vertex, "z", (ub + lb) / 2)

        shape = cls(middle, intrados, extrados, fill)
        shape.parameters = {"type": "general"}
        return shape

    @classmethod
    def from_pointcloud(cls, intrados_pts, extrados_pts, middle_pts=None) -> "Shape":
        """Construct a Shape from a pointcloud.

        Parameters
        ----------
        intrados_pts : list
            List of points collected in the intrados.
        extrados_pts : list
            List of points collected in the extrados.
        middle_pts : list, optional
            List of points collected in the middle.
            The default value is ``None``, in which case no middle surface is connstructed

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        extrados = mesh_from_pointcloud(extrados_pts)
        intrados = mesh_from_pointcloud(intrados_pts)
        middle = None
        if middle_pts:
            middle = mesh_from_pointcloud(middle_pts)

        shape = cls.from_meshes(intrados, extrados, middle=middle)
        shape.parameters = {"type": "general"}
        return shape

    @classmethod
    def from_middle(cls, middle: Mesh, thk=0.50, treat_creases=False, xy_span=[[0.0, 10.0], [0.0, 10.0]]) -> "Shape":
        """Construct a Shape from a middle mesh and a thickness (offset happening to compute new extrados).

        Parameters
        ----------
        middle : :class:`~compas.datastructures.Mesh`
            Mesh with middle. Topology will be considered.

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        middle: MeshDos = MeshDos.from_mesh(middle)
        if treat_creases:
            middle.identify_creases_at_diagonals(xy_span=xy_span)
            middle.store_normals(correct_creases=True)
        else:
            middle.store_normals()

        extrados_mesh = middle.offset_mesh(thk / 2, direction="up")
        intrados_mesh = middle.offset_mesh(thk / 2, direction="down")

        shape = cls.from_meshes(intrados_mesh, extrados_mesh, middle=middle)
        shape.parameters = {"type": "general", "thk": thk, "xy_span": xy_span, "treat_creases": treat_creases}
        return shape

    @classmethod
    def from_intrados(cls, intrados: Mesh, thk=0.50, treat_creases=False, xy_span=[[0.0, 10.0], [0.0, 10.0]]) -> "Shape":
        """Construct a Shape from an intrados mesh and a thickness (offset happening to both sides).

        Parameters
        ----------
        intrados : :class:`~compas.datastructures.Mesh`
            Mesh with middle. Topology will be considered.

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        intrados = MeshDos.from_mesh(intrados)
        if treat_creases:
            intrados.identify_creases_at_diagonals(xy_span=xy_span)
            intrados.store_normals(correct_creases=True)
        else:
            intrados.store_normals()
        extrados = intrados.offset_mesh(thk, direction="up")

        # bbox = intrados.bounding_box_xy()
        # data['xy_span'] = TAKE THIS FROM BBOX

        shape = cls.from_meshes(intrados, extrados)
        shape.parameters = {"type": "general", "thk": thk, "xy_span": xy_span, "treat_creases": treat_creases}
        return shape

    @classmethod
    def from_pointcloud_and_topology(cls, meshtopology: Mesh, intrados_pts, extrados_pts, middle_pts=None, fill_pts=None) -> "Shape":
        """Construct a Shape from a pointcloud and a given topology.

        Parameters
        ----------
        meshtopology: :class:`~compas.datastructures.Mesh`
            Mesh with the topology to be used.
        intrados_pts : list
            List of points collected in the intrados.
        extrados_pts : list
            List of points collected in the extrados.
        middle_pts : list, optional
            List of points collected in the middle.
            The default value is ``None``, in which case no middle surface is constructed
        fill_pts : list, optional
            List of points collected in the fill.
            The default value is ``None``, in which case no fill surface is constructed

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        intrados_mesh = create_mesh_from_topology_and_pointcloud(meshtopology, intrados_pts)
        extrados_mesh = create_mesh_from_topology_and_pointcloud(meshtopology, extrados_pts)
        middle_mesh = None

        if middle_pts:
            middle_mesh = create_mesh_from_topology_and_pointcloud(meshtopology, middle_pts)

        if fill_pts:
            fill_mesh = create_mesh_from_topology_and_pointcloud(meshtopology, fill_pts)

        shape = cls.from_meshes(intrados_mesh, extrados_mesh, middle=middle_mesh)
        shape.fill = fill_mesh  # this is not very consistent
        shape.parameters = {"type": "general"}
        return shape

    @classmethod
    def from_meshes_and_formdiagram(cls, form: "FormDiagram", intrados: Mesh, extrados: Mesh, middle: Optional[Mesh] = None) -> "Shape":
        """Construct a Shape from a pair of meshes and a formdiagram that will have its topology copied and normals copied.

        Parameters
        ----------
        form : :class:`FormDiagram`
            Form Diagram with the topology to be used
        intrados : :class:`~compas.datastructures.Mesh`
            Mesh for intrados
        extrados : :class:`~compas.datastructures.Mesh`
            Mesh for extrados
        middle : :class:`~compas.datastructures.Mesh`, optional
            Mesh for middle.
        data : dict,  optional
            Dictionary with the data in required.

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        # this would be better
        # intrados_mesh = MeshDos.from_topology_and_mesh(form, intrados, keep_normals=True)
        # extrados_mesh = MeshDos.from_topology_and_mesh(form, extrados, keep_normals=True)

        intrados_mesh = create_mesh_from_topology_and_basemesh(form, intrados, cls=MeshDos)
        extrados_mesh = create_mesh_from_topology_and_basemesh(form, extrados, cls=MeshDos)

        shape = cls.from_meshes(intrados_mesh, extrados_mesh, middle=middle)
        shape.parameters = {"type": "general"}

        return shape

    @classmethod
    def from_formdiagram_and_attributes(cls, form: "FormDiagram") -> "Shape":
        """Construct a Shape from the form diagram and its attributes 'ub' and 'lb'.

        Parameters
        ----------
        form : :class:`FormDiagram`
            Form Diagram with the topology to be used.

        Returns
        -------
        :class:`Shape`
            A Shape object.

        """
        vertices, faces = form.to_vertices_and_faces()
        intra = MeshDos.from_vertices_and_faces(vertices, faces)
        extra = MeshDos.from_vertices_and_faces(vertices, faces)
        middle = MeshDos.from_vertices_and_faces(vertices, faces)

        for key in form.vertices():
            intra.vertex_attribute(key, "z", form.vertex_attribute(key, "lb"))
            extra.vertex_attribute(key, "z", form.vertex_attribute(key, "ub"))
            middle.vertex_attribute(key, "z", form.vertex_attribute(key, "target"))

        shape = cls.from_meshes(intra, extra, middle=middle)
        shape.parameters = {"type": "general"}
        return shape

    @classmethod
    def from_assembly(self):
        """Create a TNO Shape from an assembly

        Returns
        -------
        Shape
            The shape object
        """

        raise NotImplementedError

    # @classmethod
    # def from_middle_pointcloud(cls, middle_pts, topology=None, thk=0.50, treat_creases=False, printout=False,
    #                            data={'type': 'general', 't': 0.0, 'xy_span': [[0.0, 10.0], [0.0, 10.0]]}):
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

    # =============================================================================
    # Factories
    # =============================================================================

    @classmethod
    def create_dome(cls, center=[5.0, 5.0, 0.0], radius=5.0, thk=0.5, discretisation=[16, 40], t=0.0, **kwargs) -> "Shape":
        """Create the shape representing a Hemispheric Dome

        Parameters
        ----------
        center : [float, float, float], optional
            Center of the dome, by default [5.0, 5.0, 0.0]
        radius : float, optional
            Central radius of the dome, by default 5.0
        thk : float, optional
            Thickness of the dome, by default 0.5
        discretisation : [float, float], optional
            Discretisation of the parallel and meridians, by default [16, 40]
        t : float, optional
            Negative thickness to consider if the intrados can not be interpolated in the point, by default 0.5

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        intrados, extrados, middle = set_dome_heighfield(center, radius=radius, thk=thk, discretisation=discretisation, t=t, expanded=False)
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "dome",
            "thk": thk,
            "discretisation": discretisation,
            "center": center,
            "radius": radius,
            "t": t,
            "expanded": False,
        }
        return shape

    @classmethod
    def create_dome_polar(cls, center=[5.0, 5.0, 0.0], radius=5.0, thk=0.5, discretisation=[16, 40], t=0.0, **kwargs) -> "Shape":
        """Create the shape representing a Hemispheric Dome

        Parameters
        ----------
        center : [float, float, float], optional
            Center of the dome, by default [5.0, 5.0, 0.0]
        radius : float, optional
            Central radius of the dome, by default 5.0
        thk : float, optional
            Thickness of the dome, by default 0.5
        discretisation : [float, float], optional
            Discretisation of the parallel and meridians, by default [16, 40]
        t : float, optional
            Negative thickness to consider if the intrados can not be interpolated in the point, by default 0.0

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        intrados, extrados, middle = set_dome_polar_coord(center, radius=radius, thk=thk, discretisation=discretisation, t=t)
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "dome_polar",
            "thk": thk,
            "discretisation": discretisation,
            "center": center,
            "radius": radius,
            "t": t,
            "expanded": False,
        }
        return shape

    @classmethod
    def create_arch(cls, H=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, t=0.5, discretisation=100, **kwargs) -> "Shape":
        """Create the shape representing a circular arch.

        Parameters
        ----------
        H : float, optional
            Height of the arch, by default 1.00
        L : float, optional
            Span of the arch, by default 2.0
        x0 : float, optional
            Starting coordinate of the arch , by default 0.0
        thk : float, optional
            Thickness of the arch, by default 0.20
        b : float, optional
            Out-of-plane measure of the arch, by default 0.5
        t : float, optional
            Parameter for lower bound in nodes in the boundary, by default 0.0
        discretisation : int, optional
            Density of the shape, by default 100

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        intrados, extrados, middle = arch_shape(H=H, L=L, x0=x0, thk=thk, total_nodes=discretisation, b=b, t=t)
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "arch",
            "thk": thk,
            "discretisation": discretisation,
            "H": H,
            "L": L,
            "x0": x0,
            "b": b,
            "t": t,
        }
        return shape

    @classmethod
    def create_arch_polar(cls, H=1.00, L=2.0, x0=0.0, thk=0.20, b=0.5, t=0.0, discretisation=100, **kwargs) -> "Shape":
        """Create the shape representing a circular arch.

        Parameters
        ----------
        H : float, optional
            Height of the arch, by default 1.00
        L : float, optional
            Span of the arch, by default 2.0
        x0 : float, optional
            Starting coordinate of the arch , by default 0.0
        thk : float, optional
            Thickness of the arch, by default 0.20
        b : float, optional
            Out-of-plane measure of the arch, by default 0.5
        t : float, optional
            Parameter for lower bound in nodes in the boundary, by default 0.0
        discretisation : int, optional
            Density of the shape, by default 100

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        intrados, extrados, middle = arch_shape_polar(H=H, L=L, x0=x0, thk=thk, total_nodes=discretisation, b=b)
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "arch_polar",
            "thk": thk,
            "discretisation": discretisation,
            "H": H,
            "L": L,
            "x0": x0,
            "b": b,
            "t": t,
        }
        return shape

    @classmethod
    def create_crossvault(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.50, t=0.0, spr_angle=None, discretisation=[100, 100], **kwargs) -> "Shape":
        """Create the shape representing a Crossvault

        Parameters
        ----------
        xy_span : [list[float], list[float]], optional
            xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
        thk : float, optional
            The thickness of the vault, by default 0.5
        t : float, optional
            Parameter for lower bound in nodes in the boundary, by default 0.0
        spr_angle : float, optional
            Springing angle, by default None
        discretisation : list|int, optional
            Level of discretisation of the shape, by default [100, 100]

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        if spr_angle:  # amplify the bounds if springing angle is present
            alpha = 1 / math.cos(math.radians(spr_angle))
            x0, xf = xy_span[0]
            Ldiff = (xf - x0) * (alpha - 1)
            xy_span = [[x0 - Ldiff / 2, xf + Ldiff / 2], [x0 - Ldiff / 2, xf + Ldiff / 2]]

        intrados, extrados, middle = cross_vault_highfields(
            xy_span,
            thk=thk,
            discretisation=discretisation,
            t=t,
            expanded=False,
        )
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "crossvault",
            "thk": thk,
            "discretisation": discretisation,
            "xy_span": xy_span,
            "t": t,
            "expanded": False,
        }
        return shape

    @classmethod
    def create_pointedcrossvault(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.5, discretisation=[10, 10], hc=8.0, he=None, hm=None, t=0.0, **kwargs) -> "Shape":
        """Create the shape representing a Pointed Crossvault

        Parameters
        ----------
        xy_span : list, optional
            [description], by default [[0.0, 10.0], [0.0, 10.0]]
        thk : float, optional
            [description], by default 0.5
        discretisation : list, optional
            [description], by default [10, 10]
        hc : float, optional
            Height in the middle point of the vault, by default 8.0
        he : [float, float, float, float], optional
            Height of the opening mid-span for each of the quadrants, by default None
        hm : [float, float, float, float], optional
            Height of each quadrant center (spadrel), by default None
        t : float, optional
            Parameter for lower bound in nodes in the boundary, by default 0.0

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        intrados, extrados, middle = pointed_vault_heightfields(xy_span=xy_span, discretisation=discretisation, t=t, hc=hc, he=he, hm=hm, thk=thk)
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "pointed_crossvault",
            "xy_span": xy_span,
            "thk": thk,
            "discretisation": discretisation,
            "t": t,
            "hc": hc,
            "he": he,
            "hm": hm,
        }
        return shape

    @classmethod
    def create_pavillionvault(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], thk=0.50, t=0.0, discretisation=[100, 100], spr_angle=0.0, expanded=False, **kwargs) -> "Shape":
        """Create the shape representing a Pavillion Vault

        Parameters
        ----------
        xy_span : [list[float], list[float]], optional
            xy-span of the shape, by default [[0.0, 10.0], [0.0, 10.0]]
        thk : float, optional
            The thickness of the vault, by default 0.5
        t : float, optional
            Parameter for lower bound in nodes in the boundary, by default 0.0
        discretisation : list|int, optional
            Level of discretisation of the shape, by default [100, 100]
        spr_angle : float, optional
            Springing angle, by default 0.0
        expanded : bool, optional
            If the extrados should extend beyond the floor plan, by default False

        Returns
        -------
        shape : Shape
            The shape of the dome.

        """
        intrados, extrados, middle = pavillion_vault_highfields(xy_span, thk=thk, discretisation=discretisation, t=t, spr_angle=spr_angle, expanded=expanded)
        shape = cls(middle, intrados, extrados)
        shape.parameters = {
            "type": "pavillionvault",
            "thk": thk,
            "discretisation": discretisation,
            "xy_span": xy_span,
            "t": t,
            "spr_angle": spr_angle,
            "expanded": expanded,
        }
        return shape

    # =============================================================================
    # Methods
    # =============================================================================

    def interpolate_middle_from_ub_lb(self, intrados=None, extrados=None) -> None:
        """Interpolate the middle surface based on intrados and extrados

        Parameters
        ----------
        intrados : Mesh, optional
            Intrados to consider, by default None
        extrados : Mesh, optional
            Extrados to consider, by default None

        Returns
        -------
        None
            Middle surface is added to the Shape.

        """

        if not intrados:
            intrados = self.intrados

        if not extrados:
            extrados = self.extrados

        self.middle = MeshDos.from_mesh(self.intrados.copy())

        for key in self.intrados.vertices():
            ub = self.extrados.vertex_attribute(key, "z")
            lb = self.intrados.vertex_attribute(key, "z")
            self.middle.vertex_attribute(key, "z", (ub + lb) / 2)

    def store_normals(self, mark_fixed_LB=True) -> None:
        """Store the normals of the shape

        Parameters
        ----------
        mark_fixed_LB : bool, optional
            If vertices in the boundary are taken specially, by default True

        Returns
        -------
        None

        """
        intrados = self.intrados
        extrados = self.extrados

        intrados.store_normals()
        extrados.store_normals()

        intrados.vertices_attribute("is_outside", False)
        extrados.vertices_attribute("is_outside", False)

        if mark_fixed_LB:  # Look for the vertices with LB == 't'
            try:
                t = self.parameters["t"]
            except BaseException:
                print("No t is assigned to vertices in intrados.")
                t = None
            if t is not None:
                for key in intrados.vertices():
                    if abs(intrados.vertex_attribute(key, "z") - t) < 10e-3:
                        intrados.vertex_attribute(key, "is_outside", True)
                    if abs(extrados.vertex_attribute(key, "z") - t) < 10e-3:
                        extrados.vertex_attribute(key, "is_outside", True)

    def analytical_normals(self, assume_shape=None, mark_fixed_LB=True):
        """Compute normals analytically for some shapes

        Parameters
        ----------
        assume_shape : dict, optional
            A dictionary with data of the shape to assume, by default None
        mark_fixed_LB : bool, optional
            If vertices in the boundary should be considered specially, by default True

        Returns
        -------
        None

        """

        if assume_shape:
            data = assume_shape
        else:
            data = self.parameters

        if data["type"] == "dome":
            xc, yc = data["center"]
        else:
            raise Exception

        intrados = self.intrados
        extrados = self.extrados

        for key in intrados.vertices():
            x, y, zlb = intrados.vertex_coordinates(key)
            _, _, zub = extrados.vertex_coordinates(key)
            nlb = normalize_vector([(x - xc), (y - yc), zlb])
            nub = normalize_vector([(x - xc), (y - yc), zub])
            intrados.vertex_attribute(key, "n", nlb)
            extrados.vertex_attribute(key, "n", nub)

        intrados.vertices_attribute("is_outside", False)
        extrados.vertices_attribute("is_outside", False)

        if mark_fixed_LB:  # Look for the vertices with LB == 't'
            try:
                t = self.parameters["t"]
            except BaseException:
                print("No t is assigned to vertices in intrados.")
                t = None
            if t is not None:
                for key in intrados.vertices():
                    if abs(intrados.vertex_attribute(key, "z") - t) < 10e-3:
                        intrados.vertex_attribute(key, "is_outside", True)
                    if abs(extrados.vertex_attribute(key, "z") - t) < 10e-3:
                        extrados.vertex_attribute(key, "is_outside", True)

    def compute_selfweight(self) -> float:
        """Compute and returns the total selfweight of the structure based on the area and thickness in the data.

        Returns
        -------
        float

        """

        middle = self.middle
        ro = self.ro
        try:
            thk = self.parameters["thk"]
            area = middle.area()
            total_selfweight = thk * area * ro
        except BaseException:
            intrados = self.intrados
            extrados = self.extrados
            total_selfweight = 0
            for key in intrados.vertices():
                h = extrados.vertex_attribute(key, "z") - intrados.vertex_attribute(key, "z")
                vol = h * middle.vertex_projected_area(key)
                total_selfweight += vol * ro

        return total_selfweight

    def compute_volume(self) -> float:
        """Compute and returns the vollume of the structure based on the area and thickness in the data.

        Returns
        -------
        float

        """

        middle = self.middle
        thk = self.parameters["thk"]
        area = middle.area()
        volume = thk * area

        return volume

    def add_fill_with_height(self, height, fill_ro=20.0):
        """Compute the weight of the fill applied based on the height.

        Parameters
        ----------
        height : float
            The height until which infill should be observed.
        fill_ro : float, optional
            The density of the infill, the default is 20.0

        """

        self.extrados_fill = deepcopy(self.extrados)
        volume = 0.0
        proj_area_total = 0.0

        for key in self.extrados_fill.vertices():
            zi = self.extrados_fill.vertex_attribute(key, "z")
            if zi < height:
                self.extrados_fill.vertex_attribute(key, "z", height)
                # this should become a method
                area = 0.0
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
                vol_i = proj_area * (height - zi)
                volume += vol_i

        self.fill = True
        self.fill_volume = volume
        print("Proj area total of shape", proj_area_total)
        self.fill_ro = fill_ro

    def compute_fill_weight(self) -> float:
        """Compute and returns the volume of fill in the structure.

        Returns
        -------
        float

        """
        return self.fill_volume * self.fill_ro

    # =============================================================================
    # Moved methods
    # =============================================================================

    def get_ub(self, x, y):
        """Get the height of the extrados in the point.

        Parameters
        ----------
        x : float
            x-coordinate of the point to evaluate.
        y : float
            y-coordinate of the point to evaluate.

        Returns
        -------
        z : float
            The extrados evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.extrados.vertices_attributes("xyz"), [x, y], method=method)

    def get_ub_pattern(self, XY):
        """Get the height of the extrados in a list of points.

        Parameters
        ----------
        XY : list or array
            list of the x-coordinate and y-coordinate of the points to evaluate.

        Returns
        -------
        z : float
            The extrados evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.extrados.vertices_attributes("xyz"), XY, method=method)

    def get_ub_fill(self, x, y):
        """Get the height of the fill in the point.

        Parameters
        ----------
        x : float
            x-coordinate of the point to evaluate.
        y : float
            y-coordinate of the point to evaluate.

        Returns
        -------
        z : float
            The extrados evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.extrados_fill.vertices_attributes("xyz"), [x, y], method=method)

    def get_lb(self, x, y):
        """Get the height of the intrados in the point.

        Parameters
        ----------
        x : float
            x-coordinate of the point to evaluate.
        y : float
            y-coordinate of the point to evaluate.

        Returns
        -------
        z : float
            The intrados evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.intrados.vertices_attributes("xyz"), [x, y], method=method)

    def get_lb_pattern(self, XY):
        """Get the height of the intrados in a list of points.

        Parameters
        ----------
        XY : list or array
            list of the x-coordinate and y-coordinate of the points to evaluate.

        Returns
        -------
        z : float
            The extrados evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.intrados.vertices_attributes("xyz"), XY, method=method)

    def get_middle(self, x, y):
        """Get the height of the target/middle surface in the point.

        Parameters
        ----------
        x : float
            x-coordinate of the point to evaluate.
        y : float
            y-coordinate of the point to evaluate.

        Returns
        -------
        z : float
            The middle surface evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.middle.vertices_attributes("xyz"), [x, y], method=method)

    def get_middle_pattern(self, XY):
        """Get the height of the target/middle surface in a list of points.

        Parameters
        ----------
        XY : list or array
            list of the x-coordinate and y-coordinate of the points to evaluate.

        Returns
        -------
        z : float
            The extrados evaluated in the point.

        """
        method = self.parameters.get("interpolation", "linear")
        return interpolate_from_pointcloud(self.middle.vertices_attributes("xyz"), XY, method=method)
