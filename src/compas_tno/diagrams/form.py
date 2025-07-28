import math
from typing import Annotated
from typing import Optional

import numpy.typing as npt
from numpy import vstack
from numpy import zeros

from compas.datastructures import Mesh
from compas.geometry import bounding_box_xy
from compas_tna.diagrams import FormDiagram as TNAFormDiagram
from compas_tno.shapes import Shape

# this needs to be internalised by the shape class
from compas_tno.shapes import arch_ub_lb_update
from compas_tno.shapes import crossvault_middle_update
from compas_tno.shapes import crossvault_ub_lb_update
from compas_tno.shapes import dome_ub_lb_update
from compas_tno.shapes import dome_zt_update
from compas_tno.shapes import pavillionvault_ub_lb_update
from compas_tno.shapes import pointed_arch_ub_lb_update
from compas_tno.shapes import pointed_vault_middle_update
from compas_tno.shapes import pointed_vault_ub_lb_update

from .diagram_arch import create_arch_form_diagram
from .diagram_arch import create_linear_form_diagram
from .diagram_arch import create_linear_form_diagram_sp_ep
from .diagram_circular import create_circular_radial_form
from .diagram_circular import create_circular_radial_spaced_form
from .diagram_circular import create_circular_spiral_form
from .diagram_rectangular import create_cross_diagonal
from .diagram_rectangular import create_cross_form
from .diagram_rectangular import create_cross_with_diagonal
from .diagram_rectangular import create_delta_form
from .diagram_rectangular import create_fan_form
from .diagram_rectangular import create_ortho_form
from .diagram_rectangular import create_parametric_form


class FormDiagram(TNAFormDiagram):
    """The FormDiagram represents the projection of the internal force network within the structure.

    Attributes
    ----------
    args : list
    parameters : dict

    Notes
    -----
    Vertices have the following default attributes.

    * "x" (float): x-coordinate of the vertex.
    * "y" (float): y-coordinate of the vertex.
    * "z" (float): z-coordinate of the vertex.
    * "is_roller" (bool): Whether or not the attibute is set to be a roller. Takes only horizontal reaction.
    * "tubmax" (float): Maximum allowable increase in the extrados
    * "tlbmax" (float): Maximum allowable increase in the intradoss
    * "tub_reac" (float): Maximum allowable increase in the horizontal projection of the reaction.

    Edges have the following default attributes.

    * q (float): Force density of an edge.
    * f (float): Force magnitude of an edge.
    * is_ind (bool): Whether or not the edge is an independent edge.
    * qmin (float): The minimum force density magnitude.
    * qmax (float): The maximum force density magnitude.

    Examples
    --------
    >>> from compas_tno.diagrams import FormDiagram
    >>> data = {"type": "fan_fd", "xy_span": [[0, 10], [0, 10]], "discretisation": [10, 10], "fix": "corners"}
    >>> form = FormDiagram.from_library(data)

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.update_default_vertex_attributes(
            {
                "is_roller": False,
                "tubmax": None,
                "tlbmax": None,
                "tub_reac": None,
            }
        )
        self.update_default_edge_attributes(
            {
                "q": -1.0,
                "f": -1.0,
                "is_symmetry": False,
                "is_ind": False,
                "qmin": -1e4,
                "qmax": 1e-8,
            }
        )
        self.args = None
        self.attributes["loadpath"] = 0.0
        self.attributes["indset"] = None
        self.parameters = {
            # "type": None,
            # "discretisation": None,
            # "E": None,
            # "Ah": None,
            # # 'x0': None,
            # # 'discretisation': 20,
            # # 'xy_span': [[0,10],[0,10]],
            # # 'discretisation': [10,10],
            # # 'fix': 'corners',
            # # 'r_oculus:' None,
        }

    # =============================================================================
    # Constructors
    # =============================================================================

    @classmethod
    def from_library(cls, data: dict) -> "FormDiagram":
        """Construct a FormDiagram the library of FormDiagrams available.

        Parameters
        ----------
        data : dict
            Dictionary with information to generate the FormDiagram.

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form_type = data.pop("type")

        if form_type == "arch":
            form = cls.create_arch(**data)
        elif form_type == "linear_arch" or form_type == "pointed_arch":
            form = cls.create_linear_form_diagram(**data)
        elif form_type == "linear_arch_sp_ep" or form_type == "pointed_arch":
            form = cls.create_linear_form_diagram_sp_ep(**data)
        elif form_type == "cross_fd":
            form = cls.create_cross_form(**data)
        elif form_type == "cross_diagonal":
            form = cls.create_cross_diagonal(**data)
        elif form_type == "cross_with_diagonal":
            form = cls.create_cross_with_diagonal(**data)
        elif form_type == "fan_fd":
            form = cls.create_fan_form(**data)
        elif form_type == "ortho" or form_type == "ortho_fd":
            form = cls.create_ortho_form(**data)
        elif form_type == "radial_fd":
            form = cls.create_circular_radial_form(**data)
        elif form_type == "radial_spaced_fd":
            form = cls.create_circular_radial_spaced_form(**data)
        elif form_type == "spiral_fd":
            form = cls.create_circular_spiral_form(**data)
        elif form_type == "parametric_form":
            form = cls.create_parametric_form(**data)
        elif form_type == "delta_form":
            form = cls.create_delta_form(**data)
        else:
            raise NotImplementedError

        return form

    @classmethod
    def create_arch(cls, H: float = 1.00, L: float = 2.00, x0: float = 0.0, discretisation: int = 100) -> "FormDiagram":
        """Construct a FormDiagram based on an arch linear discretisation.
        Note: The nodes of the form diagram are spaced following a projection in a semicircular arch.

        Parameters
        ----------
        H : float, optional
            Height of the arch, by default 1.00
        L : float, optional
            Span of the arch, by default 2.00
        x0 : float, optional
            Initial coordiante of the arch, by default 0.0
        discretisation : int, optional
            Numbers of nodes to be considered in the form diagram, by default 100

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_arch_form_diagram(cls, H=H, L=L, x0=x0, discretisation=discretisation)
        form.parameters = {
            "type": "arch",
            "L": L,
            "H": H,
            "x0": x0,
            "discretisation": discretisation,
        }
        return form

    @classmethod
    def create_linear_form_diagram(cls, L: float = 2.0, x0: float = 0.0, discretisation: int = 100) -> "FormDiagram":
        """Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes.

        Parameters
        ----------
        L : float, optional
            Span of the arch, by default 2.00
        x0 : float, optional
            Initial coordiante of the arch, by default 0.0
        discretisation : int, optional
            Numbers of nodes to be considered in the form diagram, by default 100

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            FormDiagram generated according to the parameters.
        """
        form: FormDiagram = create_linear_form_diagram(cls, L=L, x0=x0, discretisation=discretisation)
        form.parameters = {
            "type": "linear_arch",
            "L": L,
            "x0": x0,
            "discretisation": discretisation,
        }
        return form

    @classmethod
    def create_linear_form_diagram_sp_ep(cls, sp=[0, 0, 0], ep=[2, 0, 0], discretisation=100) -> "FormDiagram":
        """Helper to create a arch linear form-diagram with equaly spaced (in 2D) nodes.

        Parameters
        ----------
        sp : list, optional
            Starting point coordinates, by default [0, 0, 0]
        sp : list, optional
            End point coordinates, by default [2, 0, 0]
        discretisation : int, optional
            Numbers of nodes to be considered in the form diagram, by default 100

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            FormDiagram generated according to the parameters.
        """
        form: FormDiagram = create_linear_form_diagram_sp_ep(cls, sp=sp, ep=ep, discretisation=discretisation)
        form.parameters = {
            "type": "linear_arch_sp_ep",
            "ep": ep,
            "sp": sp,
            "discretisation": discretisation,
        }
        return form

    @classmethod
    def create_cross_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=10, fix="corners") -> "FormDiagram":
        """Construct a FormDiagram based on cross discretiastion with orthogonal arrangement and diagonal.

        Parameters
        ----------
        xy_span : list, optional
            List with initial- and end-points of the vault, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted., by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_cross_form(cls, xy_span=xy_span, discretisation=discretisation, fix=fix)
        form.parameters = {
            "type": "cross_fd",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "fix": fix,
        }
        return form

    @classmethod
    def create_cross_diagonal(
        cls,
        xy_span=[[0.0, 10.0], [0.0, 10.0]],
        discretisation=10,
        partial_bracing_modules=None,
        fix="corners",
    ) -> "FormDiagram":
        """Construct a FormDiagram based on a mixture of cross and fan discretiastion

        Parameters
        ----------
        xy_span : list, optional
            List with initial- and end-points of the vault, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        partial_bracing_modules : str, optional
            If partial bracing modules are included, by default None
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted, by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_cross_diagonal(cls, xy_span=xy_span, partial_bracing_modules=partial_bracing_modules, discretisation=discretisation, fix=fix)
        form.parameters = {
            "type": "cross_diagonal",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "partial_bracing_modules": partial_bracing_modules,
            "fix": fix,
        }
        return form

    @classmethod
    def create_cross_with_diagonal(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=10, fix="corners") -> "FormDiagram":
        """Construct a FormDiagram based on cross discretiastion with diagonals.

        Parameters
        ----------
        xy_span : list, optional
            List with initial- and end-points of the vault, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted, by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_cross_with_diagonal(cls, xy_span=xy_span, discretisation=discretisation, fix=fix)
        form.parameters = {
            "type": "cross_with_diagonal",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "fix": fix,
        }
        return form

    @classmethod
    def create_fan_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10, 10], fix="corners") -> "FormDiagram":
        """Helper to construct a FormDiagram based on fan discretiastion with straight lines to the corners.

        Parameters
        ----------
        xy_span : list, optional
            List with initial- and end-points of the vault, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted, by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_fan_form(cls, xy_span=xy_span, discretisation=discretisation, fix=fix)
        form.parameters = {
            "type": "fan_fd",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "fix": fix,
        }
        return form

    @classmethod
    def create_parametric_form(
        cls,
        xy_span=[[0.0, 10.0], [0.0, 10.0]],
        discretisation=10,
        lambd=0.5,
        fix="corners",
    ) -> "FormDiagram":
        """Helper to construct a FormDiagram envelop that includes fan and cross diagrams. The diagram is modified by the parameter lambda.

        Parameters
        ----------
        xy_span : [[float, float], [float, float]], optional
            List with initial- and end-points of the vault, by default, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        lambd : float, optional
            Inclination of the arches in the diagram (0.0 will result in cross and 1.0 in fan diagrams), by default 0.5
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted, by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_parametric_form(cls, xy_span=xy_span, discretisation=discretisation, lambd=lambd, fix=fix)
        form.parameters = {
            "type": "parametric_form",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "fix": fix,
            "lambd": lambd,
        }
        return form

    @classmethod
    def create_delta_form(
        cls,
        xy_span=[[0.0, 10.0], [0.0, 10.0]],
        discretisation=10,
        delta=0.5,
        fix="corners",
    ) -> "FormDiagram":
        """Helper to construct a FormDiagram based on fan discretiastion with straight lines to the corners.

        Parameters
        ----------
        xy_span : [[float, float], [float, float]], optional
            List with initial- and end-points of the vault, by default, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        lambd : float, optional
            Inclination of the arches in the diagram (0.0 will result in cross and 1.0 in fan diagrams), by default 0.5
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted, by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_delta_form(cls, xy_span=xy_span, discretisation=discretisation, delta=delta, fix=fix)
        form.parameters = {
            "type": "delta_form",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "fix": fix,
            "delta": delta,
        }
        return form

    @classmethod
    def create_ortho_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10, 10], fix="corners") -> "FormDiagram":
        """Helper to construct a FormDiagram based on a simple orthogonal discretisation.

        Parameters
        ----------
        xy_span : list, optional
            List with initial- and end-points of the vault, by default [[0.0, 10.0], [0.0, 10.0]]
        discretisation : int, optional
            Set the density of the grid in x and y directions, by default 10
        fix : str, optional
            Option to select the constrained nodes: 'corners', 'all' are accepted, by default 'corners'

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_ortho_form(cls, xy_span=xy_span, discretisation=discretisation, fix=fix)
        form.parameters = {
            "type": "ortho_fd",
            "xy_span": xy_span,
            "discretisation": discretisation,
            "fix": fix,
        }
        return form

    @classmethod
    def create_circular_radial_form(
        cls,
        center=[5.0, 5.0],
        radius=5.0,
        discretisation=[8, 20],
        r_oculus=0.0,
        diagonal=False,
        partial_diagonal=False,
    ) -> "FormDiagram":
        """Construct a circular radial FormDiagram with hoops not equally spaced in plan.

        Parameters
        ----------
        center : list, optional
            Planar coordinates of the form-diagram [xc, yc], by default [5.0, 5.0]
        radius : float, optional
            Radius of the form diagram, by default 5.0
        discretisation : list, optional
            Number of hoops, and of parallels of the dome form diagram], by default [8, 20]
        r_oculus : float, optional
            Value of the radius of the oculus, if no oculus is present should be set to zero, by default 0.0
        diagonal : bool, optional
            Activate diagonal in the quads, by default False
        partial_diagonal : bool, optional
            Activate partial diagonal in the quads, by default False

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_circular_radial_form(
            cls, center=center, radius=radius, discretisation=discretisation, r_oculus=r_oculus, diagonal=diagonal, partial_diagonal=partial_diagonal
        )
        form.parameters = {
            "type": "radial_fd",
            "center": center,
            "radius": radius,
            "discretisation": discretisation,
            "r_oculus": r_oculus,
            "diagonal": diagonal,
            "partial_diagonal": partial_diagonal,
        }
        return form

    @classmethod
    def create_circular_radial_spaced_form(
        cls,
        center=[5.0, 5.0],
        radius=5.0,
        discretisation=[8, 20],
        r_oculus=0.0,
        diagonal=False,
        partial_diagonal=False,
    ) -> "FormDiagram":
        """Construct a circular radial FormDiagram with hoops not equally spaced in plan,
        but equally spaced with regards to the projection on a hemisphere.

        Parameters
        ----------
        center : list, optional
            Planar coordinates of the form-diagram [xc, yc], by default [5.0, 5.0]
        radius : float, optional
            Radius of the form diagram, by default 5.0
        discretisation : list, optional
            Number of hoops, and of parallels of the dome form diagram], by default [8, 20]
        r_oculus : float, optional
            Value of the radius of the oculus, if no oculus is present should be set to zero, by default 0.0
        diagonal : bool, optional
            Activate diagonal in the quads, by default False
        partial_diagonal : bool, optional
            Activate partial diagonal in the quads, by default False

        Returns
        -------
        :class:`~compas_tno.diagrams.FormDiagram`
            The FormDiagram created.

        """
        form: FormDiagram = create_circular_radial_spaced_form(
            cls, center=center, radius=radius, discretisation=discretisation, r_oculus=r_oculus, diagonal=diagonal, partial_diagonal=partial_diagonal
        )
        form.parameters = {
            "type": "radial_spaced_fd",
            "center": center,
            "radius": radius,
            "discretisation": discretisation,
            "r_oculus": r_oculus,
            "diagonal": diagonal,
            "partial_diagonal": partial_diagonal,
        }
        return form

    @classmethod
    def create_circular_spiral_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0) -> "FormDiagram":
        """Construct a circular radial FormDiagram with hoops not equally spaced in plan,
        but equally spaced with regards to the projection on a hemisphere.

        Parameters
        ----------
        center : list, optional
            Planar coordinates of the form-diagram [xc, yc], by default [5.0, 5.0]
        radius : float, optional
            Radius of the form diagram, by default 5.0
        discretisation : list, optional
            Number of hoops, and of parallels of the dome form diagram], by default [8, 20]
        r_oculus : float, optional
            Value of the radius of the oculus, if no oculus is present should be set to zero, by default 0.0

        Returns
        -------
        FormDiagram
            The FormDiagram created.

        """
        form: FormDiagram = create_circular_spiral_form(cls, center=center, radius=radius, discretisation=discretisation, r_oculus=r_oculus)
        form.parameters = {
            "type": "spiral_fd",
            "center": center,
            "radius": radius,
            "discretisation": discretisation,
            "r_oculus": r_oculus,
        }
        return form

    @classmethod
    def from_assembly(self):
        NotImplementedError

    # @classmethod
    # def from_singular(cls, json_path):
    #     """Create form diagram from singular dense mesh.

    #     Parameters
    #     ----------
    #     json_path : str
    #         Path of the mesh saved in JSON format.

    #     Returns
    #     -------
    #     :class:`~compas_tno.diagrams.FormDiagram`
    #         The FormDiagram created.
    #     """

    #     from compas.datastructures import mesh_delete_duplicate_vertices

    #     mesh = Mesh.from_json(json_path)
    #     mesh_delete_duplicate_vertices(mesh)
    #     form = cls().from_mesh(mesh)
    #     form.set_boundary_supports()
    #     form.delete_boundary_edges()

    #     return form

    # @classmethod
    # def from_triangle(cls, boundary_points, area=0.5, angle=30):
    #     """Generates a Form Diagram from triangle

    #     Parameters
    #     ----------
    #     boundary_points : [list]
    #         Boundary point of pattern
    #     area : float, optional
    #         Average area of the trangles, by default 0.5
    #     angle : int, optional
    #         minimum angle in the triangless, by default 30

    #     Returns
    #     -------
    #     :class:`~compas_tno.diagrams.FormDiagram`
    #         The FormDiagram created.
    #     """

    #     from compas_triangle.delaunay import conforming_delaunay_triangulation

    #     vertices, faces = conforming_delaunay_triangulation(
    #         boundary_points + boundary_points[:1], angle=angle, area=area
    #     )
    #     mesh = Mesh.from_vertices_and_faces(vertices, faces)
    #     form = cls().from_mesh(mesh)

    #     return form

    # @classmethod
    # def from_skeleton(cls, data):
    #     """Create form diagram from compas_skeleton

    #     Parameters
    #     ----------
    #     data : dict
    #         Dictionary with the data for crerating the skeleton mesh.

    #     Returns
    #     -------
    #     :class:`~compas_tno.diagrams.FormDiagram`
    #         The FormDiagram created.
    #     """

    #     from compas_skeleton.datastructure import Skeleton
    #     from compas.geometry import matrix_from_scale_factors
    #     from compas.datastructures import mesh_transform

    #     skeleton_type = data["type"]

    #     if skeleton_type == "circular":

    #         point = data.get("center", [0, 0, 0])
    #         subdivisions = data.get("subdivision", 3)
    #         radius = data.get("radius", 5.0)

    #         skeleton = Skeleton.from_center_point(point)
    #         skeleton.node_width = 7.5
    #         skeleton.update_mesh_vertices_pos()
    #         skeleton.vertex_attribute(0, "z", 10.0)
    #         skeleton.subdivide(subdivisions)
    #         skeleton_mesh = skeleton.to_mesh()

    #         max_x = 0
    #         for key in skeleton_mesh.vertices():
    #             x, _, _ = skeleton_mesh.vertex_coordinates(key)
    #             if x > max_x:
    #                 max_x = x
    #         scale = radius / max_x
    #         T = matrix_from_scale_factors([scale, scale, 0.0])
    #         mesh_transform(skeleton_mesh, T)
    #         form = cls().from_mesh(skeleton_mesh)
    #         form.set_boundary_supports()
    #         form.delete_boundary_edges()

    #     return form

    # =============================================================================
    # Methods - Utility
    # =============================================================================

    def q(self) -> list[float]:
        """Return the force densities of the (real) edges of the diagram."""
        return [self.edge_attribute(edge, "q") for edge in self.edges_where(_is_edge=True)]

    def xy(self) -> list[Annotated[list[float], 2]]:
        """Return the xy coordinates of the diagram."""
        return self.vertices_attributes("xy")

    def fixed(self) -> list[int]:
        """Return the keys of the fixed vertices of the diagram."""
        return list(self.vertices_where(is_fixed=True))

    def fixed_x(self) -> list[int]:
        """Return the keys of the fixed-in-x vertices of the diagram."""
        return list(self.vertices_where(is_fixed_x=True, is_fixed=False))

    def fixed_y(self) -> list[int]:
        """Return the keys of the fixed-in-y vertices of the diagram."""
        return list(self.vertices_where(is_fixed_y=True, is_fixed=False))

    def ind(self) -> list[int]:
        """Return the identifiers for the edges selected as independents on the diagram."""
        return list(self.edges_where(is_ind=True))

    def number_of_real_edges(self) -> int:
        """Return the number of real edges (_is_edge=True) in the diagram."""
        return len(list(self.edges_where(_is_edge=True)))

    def update_f(self) -> None:
        """Update the edge forces based on force density and length."""
        for edge in self.edges():
            f = self.edge_attribute(edge, name="q") * self.edge_length(edge)
            self.edge_attribute(edge, name="f", value=f)

    def add_pz(self, key: int, force: float = 0.0) -> None:
        """Add a vertical load to the loads applied to a specific vertex.

        Parameters
        ----------
        key : int
            Edge identifier
        force : float, optional
            Load to add at the node in the attribute 'pz', by default 0.0

        """
        pz0 = self.vertex_attribute(key, name="pz")
        self.vertex_attribute(key, name="pz", value=pz0 + force)

    def assign_inds(self, edges: Optional[list[tuple[int, int]]] = None, printout: bool = True) -> None:
        """Assign manually independent edges to the pattern

        Parameters
        ----------
        edges : list[tuple[int, int]], optional
            List of edges to be assigned as independents, by default None
        printout : bool, optional
            Whether the number of assigned independent edges should appear, by default True

        """
        if not edges:
            return

        self.edges_attribute(name="is_ind", value=False)
        self.edges_attribute(name="is_ind", value=True, keys=edges)
        self.store_indset(edges)

        if printout:
            print(f"Assigned {len(edges)} independent edges to the pattern.")

    def store_indset(self, edges: Optional[list[tuple[int, int]]] = None) -> None:
        """Store the indset with the information about the independent edges (used within algorithms)

        Parameters
        ----------
        edges : list[tuple[int, int]], optional
            List of edges to be assigned as independents, by default None

        """

        if not edges:
            self.attributes["indset"] = None
            return

        points = []
        for edge in edges:
            point = self.edge_midpoint(edge)
            point.z = 0
            points.append(point)

        self.attributes["indset"] = points

    def update_indset(self) -> None:
        """Update the indset with the information about the independent edges (used within algorithms)"""
        edges = list(self.edges_where(is_ind=True))
        if not edges:
            self.attributes["indset"] = None
            return

        self.store_indset(edges)

    def tributary_dict(self) -> dict[int, dict[int, list[int]]]:
        """Make a tributary dictionary to help in the calculation of the tributary weights to the nodes afterwards.

        Returns
        -------
        dict[int, dict[int, list[int]]]
            The tributary area dictionary

        """
        tributary_dict = {}
        vertex_index = self.vertex_index()

        for key in self.vertices():
            i = vertex_index[key]
            tributary_dict[i] = {}
            for nbr in self.halfedge[key]:
                j = vertex_index[nbr]
                tributary_dict[i][j] = []
                fkey = self.halfedge[key][nbr]
                if fkey is not None:
                    facekeys = self.face_vertices(fkey)
                    faceindices = [vertex_index[v] for v in facekeys]
                    tributary_dict[i][j].append(faceindices)

                fkey = self.halfedge[nbr][key]
                if fkey is not None:
                    facekeys = self.face_vertices(fkey)
                    faceindices = [vertex_index[v] for v in facekeys]
                    tributary_dict[i][j].append(faceindices)

        return tributary_dict

    def tributary_matrices(self, sparse: bool = False) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
        """Returns the Matrices used for computation of the tributary area.

        Parameters
        ----------
        sparse : bool, optional
            If the matrices should be returned as sparse matrices (see ``scipy.linalg.sparse``), by default False

        Returns
        -------
        F : array [f x n]
            Linear transformation from ``X`` [n x 3] to ``c`` [f x 3] with the position of the centroids
        V0 : array [g x n]
            Mark the influence of the original point in the calculation
        V1 : array [g x n]
            Mark the influence of the neighbor points the calculation
        V2 : array [g x f]
            Mark the influence of the centroid points the calculation

        """
        n = self.number_of_vertices()
        f = self.number_of_faces()

        vertex_index = self.vertex_index()
        face_index = {}

        F = zeros((f, n))
        V0 = zeros((0, n))
        V1 = zeros((0, n))
        V2 = zeros((0, f))

        for i, fkey in enumerate(self.faces()):
            face_index[fkey] = i
            faceindices = [vertex_index[v] for v in self.face_vertices(fkey)]
            np = float(len(faceindices))
            for j in faceindices:
                F[i, j] = 1.0 / np

        ig = 0
        for key in self.vertices():
            i = vertex_index[key]
            v0i = zeros((1, n))
            v0i[0, i] = 1.0

            for nbr in self.halfedge[key]:
                j = vertex_index[nbr]

                fkey = self.halfedge[key][nbr]
                if fkey is not None:
                    findex = face_index[fkey]
                    facekeys = self.face_vertices(fkey)
                    faceindices = [vertex_index[v] for v in facekeys]

                    v1i = zeros((1, n))
                    v1i[0, j] = 1.0
                    V1 = vstack([V1, v1i])

                    v2i = zeros((1, f))
                    v2i[0, findex] = 1.0
                    V2 = vstack([V2, v2i])

                    V0 = vstack([V0, v0i])

                    ig += 1

                fkey = self.halfedge[nbr][key]
                if fkey is not None:
                    findex = face_index[fkey]
                    facekeys = self.face_vertices(fkey)
                    faceindices = [vertex_index[v] for v in facekeys]

                    v1i = zeros((1, n))
                    v1i[0, j] = 1.0
                    V1 = vstack([V1, v1i])

                    v2i = zeros((1, f))
                    v2i[0, findex] = 1.0
                    V2 = vstack([V2, v2i])

                    V0 = vstack([V0, v0i])

                    ig += 1

        return F, V0, V1, V2

    def update_lumped_weights(self, thickness: float = 0.5, density: float = 20.0) -> None:
        """Update the lumped weights based on the curent geometry of the thrust network
        The loads are computed based on the tributary area times the thickness times the density.
        For variable thickness the nodal attribute `thk` is considered.

        Parameters
        ----------
        thickness : float, optional
            The thickness of the problem, by default 0.50
            If None is passed, the thickness is taken from the nodal attribute `thk`
        density : float, optional
            The density assumed, by default 20.0

        Returns
        -------
        None
            Loads are updated in place.
        """
        self.apply_selfweight_from_thrust(thickness=thickness, density=density)

    def overview_forces(self) -> None:
        """Print an overview of the forces in the ``Form Diagram``."""

        self.update_f()

        f = []
        q = []
        z = []
        pz = 0

        lp = 0

        for edge in self.edges_where(_is_edge=True):
            qi = self.edge_attribute(edge, name="q")
            li = self.edge_length(edge)
            lp += qi * li**2
            q.append(qi)
            f.append(qi * li)

        print("=" * 20)
        print("Overview on forces:")

        print("q: {0:.3f} : {1:.3f}".format(float(min(q)), float(max(q))))
        print("f: {0:.3f} : {1:.3f}".format(float(min(f)), float(max(f))))
        for key in self.vertices():
            z.append(self.vertex_attribute(key, "z"))
            pz += self.vertex_attribute(key, "pz")
        print("z: {0:.3f} : {1:.3f}".format(float(min(z)), float(max(z))))
        print("pz: {0:.3f}".format(pz))
        print("lp: {0:.3f}".format(lp))
        self.attributes["loadpath"] = lp

    def thrust(self) -> float:
        """Computes the total (horizontal) thrust"""

        thrust = 0
        for key in self.vertices_where(is_fixed=True):
            rx = self.vertex_attribute(key, name="_rx")
            ry = self.vertex_attribute(key, name="_ry")
            r = math.sqrt(rx**2 + ry**2)
            thrust += r

        return thrust

    def lumped_swt(self) -> float:
        """Returns the lumped selfweight on the nodes."""

        swt = 0
        for key in self.vertices():
            pz = self.vertex_attribute(key, "pz")
            swt += pz

        return swt

    def distance_target(self) -> float:
        """Compute the squared vertical distance from a target."""

        f = 0
        for key in self.vertices():
            z = self.vertex_attribute(key, "z")
            s = self.vertex_attribute(key, "target")
            f += (z - s) ** 2

        return f

    def loadpath(self) -> float:
        """Compute loadpath in the current configuration based on the force attributes stored in the ``FormDiagram``."""

        lp = 0
        for edge in self.edges_where(_is_edge=True):
            qi = abs(self.edge_attribute(edge, "q"))
            li = self.edge_length(edge)
            lp += qi * li**2

        self.attributes["loadpath"] = lp

        return lp

    def number_of_independents(self) -> int:
        """Compute the number of independent edges."""

        return len(list(self.edges_where(is_ind=True)))

    # =============================================================================
    # Methods - Boundaries
    # =============================================================================

    def number_of_supports(self) -> int:
        """Compute the number of supports."""

        return len(list(self.vertices_where(is_fixed=True)))

    def delete_boundary_edges(self, delete_corner_vertex: bool = True) -> None:
        """
        Delete boundary edges on a diagram.

        Parameters
        ----------
        delete_corner_vertex : bool, optional
            Delete the corner vertices. The default value is ``True``.

        Returns
        -------
        None
            The ``FormDiagram`` is modified in place.

        """

        for u, v in self.edges_on_boundary():
            if self.vertex_attribute(u, "is_fixed") and self.vertex_attribute(v, "is_fixed"):
                self.edge_attribute((u, v), "_is_edge", False)

        if delete_corner_vertex:
            corners = [key for key in self.corners()]
            for key in corners:
                neighbors = self.vertex_neighbors(key)
                print(neighbors)
                n1 = self.vertex_neighbors(neighbors[0])
                n2 = self.vertex_neighbors(neighbors[1])
                intersection_vertex = list(set.intersection(set(n1), set(n2)))
                sum_vertices = neighbors + intersection_vertex
                new_face = []
                for key_ in sum_vertices:
                    if key_ != key:
                        new_face.append(key_)
                self.delete_vertex(key)
                self.add_face(new_face)
                self.edge_attribute((neighbors[0], neighbors[1]), "_is_edge", False)

    def set_boundary_supports(self) -> None:
        """Set all node on the boundary of a pattern as supports."""

        for key in self.vertices_on_boundary():
            self.vertex_attribute(key, "is_fixed", True)

    def set_boundary_rollers(self, max_rx=[0.0, 0.0], max_ry=[0.0, 0.0], total_rx=None, total_ry=None) -> None:
        """
        Set boundary vertices as rollers.

        Parameters
        ----------
        max_rx : list, optional
            Assign the maximum reaction in x to the vertices in the left and right boundary. The default value is ``[0.0, 0.0]``.
        max_ry : list, optional
            Assign the maximum reaction in y to the vertices in the top and bottom boundary. The default value is ``[0.0, 0.0]``.
        total_rx : list, optional
            Assign the maximum total reaction in x of the vertices in the right or left boundary. The default value is ``None``, in which the ``max_rx`` is used for all vertices.
        total_ry : list, optional
            Assign the maximum total reaction in y of the vertices in the top or bottom boundary.The default value is ``None``, in which the ``max_ry`` is used for all vertices.

        Returns
        -------
        None
            The ``FormDiagram`` is modified in place.

        """
        corners = bounding_box_xy(self.vertices_attributes("xyz"))
        xs = [point[0] for point in corners]
        ys = [point[1] for point in corners]
        xlimits = [min(xs), max(xs)]
        ylimits = [min(ys), max(ys)]
        if not total_rx:
            nx0 = 0
            nx1 = 0
            for key in self.vertices_where(x=xlimits[0]):
                nx0 += 1
            for key in self.vertices_where(x=xlimits[1]):
                nx1 += 1
            max_rx = [total_rx[0] / nx0, total_rx[1] / nx1]
        if not total_ry:
            ny0 = 0
            ny1 = 0
            for key in self.vertices_where(y=ylimits[0]):
                ny0 += 1
            for key in self.vertices_where(y=ylimits[1]):
                ny1 += 1
            max_ry = [total_ry[0] / ny0, total_ry[1] / ny1]
        for key in self.vertices_on_boundary():
            if self.vertex_attribute(key, "is_fixed") is False:
                x, y, _ = self.vertex_coordinates(key)
                if x == xlimits[0]:
                    self.vertex_attribute(key, "rol_x", True)
                    self.vertex_attribute(key, "max_rx", max_rx[0])
                elif x == xlimits[1]:
                    self.vertex_attribute(key, "rol_x", True)
                    self.vertex_attribute(key, "max_rx", max_rx[1])
                elif y == ylimits[0]:
                    self.vertex_attribute(key, "rol_y", True)
                    self.vertex_attribute(key, "max_ry", max_ry[0])
                elif y == ylimits[1]:
                    self.vertex_attribute(key, "rol_y", True)
                    self.vertex_attribute(key, "max_ry", max_ry[1])

    # =============================================================================
    # Methods - Symmetry
    # =============================================================================

    def number_of_sym_edges(self, printout=False) -> int:
        """Compute the symmetric edges."""

        i_sym_max = 0
        for u, v in self.edges_where(_is_edge=True):
            try:
                i_sym = self.edge_attribute((u, v), "sym_key")
            except BaseException:
                print("Warning, no symmetry relation found!")
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        return i_sym_max

    def number_of_sym_vertices(self) -> int:
        """Compute the symmetric vertices."""

        i_sym_max = 0
        for key in self.vertices():
            try:
                i_sym = self.vertex_attribute(key, "sym_key")
            except BaseException:
                print("Warning, no symmetry relation found!")
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        return i_sym_max

    def number_of_sym_supports(self) -> int:
        """Compute the symmetric supports."""

        i_sym_max = 0
        for key in self.vertices_where(is_fixed=True):
            try:
                i_sym = self.vertex_attribute(key, "sym_key")
            except BaseException:
                print("Warning, no symmetry relation found!")
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        return i_sym_max

    def build_symmetry_map(self) -> dict[int, int]:
        """Build the dictionary mapsym (i -> j) that associate one edge j per group of symmetry i."""

        mapsym = {}

        i = 0
        while True:
            j = 0
            has_i_sym = False
            for u, v in self.edges():
                i_sym = self.edge_attribute((u, v), "sym_key")
                if i_sym == i:
                    mapsym[i] = j
                    has_i_sym = True
                    break
                j += 1
            if has_i_sym is False:
                break
            i += 1

        return mapsym

    # =============================================================================
    # Methods - TNA
    # =============================================================================

    def build_dual(self) -> Mesh:
        """
        Build the dual of the form diagram.

        Returns
        -------
        dual
            The middle dual mesh.

        Notes
        -----
        Construction of the dual diagram is based on the faces around the inner, free vertices of the form diagram.
        This means not only the vertices on the boundary are ignored, but also the vertices that are anchored.

        """

        self.modify_diagram_boundary()
        return self.dual_diagram(cls=Mesh)

    def build_complete_dual(self) -> Mesh:
        """
        Build the dual of the form diagram considering half-faces in the edges on the perimeter.

        Returns
        -------
        Mesh
            The middle dual mesh.

        Notes
        -----
        Construction of the dual diagram is based on the faces around all edges.
        The edges in the boundary recieve special tratment and their midpoint is.

        """
        dual = Mesh()
        fkey_centroid = {fkey: self.face_centroid(fkey) for fkey in self.faces()}
        midpt_boundary = {edge: self.edge_midpoint(edge) for edge in self.edges()}
        fixed = self.fixed()
        boundary = list(self.vertices_on_boundary())
        vertices = {}
        faces = {}
        i = 0  # face keys
        j = 0  # vertex keys
        for key in self.vertices():
            facekeys = []
            fkeys = self.vertex_faces(key, ordered=True)
            for fkey in fkeys:
                pt = fkey_centroid[fkey]
                if pt not in vertices.values():
                    vertices[j] = pt
                    facekeys.append(j)
                    j += 1
                else:
                    j_back = list(vertices.keys())[list(vertices.values()).index(pt)]
                    facekeys.append(j_back)
            if key in boundary:
                nbrs = self.vertex_neighbors(key, ordered=True)
                nbrs.reverse()
                for nbr in nbrs:
                    if nbr in boundary:
                        if (key, nbr) in midpt_boundary.keys():
                            pt = midpt_boundary[(key, nbr)]
                        else:
                            pt = midpt_boundary[(nbr, key)]

                        if pt not in vertices.values():
                            vertices[j] = pt
                            facekeys.append(j)
                            j += 1
                        else:
                            j_back = list(vertices.keys())[list(vertices.values()).index(pt)]
                            facekeys.append(j_back)
            if key in fixed:
                pt = self.vertex_coordinates(key)
                vertices[j] = pt
                facekeys = facekeys[:-1] + [j] + [facekeys[-1]]  # add the support as last element
                j += 1

            faces[i] = facekeys
            i += 1

        for key, (x, y, z) in vertices.items():
            dual.add_vertex(key, x=x, y=y, z=z)
        for fkey, vertices in faces.items():
            dual.add_face(vertices, fkey=fkey)

        return dual

    def modify_diagram_boundary(self) -> None:
        """
        Modify the diagram boundary adding the outer perimeter around the fixed vertices with edges having ``_is_edge`` False.

        Parameters
        ----------
        delete_face : bool, optional
            Option to delete the outside face of the diagram.
            The default value is ``False``.

        Returns
        -------
        None
            The ``FormDiagram`` is modified in place.

        """

        corners = list(self.vertices_where(is_fixed=True))
        self.vertices_attributes(names=("is_anchor", "is_fixed"), values=(True, True), keys=corners)
        self.update_boundaries()

    def restore_diagram_boundary(self, tol=1e-3) -> None:
        """
        Remove faces of the diagram with area zero.


        Returns
        -------
        None
            The ``FormDiagram`` is modified in place.

        """

        faces_delete = []
        for face in self.faces():
            A = self.face_area(face)
            if A < tol:
                faces_delete.append(face)

        for face in faces_delete:
            self.delete_face(face)

    # =============================================================================
    # Methods - Shape
    # =============================================================================

    def apply_envelope_from_shape(self, shape: Shape) -> None:
        """Apply an envelope (intrados and extrados) to the FormDiagram based on the input shape.

        Parameters
        ----------
        shape : Shape
            The input Shape with intrados and extrardos.

        Returns
        -------
        None
            The formdiagram is updated in place.

        """
        x = self.vertices_attribute("x")  # check if array is necessary here
        y = self.vertices_attribute("y")

        # it is a bit silly to only use the shape type
        # and not just hand this off to the shape itself...

        if shape.parameters["type"] == "dome":
            zub, zlb = dome_ub_lb_update(
                x,
                y,
                shape.parameters["thk"],
                shape.parameters["t"],
                shape.parameters["center"],
                shape.parameters["radius"],
            )

        elif shape.parameters["type"] == "crossvault":
            zub, zlb = crossvault_ub_lb_update(
                x,
                y,
                shape.parameters["thk"],
                shape.parameters["t"],
                shape.parameters["xy_span"],
            )

        elif shape.parameters["type"] == "pointed_crossvault":
            zub, zlb = pointed_vault_ub_lb_update(
                x,
                y,
                shape.parameters["thk"],
                shape.parameters["t"],
                shape.parameters["xy_span"],
                hc=shape.parameters["hc"],
                he=shape.parameters["he"],
                hm=shape.parameters["hm"],
            )

        elif shape.parameters["type"] == "arch":
            zub, zlb = arch_ub_lb_update(
                x,
                y,
                shape.parameters["thk"],
                shape.parameters["t"],
                H=shape.parameters["H"],
                L=shape.parameters["L"],
                x0=shape.parameters["x0"],
            )

        elif shape.parameters["type"] == "pointed_arch":
            zub, zlb = pointed_arch_ub_lb_update(
                x,
                y,
                shape.parameters["thk"],
                shape.parameters["t"],
                hc=shape.parameters["hc"],
                L=shape.parameters["L"],
                x0=shape.parameters["x0"],
            )

        elif shape.parameters["type"] == "pavillionvault":
            zub, zlb = pavillionvault_ub_lb_update(
                x,
                y,
                shape.parameters["thk"],
                shape.parameters["t"],
                shape.parameters["xy_span"],
                shape.parameters["spr_angle"],
            )

        elif shape.parameters["type"] == "general":
            XY = self.vertices_attributes("xy")
            zub = shape.get_ub_pattern(XY)
            zlb = shape.get_lb_pattern(XY)

        else:
            raise Exception

        for index, vertex in enumerate(self.vertices()):
            self.vertex_attribute(vertex, name="ub", value=float(zub[index]))
            self.vertex_attribute(vertex, name="lb", value=float(zlb[index]))

    def apply_envelope_on_xy(self, c=0.5) -> None:
        """Apply an envelope to the FormDiagram in the plan (x, y) by a given distance.

        Parameters
        ----------
        c : float, optional
            The maximum allowed movement of the nodes in ``x`` or ``y``, by default 0.5

        Returns
        ----------
        None
            The formdiagram is updated in place in the attributes

        """
        for vertex in self.vertices():
            x, y = self.vertex_attributes(vertex, names=["x", "y"])

            self.vertex_attribute(vertex, name="xmin", value=x - c)
            self.vertex_attribute(vertex, name="xmax", value=x + c)
            self.vertex_attribute(vertex, name="ymin", value=y - c)
            self.vertex_attribute(vertex, name="ymax", value=y + c)

    def apply_envelope_on_xy_from_base(self, form_base: "FormDiagram", c=0.5) -> None:
        """Apply an envelope to the FormDiagram considering a given distance applied to a base form diagram.

        Parameters
        ----------
        form_base : FormDiagram
            The base FormDiagram to consider.
        c : float, optional
            The maximum allowed movement of the nodes in ``x`` or ``y``, by default 0.5

        Returns
        -------
        None
            The formdiagram is updated in place in the attributes

        """
        for vertex in form_base.vertices():
            x, y = form_base.vertex_attributes(vertex, names=["x", "y"])

            self.vertex_attribute(vertex, name="xmin", value=x - c)
            self.vertex_attribute(vertex, name="xmax", value=x + c)
            self.vertex_attribute(vertex, name="ymin", value=y - c)
            self.vertex_attribute(vertex, name="ymax", value=y + c)

    def apply_bounds_on_q(self, qmin=-1e4, qmax=1e-8) -> None:
        """Apply bounds on the magnitude of the edges'force densities.

        Parameters
        ----------
        qmin : float, optional
            The minimum allowed force density ``qmin``, by default -1e+4
        qmax : float, optional
            The maximum allowed force density ``qmax``, by default 1e-8

        Returns
        -------
        None
            The formdiagram is updated in place in the attributes.

        """
        if isinstance(qmin, list):
            for i, edge in enumerate(self.edges_where({"_is_edge": True})):
                self.edge_attribute(edge, "qmin", qmin[i])
                self.edge_attribute(edge, "qmax", qmax[i])
        else:
            for i, edge in enumerate(self.edges_where({"_is_edge": True})):
                self.edge_attribute(edge, "qmin", qmin)
                self.edge_attribute(edge, "qmax", qmax)

    def apply_bounds_reactions(self, shape: Shape, assume_shape=None):
        """Apply bounds on the magnitude of the allowed increase in thickness of the upper-bound (tub), lower-bound (tlb), and of the reaction vector (tub_reacmax).

        Parameters
        ----------
        shape : Shape
            The shape of masonry to constraint the form diagram to.
        assume_shape : dict, optional
            Whether or not consider the settings of a different shape on the process, by default None

        Returns
        ----------
        None
            The formdiagram is updated in place in the attributes.
        """
        if assume_shape:
            parameters = assume_shape
        else:
            parameters = shape.parameters

        thk = parameters["thk"]

        if parameters["type"] == "dome" or parameters["type"] == "dome_polar":
            [x0, y0] = parameters["center"][:2]
            for key in self.vertices_where({"is_fixed": True}):
                x, y, _ = self.vertex_coordinates(key)
                theta = math.atan2((y - y0), (x - x0))
                x_ = thk / 2 * math.cos(theta)
                y_ = thk / 2 * math.sin(theta)
                self.vertex_attribute(key, "b", [x_, y_])

        b_manual = parameters.get("b_manual", None)
        if parameters["type"] == "arch":
            H = parameters["H"]
            L = parameters["L"]
            thk = parameters["thk"]
            radius = H / 2 + (L**2 / (8 * H))
            zc = radius - H
            re = radius + thk / 2
            x = math.sqrt(re**2 - zc**2)
            for key in self.vertices_where({"is_fixed": True}):
                self.vertex_attribute(key, "b", [x - L / 2, 0.0])
                if b_manual:
                    self.vertex_attribute(key, "b", [b_manual, 0.0])
                    print("Applied b manual")

        if parameters["type"] == "dome_spr":
            x0 = parameters["center"][0]
            y0 = parameters["center"][1]
            radius = parameters["radius"]
            [_, theta_f] = parameters["theta"]
            r_proj_e = (radius + thk / 2) * math.sin(theta_f)
            r_proj_m = (radius) * math.sin(theta_f)
            delt = r_proj_e - r_proj_m
            for key in self.vertices_where({"is_fixed": True}):
                x, y, _ = self.vertex_coordinates(key)
                theta = math.atan2((y - y0), (x - x0))
                x_ = delt * math.cos(theta)
                y_ = delt * math.sin(theta)
                self.vertex_attribute(key, "b", [x_, y_])

        if parameters["type"] == "pavillionvault":
            x0, x1 = parameters["xy_span"][0]
            y0, y1 = parameters["xy_span"][1]
            for key in self.vertices_where({"is_fixed": True}):
                x, y, _ = self.vertex_coordinates(key)
                if x == x0:
                    self.vertex_attribute(key, "b", [-thk / 2, 0])
                elif x == x1:
                    self.vertex_attribute(key, "b", [+thk / 2, 0])
                if y == y0:
                    if self.vertex_attribute(key, "b"):
                        b = self.vertex_attribute(key, "b")
                        b[1] = -thk / 2
                        self.vertex_attribute(key, "b", b)
                    else:
                        self.vertex_attribute(key, "b", [0, -thk / 2])
                elif y == y1:
                    if self.vertex_attribute(key, "b"):
                        b = self.vertex_attribute(key, "b")
                        b[1] = +thk / 2
                        self.vertex_attribute(key, "b", b)
                    else:
                        self.vertex_attribute(key, "b", [0, +thk / 2])

    def apply_bounds_tub_tlb(self, tubmax=0.5, tlbmax=0.5):
        """Apply bounds on the magnitude of the allowed increase in thickness of the upper-bound (tub), lower-bound (tlb), and of the reaction vector (tub_reacmax).

        Parameters
        ----------
        tubmax : float, optional
            The maximum increase in thickness of the extrados. The default value is ``0.5``.
        tlbmax : float, optional
            The maximum increase in thickness of the intrados. The default value is ``0.5``.

        Returns
        -------
        None
            The formdiagram is updated in place in the attributes.

        """
        for vertex in self.vertices():
            self.vertex_attribute(vertex, "tubmax", tubmax)
            self.vertex_attribute(vertex, "tlbmax", tlbmax)

    # =============================================================================
    # Methods - Loads
    # =============================================================================

    def apply_selfweight_from_shape(self, shape: Shape, pz_negative=True, normalize=True) -> None:
        """Apply selfweight to the nodes of the form diagram based on the shape

        Parameters
        ----------
        shape : Shape
            Shape of the masonry
        pz_negative : bool, optional
            Wether or not the vertical loads are negative, by default True
        normalize : bool, optional
            Wether or not normalise the selfweight to match shape.total_weight, by default True

        Returns
        -------
        None
            The FormDiagram is modified in place

        """
        form_: "FormDiagram" = self.copy()
        total_selfweight = shape.compute_selfweight()
        ro = shape.ro
        thk = shape.parameters["thk"]

        x = self.vertices_attribute("x")
        y = self.vertices_attribute("y")

        if shape.parameters["type"] == "dome":
            zt = dome_zt_update(x, y, shape.parameters["radius"], shape.parameters["t"], shape.parameters["center"])

        elif shape.parameters["type"] == "crossvault":
            zt = crossvault_middle_update(x, y, shape.parameters["t"], xy_span=shape.parameters["xy_span"])

        elif shape.parameters["type"] == "pointed_crossvault":
            zt = pointed_vault_middle_update(
                x,
                y,
                shape.parameters["t"],
                xy_span=shape.parameters["xy_span"],
                hc=shape.parameters["hc"],
                he=shape.parameters["he"],
                hm=shape.parameters["hm"],
            )

        else:
            XY = self.vertices_attributes("xy")
            zt = shape.get_middle_pattern(XY)

        i = 0
        for key in form_.vertices():
            z = float(zt[i])
            form_.vertex_attribute(key, "z", value=z)
            self.vertex_attribute(key, "target", value=z)
            i += 1

        pzt = 0
        for key in self.vertices():
            pz = form_.vertex_area(key)
            self.vertex_attribute(key, "pz", value=pz)
            pzt += pz

        if shape.parameters["type"] == "arch" or shape.parameters["type"] == "pointed_arch":
            pzt = 0
            for key in self.vertices():
                self.vertex_attribute(key, "pz", value=1.0)
                if self.vertex_attribute(key, "is_fixed") is True:
                    self.vertex_attribute(key, "pz", value=0.5)
                pzt += self.vertex_attribute(key, "pz")

        factor = 1.0 * ro * thk  # Transform tributary area in tributary load
        if normalize:
            factor = total_selfweight / pzt
        if pz_negative:
            factor *= -1  # make loads negative

        for key in self.vertices():
            pzi = factor * self.vertex_attribute(key, "pz")
            self.vertex_attribute(key, "pz", value=pzi)

    def apply_selfweight_from_pattern(self, pattern: Mesh, pz_negative=True, tol=10e-4):
        """Apply selfweight to the nodes considering a different Form Diagram to locate loads.
        Note: the base pattern has to coincide with nodes from the original form diagram.

        Parameters
        ----------
        pattern : Mesh
            Mesh in which the forces should be based.
        plot : bool, optional
            Whether or not plot the diagrams and its overlap, by False
        pz_negative : bool, optional
            Wether or not the vertical loads are negative by default ``True``.
        tol : float, optional
            Tolerance for the nodal position match of form and pattern, by default ``10e-4``.

        Returns
        -------
        None
            The FormDiagram is modified in place.
        """
        self.vertices_attribute("pz", 0.0)
        key_real_to_key = {}

        for key in pattern.vertices():
            x, y, _ = pattern.vertex_coordinates(key)
            for key_real in self.vertices():
                x_real, y_real, _ = self.vertex_coordinates(key_real)
                if x - tol < x_real < x + tol and y - tol < y_real < y + tol:
                    key_real_to_key[key_real] = key
                    break

        pzt = 0
        for key in key_real_to_key:
            pz = pattern.vertex_attribute(key_real_to_key[key], "pz")
            # if pz_negative:
            #     pz *= -1  # make loads negative
            self.vertex_attribute(key, "pz", value=pz)
            pzt += pz

    def apply_selfweight_from_thrust(self, thickness=0.5, density=20.0):
        """Lump the selfweight in the nodes of the thrust network based on their current position.
        The loads are computed based on the tributary area times the thickness times the density.
        For variable thickness the nodal attribute `thk` is considered.

        Parameters
        ----------
        thickness : float, optional
            The thickness of the problem, by default 0.50
            If None is passed, the thickness is taken from the nodal attribute `thk`
        density : float, optional
            The density of the material, by default 20.0

        Return
        ------
        None
            The form diagram is updated in place

        """
        for key in self.vertices():
            ai = self.vertex_area(key)
            if thickness:
                load = -1 * ai * thickness * density
            else:
                thk = self.vertex_attribute(key, "thk")
                load = -1 * ai * thk * density
            self.vertex_attribute(key, "pz", load)

    def apply_horizontal_multiplier(self, lambd=1.0, direction="x"):
        """Modify the applied loads considering a load multiplier.

        Parameters
        ----------
        lambd : float, optional
            Value of the horizontal multiplier, by default ``1.0``.
        direction : str, optional
            Direction to apply the loads, ``x`` or ``y``, by default ``x``.

        Returns
        -------
        None
            The FormDiagram is modified in place.

        """
        arg = "p" + direction

        for key in self.vertices():
            pz = self.vertex_attribute(key, "pz")
            self.vertex_attribute(key, arg, -1 * pz * lambd)  # considers that swt (pz) is negative

    def apply_fill_load(self):
        """Modify the applied loads considering a fill."""

        raise NotImplementedError("Not implemented")
