from compas.utilities import geometric_key

from compas.datastructures import Mesh
from compas.datastructures import mesh_bounding_box_xy

from compas_tno.diagrams.diagram_arch import create_arch_form_diagram
from compas_tno.diagrams.diagram_arch import create_linear_form_diagram
from compas_tno.diagrams.diagram_rectangular import create_cross_form
from compas_tno.diagrams.diagram_rectangular import create_cross_diagonal
from compas_tno.diagrams.diagram_rectangular import create_cross_with_diagonal
from compas_tno.diagrams.diagram_rectangular import create_fan_form
from compas_tno.diagrams.diagram_circular import create_circular_radial_form
from compas_tno.diagrams.diagram_circular import create_circular_radial_spaced_form
from compas_tno.diagrams.diagram_circular import create_circular_spiral_form
from compas_tno.diagrams.diagram_rectangular import create_ortho_form

from compas_tna.diagrams import FormDiagram

import math


__all__ = ['FormDiagram']


class FormDiagram(FormDiagram):

    """The ``FormDiagram`` class imports the attributes and set ups from ``compas_tna.diagrams.FormDiagram`` and include some functionalities
    useful for the assessment of masonry structures.

    It defined the form-diagram that will be the layout of the forces within the structure

    Notes
    -----
    A ``FormDiagram`` has the following constructor functions

    *   ``from_library`` : Construct the shape from a dictionary with instructions, the library supports the creation of parametric arches, domes, and vaults.
    *   ``from_rhinomesh`` : Construct Extrados, Intrados and Middle surfaces from RhinoMeshes.
    *   ``from_rhinosurface`` : Construct Extrados, Intrados and Middle surfaces using the U and V isolines.

    A parametric ``FormDiagram`` contains the following information:

    *   ``data``

        *   ``type``  : The type of the Form Diafram to be constructed.
        *   ``xy_span``  : Planar range of the structure.
        *   ``radius``  : In case of a circular diagram.
        *   ``center``  : In case of a circular diagram.
        *   ``density``  : Density of the diagram.

    """

    # __module__ = 'compas_tno.diagrams'

    def __init__(self):
        super(FormDiagram, self).__init__()
        self.update_default_vertex_attributes({
            'is_roller': False,
        })
        self.update_default_edge_attributes({
            'q': 1.0,
            'f': 1.0,
            'is_symmetry': False,
            'is_ind': False,
            'qmin': -1e+4,
            'qmax': 1e-8
        })
        self.args = None
        self.attributes['loadpath'] = 0.0
        self.attributes['indset'] = None
        self.parameters = {
            'type': None,
            'discretisation': None,
            # 'x0': None,
            # 'total_nodes': 20,
            # 'xy_span': [[0,10],[0,10]],
            # 'discretisation': [10,10],
            # 'fix': 'corners',
            # 'r_oculus:' None,
        }

    # --------------------------------------------------------------- #
    # -----------------------CONSTRUCTORS---------------------------- #
    # --------------------------------------------------------------- #

    @classmethod
    def from_library(cls, data):
        """ Construct a FormDiagram the library of FormDiagrams available.

        Parameters
        ----------
        data : dict
            Dictionary with information to generate the FormDiagram.

        Returns
        -------
        obj
            FormDiagram.
        """

        form_type = data['type']

        r_oculus = data.get('r_oculus', 0.0)
        diagonal = data.get('diagonal')
        partial_diagonal = data.get('partial_diagonal')
        partial_bracing_modules = data.get('partial_bracing_modules')

        if form_type == 'arch':
            form = cls().create_arch_form_diagram(H=data['H'], L=data['L'], x0=data['x0'], total_nodes=data['total_nodes'])
        if form_type == 'pointed_arch':
            form = cls().create_pointed_arch(L=data['L'], x0=data['x0'], total_nodes=data['total_nodes'])
        if form_type == 'cross_fd':
            form = cls().create_cross_form(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data['fix'])
        if form_type == 'cross_diagonal':
            form = cls().create_cross_diagonal(xy_span=data['xy_span'], discretisation=data['discretisation'], partial_bracing_modules=partial_bracing_modules, fix=data['fix'])
        if form_type == 'cross_with_diagonal':
            form = cls().create_cross_with_diagonal(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data.get('fix', 'all'))
        if form_type == 'fan_fd':
            form = cls().create_fan_form(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data['fix'])
        if form_type == 'ortho' or form_type == 'ortho_fd':
            form = cls().create_ortho_form(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data['fix'])
        if form_type == 'radial_fd':
            form = cls().create_circular_radial_form(center=data['center'], radius=data['radius'],
                                                     discretisation=data['discretisation'], r_oculus=r_oculus, diagonal=diagonal, partial_diagonal=partial_diagonal)
        if form_type == 'radial_spaced_fd':
            form = cls().create_circular_radial_spaced_form(center=data['center'], radius=data['radius'],
                                                            discretisation=data['discretisation'], r_oculus=r_oculus, diagonal=diagonal, partial_diagonal=partial_diagonal)
        if form_type == 'spiral_fd':
            form = cls().create_circular_spiral_form(center=data['center'], radius=data['radius'], discretisation=data['discretisation'], r_oculus=r_oculus)

        form.parameters = data

        return form

    @classmethod
    def create_arch(cls, H=1.00, L=2.00, x0=0.0, total_nodes=100):
        """ Construct a FormDiagram based on an arch linear discretisation.

        Parameters
        ----------
        D : float
            Central diameter of the arch.

        x0: float
            Beginning of the linear form diagram.

        total_nodes : int
            Numbers of nodes to be considered in the form diagram.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_arch_form_diagram(cls(), L=L, H=H, x0=x0, total_nodes=total_nodes)

        return form

    @classmethod
    def create_pointed_arch(cls, L=2.00, x0=0.0, total_nodes=100):

        form = create_linear_form_diagram(cls(), L=L, x0=x0, total_nodes=total_nodes)

        return form

    @classmethod
    def create_cross_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=10, fix='corners'):
        """ Construct a FormDiagram based on cross discretiastion with orthogonal arrangement and diagonal.

        Parameters
        ----------
        xy_span : list
            List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

        discretisation: list
            Set the density of the grid in x and y directions.

        fix : string
            Option to select the constrained nodes: 'corners', 'all' are accepted.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_cross_form(cls(), xy_span=xy_span, discretisation=discretisation, fix=fix)

        return form

    @classmethod
    def create_cross_diagonal(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=10, partial_bracing_modules=None, fix='corners'):
        """ Construct a FormDiagram based on a mixture of cross and fan discretiastion.

        Parameters
        ----------
        xy_span : list
            List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

        discretisation: list
            Set the density of the grid in x and y directions.

        fix : string
            Option to select the constrained nodes: 'corners', 'all' are accepted.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_cross_diagonal(cls(), xy_span=xy_span, discretisation=discretisation, partial_bracing_modules=partial_bracing_modules, fix=fix)

        return form

    @classmethod
    def create_cross_with_diagonal(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=10, fix='all'):
        """ Construct a FormDiagram based on cross discretiastion with diagonals.

        Parameters
        ----------
        xy_span : list
            List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

        discretisation: list
            Set the density of the grid in x and y directions.

        fix : string
            Option to select the constrained nodes: 'corners', 'all' are accepted.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_cross_with_diagonal(cls(), xy_span=xy_span, discretisation=discretisation, fix=fix)

        return form

    @classmethod
    def create_fan_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10, 10], fix='corners'):
        """ Helper to construct a FormDiagram based on fan discretiastion with straight lines to the corners.

        Parameters
        ----------
        xy_span : list
            List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

        division: int
            Set the density of the grid in x and y directions.

        fix : string
            Option to select the constrained nodes: 'corners', 'all' are accepted.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_fan_form(cls(), xy_span=xy_span, discretisation=discretisation, fix=fix)

        return form

    @classmethod
    def create_ortho_form(cls, xy_span=[[0.0, 10.0], [0.0, 10.0]], discretisation=[10, 10], fix='corners'):
        """ Helper to construct a FormDiagram based on a simple orthogonal discretisation.

        Parameters
        ----------
        xy_span : list
            List with initial- and end-points of the vault [(x0,x1),(y0,y1)].

        division: int
            Set the density of the grid in x and y directions.

        fix : string
            Option to select the constrained nodes: 'corners', 'all' are accepted.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_ortho_form(cls(), xy_span=xy_span, discretisation=discretisation, fix=fix)

        return form

    @classmethod
    def create_circular_radial_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0, diagonal=False, partial_diagonal=False):
        """ Construct a circular radial FormDiagram with hoops not equally spaced in plan.

        Parameters
        ----------
        center : list
            Planar coordinates of the form-diagram [xc, yc].

        radius: float
            Radius of the form-diagram

        discretisation : list
            Number of meridians, and of spikes from the center on the dome form-diagram.

        r_oculus: float
            Value of the radius of the oculus, if no oculus is present should be set to zero.

        diagonal: float
            Activate diagonal in the quads.

        partial_diagonal: float
            Activate partial diagonal in the quads.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_circular_radial_form(cls(), center=center, radius=radius, discretisation=discretisation,
                                           r_oculus=r_oculus, diagonal=diagonal, partial_diagonal=partial_diagonal)

        return form

    @classmethod
    def create_circular_radial_spaced_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0, diagonal=False, partial_diagonal=False):
        """ Construct a circular radial FormDiagram with hoops not equally spaced in plan, but equally spaced with regards to the projection on a hemisphere.

        Parameters
        ----------
        center : list
            Planar coordinates of the form-diagram [xc, yc].
        radius: float
            Radius of the form-diagram
        discretisation : list
            Number of meridians, and of spikes from the center on the dome form-diagram.
        r_oculus: float
            Value of the radius of the oculus, if no oculus is present should be set to zero.
        diagonal: float
            Activate diagonal in the quads.
        partial_diagonal: float
            Activate partial diagonal in the quads.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_circular_radial_spaced_form(cls(), center=center, radius=radius, discretisation=discretisation,
                                                  r_oculus=r_oculus, diagonal=diagonal, partial_diagonal=partial_diagonal)

        return form

    @classmethod
    def create_circular_spiral_form(cls, center=[5.0, 5.0], radius=5.0, discretisation=[8, 20], r_oculus=0.0):
        """ Construct a circular radial FormDiagram with hoops not equally spaced in plan, but equally spaced with regards to the projection on a hemisphere.

        Parameters
        ----------
        center : list
            Planar coordinates of the form-diagram [xc, yc].

        radius: float
            Radius of the form-diagram

        discretisation : list
            Number of meridians, and of spikes from the center on the dome form-diagram.

        r_oculus: float
            Value of the radius of the oculus, if no oculus is present should be set to zero.

        diagonal: float
            Activate diagonal in the quads.

        partial_diagonal: float
            Activate partial diagonal in the quads.

        Returns
        -------
        obj
            FormDiagram.

        """

        form = create_circular_spiral_form(cls(), center=center, radius=radius, discretisation=discretisation, r_oculus=r_oculus)

        return form

    @classmethod
    def from_assembly(self):
        NotImplementedError

    @classmethod
    def from_singular(cls, json_path):
        """ W.I.P create form diagram from singular dense mesh.
        """

        from compas.datastructures import mesh_delete_duplicate_vertices
        mesh = Mesh.from_json(json_path)
        mesh_delete_duplicate_vertices(mesh)
        form = cls().from_mesh(mesh)
        form.set_boundary_supports()
        form.delete_boundary_edges()

        return form

    @classmethod
    def from_triangle(cls, boundary_points, area=0.5, angle=30):
        """ W.I.P create form diagram from triangl dense mesh.
        """
        from compas_triangle.delaunay import conforming_delaunay_triangulation
        vertices, faces = conforming_delaunay_triangulation(boundary_points + boundary_points[:1], angle=angle, area=area)
        mesh = Mesh.from_vertices_and_faces(vertices, faces)
        form = cls().from_mesh(mesh)
        return form

    @classmethod
    def from_skeleton(cls, data):
        """ W.I.P create form diagram from skeleton.
        """

        from compas_skeleton.datastructure import Skeleton
        from compas.geometry import matrix_from_scale_factors
        from compas.datastructures import mesh_transform

        skeleton_type = data['type']

        if skeleton_type == 'circular':

            point = data.get('center', [0, 0, 0])
            subdivisions = data.get('subdivision', 3)
            radius = data.get('radius', 5.0)

            skeleton = Skeleton.from_center_point(point)
            skeleton.node_width = 7.5
            skeleton.update_mesh_vertices_pos()
            skeleton.vertex_attribute(0, 'z', 10.0)
            skeleton.subdivide(subdivisions)
            skeleton_mesh = skeleton.to_mesh()

            max_x = 0
            for key in skeleton_mesh.vertices():
                x, _, _ = skeleton_mesh.vertex_coordinates(key)
                if x > max_x:
                    max_x = x
            scale = radius / max_x
            T = matrix_from_scale_factors([scale, scale, 0.0])
            mesh_transform(skeleton_mesh, T)
            form = cls().from_mesh(skeleton_mesh)
            form.set_boundary_supports()
            form.delete_boundary_edges()

        return form

    # --------------------------------------------------------------- #
    # -----------------------ULTILITIES------------------------------ #
    # --------------------------------------------------------------- #

    def q(self):
        """Return the force densities ``q`` of the diagram."""
        return [self.edge_attribute(edge, 'q') for edge in self.edges_where({'_is_edge': True})]

    def xy(self):
        """Return the ``xy`` coordinates of the diagram."""
        return self.vertices_attributes('xy')

    def fixed(self):
        """Return the keys of the fixed vertices of the diagram."""
        return list(self.vertices_where({'is_fixed': True}))

    def fixed_x(self):
        """Return the keys of the fixed-in-x vertices of the diagram."""
        return list(self.vertices_where({'is_fixed_x': True, 'is_fixed': False}))

    def fixed_y(self):
        """Return the keys of the fixed-in-y vertices of the diagram."""
        return list(self.vertices_where({'is_fixed_y': True, 'is_fixed': False}))

    def ind(self):
        """Return the identifiers for the edges selected as independents on the diagram."""
        return list(self.edges_where({'is_ind': True}))

    def number_of_real_edges(self):
        """Return the identifiers for the edges selected as independents on the diagram."""
        return len(list(self.edges_where({'_is_edge': True})))

    def update_f(self):
        """Update attributes force ``f`` based on force densities ``q``and lengths ``l``."""

        for u, v in self.edges():
            f = self.edge_attribute((u, v), 'q') * self.edge_length(u, v)
            self.edge_attribute((u, v), 'f', f)

        return

    def overview_forces(self):
        """Print an overview of the forces in the ``Form Diagram``."""

        self.update_f()

        f = []
        q = []
        z = []
        pz = 0

        lp = 0

        for u, v in self.edges_where({'_is_edge': True}):
            qi = self.edge_attribute((u, v), 'q')
            li = self.edge_length(u, v)
            lp += qi*li**2
            q.append(qi)
            f.append(qi*li)

        print('='*20)
        print('Overview on forces:')

        print('q: {0:.3f} : {1:.3f}'.format(float(min(q)), float(max(q))))
        print('f: {0:.3f} : {1:.3f}'.format(float(min(f)), float(max(f))))
        for key in self.vertices():
            z.append(self.vertex_attribute(key, 'z'))
            pz += self.vertex_attribute(key, 'pz')
        print('z: {0:.3f} : {1:.3f}'.format(float(min(z)), float(max(z))))
        print('pz: {0:.3f}'.format(pz))
        print('lp: {0:.3f}'.format(lp))
        self.attributes['loadpath'] = lp

        return

    def thrust(self):
        """Returns the (horizontal) thrust."""

        thrust = 0
        for key in self.vertices_where({'is_fixed': True}):
            rx = self.vertex_attribute(key, '_rx')
            ry = self.vertex_attribute(key, '_ry')
            r = math.sqrt(rx**2 + ry**2)
            thrust += r

        return thrust

    def lumped_swt(self):
        """Returns the lumped selfweight on the nodes."""

        swt = 0
        for key in self.vertices():
            pz = self.vertex_attribute(key, 'pz')
            swt += pz

        return swt

    def distance_target(self):
        """ Compute the squared vertical distance from a target."""

        f = 0
        for key in self.vertices():
            z = self.vertex_attribute(key, 'z')
            s = self.vertex_attribute(key, 'target')
            f += (z - s)**2

        return f

    def loadpath(self):
        """ Compute loadpath in the current configuration based on the force attributes stored in the ``FormDiagram``. """

        lp = 0
        for u, v in self.edges_where({'_is_external': False}):
            if self.edge_attribute((u, v), '_is_edge') is True and self.edge_attribute((u, v), 'is_symmetry') is False:
                qi = self.edge_attribute((u, v), 'q')
                li = self.edge_length(u, v)
                lp += qi*li**2

        self.attributes['loadpath'] = lp

        return lp

    def number_of_independents(self):
        """ Compute the number of independent edges."""
        return len(list(self.vertices_where({'is_ind': True})))

    # --------------------------------------------------------------- #
    # -----------------------BOUNDARIES------------------------------ #
    # --------------------------------------------------------------- #

    def number_of_supports(self, printout=False):
        """ Compute the number of supports."""
        return len(list(self.vertices_where({'is_fixed': True})))

    def delete_boundary_edges(self, delete_corner_vertex=True):  # WIP: Fix this for a pavillion vault geometry with 'cross_fd' form diagram.
        """
        Delete boundary edges on a diagram.

        Parameters
        ----------
        delete_corner_vertex : bool, optional
            Delete the corner vertices.
            The default value is ``True``.

        Returns
        -------
        None
            The ``FormDiagram`` is modified in place.

        """

        for u, v in self.edges_on_boundary():
            self.edge_attribute((u, v), '_is_edge', False)

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
                self.edge_attribute((neighbors[0], neighbors[1]), '_is_edge', False)

        return self

    def set_boundary_supports(self):
        """ Set all node on the boundary of a pattern as supports."""

        for key in self.vertices_on_boundary():
            self.vertex_attribute(key, 'is_fixed', True)

        return

    def set_boundary_rollers(self, max_rx=[0.0, 0.0], max_ry=[0.0, 0.0], total_rx=None, total_ry=None):
        """
        Set boundary vertices as rollers.

        Parameters
        ----------
        max_rx : list, optional
            Assign the maximum reaction in x to the vertices in the left and right boundary.
            The default value is ``[0.0, 0.0]``.
        max_ry : list, optional
            Assign the maximum reaction in y to the vertices in the top and bottom boundary.
            The default value is ``[0.0, 0.0]``.
        total_rx : list, optional
            Assign the maximum total reaction in x of the vertices in the right or left boundary.
            The default value is ``None``, in which the ``max_rx`` is used for all vertices.
        total_ry : list, optional
            Assign the maximum total reaction in y of the vertices in the top or bottom boundary.
            The default value is ``None``, in which the ``max_ry`` is used for all vertices.

        Returns
        -------
        None
            The ``FormDiagram`` is modified in place.

        """
        corners = mesh_bounding_box_xy(self)
        xs = [point[0] for point in corners]
        ys = [point[1] for point in corners]
        xlimits = [min(xs), max(xs)]
        ylimits = [min(ys), max(ys)]
        if not total_rx:
            nx0 = 0
            nx1 = 0
            for key in self.vertices_where({'x': xlimits[0]}):
                nx0 += 1
            for key in self.vertices_where({'x': xlimits[1]}):
                nx1 += 1
            max_rx = [total_rx[0]/nx0, total_rx[1]/nx1]
        if not total_ry:
            ny0 = 0
            ny1 = 0
            for key in self.vertices_where({'y': ylimits[0]}):
                ny0 += 1
            for key in self.vertices_where({'y': ylimits[1]}):
                ny1 += 1
            max_ry = [total_ry[0]/ny0, total_ry[1]/ny1]
        for key in self.vertices_on_boundary():
            if self.vertex_attribute(key, 'is_fixed') is False:
                x, y, _ = self.vertex_coordinates(key)
                if x == xlimits[0]:
                    self.vertex_attribute(key, 'rol_x', True)
                    self.vertex_attribute(key, 'max_rx', max_rx[0])
                elif x == xlimits[1]:
                    self.vertex_attribute(key, 'rol_x', True)
                    self.vertex_attribute(key, 'max_rx', max_rx[1])
                elif y == ylimits[0]:
                    self.vertex_attribute(key, 'rol_y', True)
                    self.vertex_attribute(key, 'max_ry', max_ry[0])
                elif y == ylimits[1]:
                    self.vertex_attribute(key, 'rol_y', True)
                    self.vertex_attribute(key, 'max_ry', max_ry[1])

        return

    # --------------------------------------------------------------- #
    # --------------------------SYMMETRY----------------------------- #
    # --------------------------------------------------------------- #

    def number_of_sym_edges(self, printout=False):
        """ Compute the symmetric edges."""

        i_sym_max = 0
        for u, v in self.edges_where({'_is_edge': True}):
            try:
                i_sym = self.edge_attribute((u, v), 'sym_key')
            except BaseException:
                print('Warning, no symmetry relation found!')
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        if printout:
            print('Form has {0} unique sym edges'.format(i_sym_max))

        return i_sym_max

    def number_of_sym_vertices(self, printout=False):
        """ Compute the symmetric vertices."""

        i_sym_max = 0
        for key in self.vertices():
            try:
                i_sym = self.vertex_attribute(key, 'sym_key')
            except BaseException:
                print('Warning, no symmetry relation found!')
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        if printout:
            print('Form has {0} sym vertices'.format(i_sym_max))

        return i_sym_max

    def number_of_sym_supports(self, printout=False):
        """ Compute the symmetric supports."""

        i_sym_max = 0
        for key in self.vertices_where({'is_fixed': True}):
            try:
                i_sym = self.vertex_attribute(key, 'sym_key')
            except BaseException:
                print('Warning, no symmetry relation found!')
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        if printout:
            print('Form has {0} unique supports'.format(i_sym_max))

        return i_sym_max

    def build_symmetry_map(self):
        """Build the dictionary mapsym (i -> j) that associate one edge j per group of symmetry i."""

        mapsym = {}

        i = 0
        while True:
            j = 0
            has_i_sym = False
            for u, v in self.edges():
                i_sym = self.edge_attribute((u, v), 'sym_key')
                if i_sym == i:
                    mapsym[i] = j
                    has_i_sym = True
                    break
                j += 1
            if has_i_sym is False:
                break
            i += 1

        return mapsym

    # --------------------------------------------------------------- #
    # -----------------------TNA-CONECTION--------------------------- #
    # --------------------------------------------------------------- #

    def add_feet_(self, delete_face=False):
        """
        Add feet to the support as in ``compas_tna``.

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

        if delete_face:
            self.delete_face(0)
        corners = list(self.vertices_where({'is_fixed': True}))
        self.vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
        self.update_boundaries()
        # self.update_boundaries(self, feet=2)

        return

    def remove_feet(self, openings=None):
        """
        Remove the feet of the support as in ``compas_tna``.

        Parameters
        ----------
        openings : bool, optional
            Whether openings are present.
            The default value is ``None``.

        Returns
        -------
        form_: FormDiagram
            The new ``FormDiagram``.

        """

        lines = []
        qs = {}

        for u, v in self.edges_where({'_is_edge': True, '_is_external': False}):
            s = self.vertex_coordinates(u)
            e = self.vertex_coordinates(v)
            lines.append([s, e])
            qs[geometric_key(self.edge_midpoint(u, v))] = self.edge_attribute((u, v), 'q')

        fixed = [geometric_key(self.vertex_coordinates(key)) for key in self.vertices_where({'is_anchor': True})]
        zs = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_coordinates(key)[2] for key in self.vertices_where({'_is_external': False})}
        pz = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'pz') for key in self.vertices()}
        target = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'target') for key in self.vertices()}
        lb = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'lb') for key in self.vertices()}
        ub = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'ub') for key in self.vertices()}

        form_ = FormDiagram.from_lines(lines)
        form_.update_default_edge_attributes({'q': 1, 'is_symmetry': False, '_is_edge': True})
        form_.update_default_vertex_attributes({'is_roller': False})

        if openings:
            for key in self.faces():
                if self.face_area(key) > openings - 1.0 and self.face_area(key) < openings + 1.0:
                    self.delete_face(key)
                    print('Deleted area of face {0}'.format(key))
                    break
        gkey_key = form_.gkey_key()

        for pt in fixed:
            form_.vertex_attribute(gkey_key[pt], name='is_fixed', value=True)

        for key, attr in form_.vertices(True):
            pzi = pz[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
            zi = zs[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
            ti = target[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
            ub_i = ub[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
            lb_i = lb[geometric_key(form_.vertex_coordinates(key)[:2] + [0])]
            attr['pz'] = pzi
            attr['z'] = zi
            attr['target'] = ti
            attr['lb'] = lb_i
            attr['ub'] = ub_i

        for u, v in form_.edges():
            qi = qs[geometric_key(form_.edge_midpoint(u, v))]
            form_.edge_attribute((u, v), name='q', value=qi)

        return form_
