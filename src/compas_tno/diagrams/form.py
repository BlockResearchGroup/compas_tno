
from compas.utilities import geometric_key
from compas.numerical import connectivity_matrix

from random import shuffle
from numpy import max
from numpy import min
from numpy import zeros
from numpy import newaxis
from numpy import array

from numpy import vstack
from numpy import hstack

import math

import matplotlib.pyplot as plt

from compas_tno.diagrams.diagram_arch import create_arch
from compas_tno.diagrams.diagram_rectangular import create_cross_form
from compas_tno.diagrams.diagram_rectangular import create_cross_diagonal
from compas_tno.diagrams.diagram_rectangular import create_cross_with_diagonal
from compas_tno.diagrams.diagram_rectangular import create_fan_form
from compas_tno.diagrams.diagram_circular import create_circular_radial_form
from compas_tno.diagrams.diagram_circular import create_circular_radial_spaced_form
from compas_tno.diagrams.diagram_circular import create_circular_spiral_form
from compas_tno.diagrams.diagram_rectangular import create_ortho_form

from compas_tno.shapes.dome import dome_zt_update
from compas_tno.shapes.dome import dome_ub_lb_update
from compas_tno.shapes.crossvault import crossvault_middle_update
from compas_tno.shapes.crossvault import crossvault_ub_lb_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_middle_update
from compas_tno.shapes.pointed_crossvault import pointed_vault_ub_lb_update

from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_force

from compas_tno.diagrams import ForceDiagram
from compas_tna.diagrams import FormDiagram

from compas.datastructures import Mesh
from compas.datastructures import mesh_bounding_box_xy
from compas.geometry import distance_point_point_xy
from compas.geometry import distance_point_line_xy

from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax

from compas.geometry import distance_point_point_xy

from compas_plotters import MeshPlotter

__all__ = [
    'FormDiagram'
]


class FormDiagram(FormDiagram):

    """The ``FormDiagram`` class imports the attributes and set ups from compas_tna.diagrams.FormDiagram and include some functionalities useful for the assessment of masonry structures.

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

    __module__ = 'compas_tna.diagrams'

    def __init__(self):
        super(FormDiagram, self).__init__()
        self.update_default_vertex_attributes({
            'is_roller': False,
        })
        self.update_default_edge_attributes({
            'q': 1.0,
            'is_symmetry': False,
            'is_ind': False,
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

        form_type = data['type']

        r_oculus = data.get('r_oculus', 0.0)
        diagonal = data.get('diagonal')
        partial_diagonal = data.get('partial_diagonal')
        partial_bracing_modules = data.get('partial_bracing_modules')

        if form_type == 'arch':
            form = cls().create_arch(H=data['H'], L=data['L'], x0=data['x0'], total_nodes=data['total_nodes'])
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

        form = create_arch(cls(), L=L, H=H, x0=x0, total_nodes=total_nodes)

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
    def from_assembly():
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
                x, y, z = skeleton_mesh.vertex_coordinates(key)
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

    def update_f(self):
        """ Update 'f' attribute.
        """

        for u, v in self.edges():
            f = self.edge_attribute((u, v), 'q') * self.edge_length(u, v)
            self.edge_attribute((u, v), 'f', f)

        return


    def overview_forces(self):
        """ Quick overview of the important parameters of a network.
        """

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

        thrust = 0
        for key in self.vertices_where({'is_fixed': True}):
            rx = self.vertex_attribute(key, '_rx')
            ry = self.vertex_attribute(key, '_ry')
            r = math.sqrt(rx**2 + ry**2)
            thrust += r

        return thrust

    def lumped_swt(self):

        swt = 0
        for key in self.vertices():
            pz = self.vertex_attribute(key, 'pz')
            swt += pz

        return swt

    def shuffle_diagram(self, keep_q=False):
        """ Modify the FormDiagram by shuffling the edges.

        Parameters
        ----------
        form : obj
            Original FormDiagram.

        Returns
        -------
        obj
            Shuffled FormDiagram.

        """

        # Edges

        edges = [self.edge_coordinates(u, v) for u, v in self.edges()]
        edges = [[sp[:2] + [0], ep[:2] + [0]] for sp, ep in edges]
        qs = {geometric_key(self.edge_midpoint(u, v)[:2] + [0]): self.edge_attribute((u, v), 'q') for u, v in self.edges()}
        shuffle(edges)

        form_ = FormDiagram.from_lines(edges, delete_boundary_face=True)
        form_.update_default_edge_attributes({'is_symmetry': False})
        sym = [geometric_key(self.edge_midpoint(u, v)[:2] + [0])for u, v in self.edges_where({'is_symmetry': True})]
        for u, v in form_.edges():
            if geometric_key(form_.edge_midpoint(u, v)) in sym:
                form_.edge_attribute((u, v), 'is_symmetry', True)
            if keep_q:
                form_.edge_attribute((u, v), 'q', qs[geometric_key(form_.edge_midpoint(u, v)[:2] + [0])])

        # Vertices

        gkey_key = form_.gkey_key()
        for key, vertex in self.vertex.items():
            gkey = geometric_key(self.vertex_coordinates(key)[:2] + [0])
            form_.vertex[gkey_key[gkey]] = vertex

        form_.attributes['indset'] = []

        return form_

    def distance_target(self):
        """ Compute the squared distance from a target.
        """

        f = 0
        for key, vertex in self.vertex.items():
            # if vertex.get('_is_external') == False:
            z = vertex.get('z')
            s = vertex.get('target')
            w = vertex.get('weight', 1.0)
            f += w * (z - s)**2

        return f

    def loadpath(self):
        """ Compute loadpath in the current configuration.
        """

        lp = 0
        for u, v in self.edges_where({'_is_external': False}):
            if self.edge_attribute((u, v), '_is_edge') is True and self.edge_attribute((u, v), 'is_symmetry') is False:
                qi = self.edge_attribute((u, v), 'q')
                li = self.edge_length(u, v)
                lp += qi*li**2

        self.attributes['loadpath'] = lp

        return lp

    def delete_boundary_edges(self, delete_corner_vertex=True):
        """ Delete boundary edges on a diagram.
        """

        for u, v in self.edges_on_boundary():
            self.edge_attribute((u, v), '_is_edge', False)

        if delete_corner_vertex:
            corners = [key for key in self.corners()]
            for key in corners:
                self.delete_vertex(key)

        return self

    def set_boundary_supports(self):
        """ Set all node on the boundary of a pattern as supports.
        """

        for key in self.vertices_on_boundary():
            self.vertex_attribute(key, 'is_fixed', True)

        return self

    def set_boundary_rollers(self, max_rx=[0.0, 0.0], max_ry=[0.0, 0.0], total_rx=None, total_ry=None):
        corners = mesh_bounding_box_xy(self)
        xs = [point[0] for point in corners]
        ys = [point[1] for point in corners]
        xlimits = [min(xs), max(xs)]
        ylimits = [min(ys), max(ys)]
        if total_rx is not None:
            nx0 = 0
            nx1 = 0
            for key in self.vertices_where({'x': xlimits[0]}):
                nx0 += 1
            for key in self.vertices_where({'x': xlimits[1]}):
                nx1 += 1
            max_rx = [total_rx[0]/nx0, total_rx[1]/nx1]
        if total_ry is not None:
            ny0 = 0
            ny1 = 0
            for key in self.vertices_where({'y': ylimits[0]}):
                ny0 += 1
            for key in self.vertices_where({'y': ylimits[1]}):
                ny1 += 1
            max_ry = [total_ry[0]/ny0, total_ry[1]/ny1]
        for key in self.vertices_on_boundary():
            if self.vertex_attribute(key, 'is_fixed') is False:
                x, y, z = self.vertex_coordinates(key)
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

    def residual(self, plot=False):
        """ Compute residual forces.
        """

        # Mapping

        k_i = self.key_index()

        # Vertices and edges

        n = self.number_of_vertices()
        fixed = [k_i[key] for key in self.fixed()]
        rol = [k_i[key] for key in self.vertices_where({'is_roller': True})]
        edges = [(k_i[u], k_i[v]) for u, v in self.edges()]
        free = list(set(range(n)) - set(fixed) - set(rol))

        # Co-ordinates and loads

        xyz = zeros((n, 3))
        px = zeros((n, 1))
        py = zeros((n, 1))
        pz = zeros((n, 1))

        for key, vertex in self.vertex.items():
            i = k_i[key]
            xyz[i, :] = self.vertex_coordinates(key)
            px[i] = vertex.get('px', 0)
            py[i] = vertex.get('py', 0)
            pz[i] = vertex.get('pz', 0)

        px = px[free]
        py = py[free]
        pz = pz[free]

        # C and E matrices

        C = connectivity_matrix(edges, 'csr')
        Ci = C[:, free]
        Cit = Ci.transpose()
        uvw = C.dot(xyz)
        U = uvw[:, 0]
        V = uvw[:, 1]
        q = array([self.edge_attribute((u, v), 'q') for u, v in self.edges_where({'_is_edge': True})])[:, newaxis]

        # Horizontal checks

        Rx = Cit.dot(U * q.ravel()) - px.ravel()
        Ry = Cit.dot(V * q.ravel()) - py.ravel()
        R = math.sqrt(Rx**2 + Ry**2)
        Rmax = max(R)
        Rs = sum(R)

        if plot:
            print('Residual Total: {0}'.format(Rs))
            print('Residual Max: {0}'.format(Rmax))

        return Rmax

    def evaluate_a(self, plot=False):
        """ Evaluate angle deviations in a form diagram.
        """

        a_total = 0
        a_max = 0
        for u, v in self.edges_where({'_is_edge': True}):
            a = self.edge_attribute((u, v), 'a')
            a_total += a
            l = self.edge_length(u, v)
            a = a*l
            if a > a_max:
                a_max = a
        if plot is True:
            print('Angle Deviation  Max: {0}'.format(a_max))
            print('Angle Deviation  Total: {0}'.format(a_total))

        return a_max

    def optimise_loadpath(self, find_inds=True, qmax=3000, printout=False):

        from compas_tno.solvers.solver_MATLAB import run_loadpath_from_form_MATLAB

        self = run_loadpath_from_form_MATLAB(self, find_inds=find_inds, qmax=qmax, printout=printout)

        return self

    def initialise_loadpath(self, find_inds=True, qmax=10000, printout=False):

        self = self.optimise_loadpath(find_inds=find_inds, qmax=qmax, printout=printout)

        return self

    def envelope_from_shape(self, shape):

        XY = array(self.vertices_attributes('xy'))

        if shape.data['type'] == 'dome':
            zub, zlb = dome_ub_lb_update(XY[:, 0], XY[:, 1], shape.data['thk'], shape.data['t'], shape.data['center'], shape.data['radius'])
        elif shape.data['type'] == 'crossvault':
            zub, zlb = crossvault_ub_lb_update(XY[:, 0], XY[:, 1], shape.data['thk'], shape.data['t'], shape.data['xy_span'])
        elif shape.data['type'] == 'pointed_crossvault':
            zub, zlb = pointed_vault_ub_lb_update(XY[:, 0], XY[:, 1], shape.data['thk'], shape.data['t'], shape.data['xy_span'], hc=shape.data['hc'], he=shape.data['he'], hm=shape.data['hm'])
        elif shape.data['type'] == 'general':
            zub = shape.get_ub_pattern(XY)
            zlb = shape.get_lb_pattern(XY)
        else:
            raise Exception

        keysnan = []
        i = 0

        for key in self.vertices():
            ub_ = float(zub[i])
            lb_ = float(zlb[i])
            # print('(x,y): ({0:.3f}, {1:.3f})'.format(x, y), ' ub:', ub_)
            # print('(x,y): ({0:.3f}, {1:.3f})'.format(x, y), ' lb:', lb_)
            # if math.isnan(lb_):  # This check seems deprecated because the get_from_pattern seems to never return NaN...
            #     lb_ = -1 * shape.data['t']
            #     keysnan.append(key)
            #     x, y, _ = self.vertex_coordinates(key)
            #     print('Shape interpolation got NaN, check results (x,y): ({0:.3f}, {1:.3f})'.format(x, y))
            self.vertex_attribute(key, 'ub', value=ub_)
            self.vertex_attribute(key, 'lb', value=lb_)
            i += 1

        if keysnan:
            print(keysnan)
            plt = MeshPlotter(self)
            plt.draw_edges()
            plt.draw_vertices(keys=keysnan, radius=0.1, text='error', facecolor='FF0000')
            plt.show()

        return

    def selfweight_from_shape(self, shape):
        """Apply selfweight to the nodes of the form diagram based on the shape"""

        form_ = self.copy()
        total_selfweight = shape.compute_selfweight()

        XY = array(self.vertices_attributes('xy'))
        if shape.data['type'] == 'dome':
            zt = dome_zt_update(XY[:, 0], XY[:, 1], shape.data['radius'], shape.data['t'], shape.data['center'])
        elif shape.data['type'] == 'crossvault':
            zt = crossvault_middle_update(XY[:, 0], XY[:, 1],  shape.data['t'],  xy_span=shape.data['xy_span'])
        elif shape.data['type'] == 'pointed_crossvault':
            zt = pointed_vault_middle_update(XY[:, 0], XY[:, 1],  shape.data['t'],  xy_span=shape.data['xy_span'], hc=shape.data['hc'], he=shape.data['he'], hm=shape.data['hm'])
        else:
            zt = shape.get_middle_pattern(XY)

        i = 0
        for key in form_.vertices():
            z = zt[i]
            form_.vertex_attribute(key, 'z', value=z)
            self.vertex_attribute(key, 'target', value=z)
            i += 1

        pzt = 0
        for key in self.vertices():
            pz = form_.vertex_area(key)
            self.vertex_attribute(key, 'pz', value=pz)
            pzt += pz

        if shape.data['type'] == 'arch':
            pzt = 0
            for key in self.vertices():
                self.vertex_attribute(key, 'pz', value=1.0)
                if self.vertex_attribute(key, 'is_fixed') == True:
                    self.vertex_attribute(key, 'pz', value=0.5)
                pzt += self.vertex_attribute(key, 'pz')

        factor = total_selfweight/pzt

        for key in self.vertices():
            pzi = factor * self.vertex_attribute(key, 'pz')
            self.vertex_attribute(key, 'pz', value=pzi)

        return


    def selfweight_from_pattern(self, pattern, plot=False, tol=10e-4):
        """Apply selfweight to the nodes considering a different Form Diagram to locate loads. Warning, the base pattern has to coincide with nodes from the original form diagram"""

        form_ = pattern

        self.vertices_attribute('pz', 0.0)
        key_real_to_key = {}

        for key in form_.vertices():
            x, y, _ = form_.vertex_coordinates(key)
            for key_real in self.vertices():
                x_real, y_real, _ = self.vertex_coordinates(key_real)
                if x - tol < x_real < x + tol and y - tol < y_real < y + tol:
                    key_real_to_key[key_real] = key
                    break

        pzt = 0
        for key in key_real_to_key:
            pz = form_.vertex_attribute(key_real_to_key[key], 'pz')
            self.vertex_attribute(key, 'pz', value=pz)
            pzt += pz
        print('total load applied:', pzt)

        if plot:
            plotter = MeshPlotter(self, figsize=(10, 10))
            plotter.draw_edges()
            plotter.draw_vertices(text=key_real_to_key)
            plotter.show()

            plotter = MeshPlotter(form_, figsize=(10, 10))
            plotter.draw_edges()
            plotter.draw_vertices(text={key: key for key in form_.vertices()})
            plotter.show()

            plotter = MeshPlotter(self, figsize=(10, 10))
            plotter.draw_edges()
            plotter.draw_vertices(text={key: round(self.vertex_attribute(key, 'pz'), 1) for key in self.vertices()})
            plotter.show()

        return

    def apply_symmetry(self, center=[5.0, 5.0, 0.0], horizontal_only=False, vertical_only=False):

        self.edges_attribute('sym_key', None)
        self.vertices_attribute('sym_key', None)
        dist_checked = []
        dist_dict = {}

        # Symmetry on the independent edges

        if horizontal_only == False and vertical_only == False:
            for u, v in self.edges_where({'is_ind': True}):
                midpoint = self.edge_midpoint(u, v)
                dist = round(distance_point_point_xy(center, midpoint), 10)
                dist_dict[(u, v)] = dist
                if dist not in dist_checked:
                    dist_checked.append(dist)
        else:
            corners = mesh_bounding_box_xy(self)
            xs = [point[0] for point in corners]
            ys = [point[1] for point in corners]
            xc = (min(xs) + max(xs))/2
            yc = (min(ys) + max(ys))/2
            if horizontal_only is True:
                line = [[0, yc, 0.0], [10.0, yc, 0.0]]
                j = 0
            else:
                line = [[xc, 0, 0.0], [xc, 10.0, 0.0]]
                j = 1
            for u, v in self.edges_where({'is_ind': True}):
                midpoint = self.edge_midpoint(u, v)
                dist_line = round(distance_point_line_xy(midpoint, line), 10)
                dist_point = round(midpoint[j], 10)
                dist = [dist_point, dist_line]
                dist_dict[(u, v)] = dist
                if dist not in dist_checked:
                    dist_checked.append(dist)
        i = 0
        for dist in dist_checked:
            for u, v in dist_dict:
                if dist_dict[(u, v)] == dist:
                    self.edge_attribute((u, v), 'sym_key', i)
            i += 1

        # Symmetry on the Support's position

        dist_checked = []
        dist_dict = {}

        if horizontal_only == False and vertical_only == False:
            for key in self.vertices_where({'is_fixed': True}):
                point = self.vertex_coordinates(key)
                dist = round(distance_point_point_xy(center, point), 10)
                dist_dict[key] = dist
                if dist not in dist_checked:
                    dist_checked.append(dist)
        else:
            for key in self.vertices_where({'is_fixed': True}):
                point = self.vertex_coordinates(key)
                dist_line = round(distance_point_line_xy(point, line), 10)
                dist_point = round(point[j], 10)
                dist = [dist_point, dist_line]
                dist_dict[key] = dist
                if dist not in dist_checked:
                    dist_checked.append(dist)
        i = 0
        for dist in dist_checked:
            for key in dist_dict:
                if dist_dict[key] == dist:
                    self.vertex_attribute(key, 'sym_key', i)
            i += 1

        return

    def number_of_independents(self, printout=False):

        total = 0
        for key in self.edges_where({'is_ind': True}):
            total += 1
        if total == 0:
            print('Warning, no independent edges found!')

        if printout:
            print('Form has {0} independents'.format(total))

        return total

    def number_of_supports(self, printout=False):

        total = 0
        for key in self.vertices_where({'is_fixed': True}):
            total += 1
        if total == 0:
            print('Warning, no fixed points were found!')

        if printout:
            print('Form has {0} supports'.format(total))

        return total

    def number_of_sym_independents(self, printout=False):

        i_sym_max = 0
        for u, v in self.edges_where({'is_ind': True}):
            try:
                i_sym = self.edge_attribute((u, v), 'sym_key')
            except:
                print('Warning, no symmetry relation found!')
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        if printout:
            print('Form has {0} unique independents'.format(i_sym_max))

        return i_sym_max

    def number_of_sym_supports(self, printout=False):

        i_sym_max = 0
        for key in self.vertices_where({'is_fixed': True}):
            try:
                i_sym = self.vertex_attribute(key, 'sym_key')
            except:
                print('Warning, no symmetry relation found!')
            if i_sym > i_sym_max:
                i_sym_max = i_sym
        i_sym_max += 1

        if printout:
            print('Form has {0} unique supports'.format(i_sym_max))

        return i_sym_max

    def build_symmetry_matrix(self, printout=False):

        n = self.number_of_independents(printout=printout)
        k_unique = self.number_of_sym_independents(printout=printout)
        Asym = zeros((n - k_unique, n))
        uv_i_ind = {}

        i = 0
        for u, v in self.edges_where({'is_ind': True}):
            uv_i_ind[(u, v)] = i
            i += 1

        line = 0
        for id_sym in range(k_unique):
            i = 0
            for u, v in self.edges_where({'is_ind': True}):
                if self.edge_attribute((u, v), 'sym_key') == id_sym:
                    index = uv_i_ind[(u, v)]
                    if i == 0:
                        index0 = index
                    else:
                        Asym[line, index0] = 1
                        Asym[line, index] = -1
                        line += 1
                    i += 1
        if printout:
            plt.matshow(Asym)
            plt.show()

        return Asym

    def build_symmetry_matrix_supports(self, printout=False):

        n = self.number_of_supports()
        k_unique = self.number_of_sym_supports(printout=printout)
        Asym = zeros((n - k_unique, n))
        key_i_sup = {}

        i = 0
        for key in self.vertices_where({'is_fixed': True}):
            key_i_sup[key] = i
            i += 1

        line = 0
        for id_sym in range(k_unique):
            i = 0
            for key in self.vertices_where({'is_fixed': True}):
                if self.vertex_attribute(key, 'sym_key') == id_sym:
                    index = key_i_sup[key]
                    if i == 0:
                        index0 = index
                    else:
                        Asym[line, index0] = 1
                        Asym[line, index] = -1
                        line += 1
                    i += 1
        if printout:
            plt.matshow(Asym)
            plt.show()

        return Asym

    def assemble_symmetry_matrix(self, independents=True, supports=True, printout=False):
        """ Assemble Symmetry Matrix.

        Parameters
        ----------
        independents : bool (True)
            If independents should be considered in the symmetry matrix.
        supports : bool (True)
            If supports should be considered in the symmetry matrix.
        printout : bool (False)
            If Matrix will be visualised
        """

        if independents:
            Aind = self.build_symmetry_matrix(printout=printout)
        if supports:
            Ab = self.build_symmetry_matrix_supports(printout=printout)

        Aind_i, Aind_j = Aind.shape
        Ab_i, Ab_j = Ab.shape

        Aind0 = zeros((Aind_i, Ab_j))
        Ab0 = zeros((Ab_i, Aind_j))

        A1 = hstack([Aind, Aind0])
        A2 = hstack([Ab0, Ab])
        A = vstack([A1, A2])
        # A = vstack([hstack([Aind, Aind0]), hstack([Ab0, Ab])])

        if printout:
            plt.matshow(A)
            plt.show()

        return A

    # --------------------------------------------------------------- #
    # -----------------------TNA-CONECTION--------------------------- #
    # --------------------------------------------------------------- #

    def add_feet_(self, delete_face=False):
        """ Add feet to the support as in compas_tna.
        """

        if delete_face:
            self.delete_face(0)
        corners = list(self.vertices_where({'is_fixed': True}))
        self.vertices_attributes(('is_anchor', 'is_fixed'), (True, True), keys=corners)
        self.update_boundaries(self, feet=2)

        return self

    def remove_feet(self, openings=None, rmax=0.01):
        """ Remove feet from the support as in compas_tna.
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

    def evaluate_scale(self, function, bounds, n=100, plot=True):
        """ Evaluate a given objective function by scaling the form-diagram in the bounds specified.

        Parameters
        ----------
        form : obj
            The FormDiagram.
        function : method
            The objective function.
        bounds : list
            the lower and upper bound of the force densities.
        n : int
            Thenumbers of divisions inside the interval (bounds).
        plot : bool
            Plot form and force if desired.

        Returns
        -------
        form : obj
            The scaled form diagram.

        """

        r0 = bounds[0]
        stp = (bounds[1]-bounds[0])/n
        x = []
        y = []
        q0 = array([self.edge_attribute((u, v), 'q') for u, v in self.edges_where({'_is_edge': True})])[:, newaxis]

        form_ = deepcopy(self)

        k_i = form_.key_index()
        uv_i = form_.uv_index()

        for k in range(n):
            r = r0 + stp * k
            q = q0 * r
            x.append(r)

            for u, v in form_.edges_where({'_is_edge': True}):
                i = uv_i[(u, v)]
                [qi] = q[i]
                form_.edge_attribute((u, v), 'q', value=qi)

            form_ = z_from_form(form_)
            y.append(function(form_))

        import matplotlib
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(x, y)

        ax.set(xlabel='Scale (r)', ylabel='Energy (f)',
               title='Evaluate Energy by Scaling')
        ax.grid()

        pos = argmin(y)
        xmin = x[pos]
        ymin = y[pos]
        text = "x={:.3f}, y={:.3f}".format(xmin, ymin)
        if not ax:
            ax = plt.gca()
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        arrowprops = dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=60")
        kw = dict(xycoords='data', textcoords="axes fraction",
                  arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
        ax.annotate(text, xy=(xmin, ymin), xytext=(0.94, 0.96), **kw)

        if plot:
            plt.show()

        return xmin

    def scale_fdm(self, r):
        """ scale the FormDiagram of a factor r using FDM (all coordinates can change).

        Parameters
        ----------
        form : obj
            The FormDiagram.
        r : float
            The scaling factor on force densities.

        Returns
        -------
        form : obj
            The scaled form diagram.

        """

        uv_i = self.uv_index()
        q = array([self.edge_attribute((u, v), 'q') for u, v in self.edges_where({'_is_edge': True, '_is_external': False})])[:, newaxis]
        q = q * r

        for u, v in self.edges_where({'_is_external': False}):
            if self.edge_attribute((u, v), '_is_edge') is True:
                i = uv_i[(u, v)]
                [qi] = q[i]
                self.edge_attribute((u, v), 'q', value=qi)

        self = z_from_form(self)

        return self

    def scale_form(self, r):
        """ Scale the FormDiagram of a factor r using built-in FDM (only z-coordinates can change).

        Parameters
        ----------
        form : obj
            The FormDiagram.
        r : float
            The scaling factor on force densities.

        Returns
        -------
        form : obj
            The scaled form diagram.

        """

        from numpy import float64
        from scipy.sparse import diags
        from scipy.sparse.linalg import spsolve

        k_i = self.key_index()
        uv_i = self.uv_index()
        vcount = len(self.vertex)
        anchors = list(self.anchors())
        fixed = list(self.fixed())
        fixed = set(anchors + fixed)
        fixed = [k_i[key] for key in fixed]
        free = list(set(range(vcount)) - set(fixed))
        edges = [(k_i[u], k_i[v]) for u, v in self.edges_where({'_is_edge': True})]
        xyz = array(self.vertices_attributes('xyz'), dtype=float64)
        p = array(self.vertices_attributes(('px', 'py', 'pz')), dtype=float64)
        q = [self.edge_attribute((u, v), 'q') for u, v in self.edges_where({'_is_edge': True})]
        q = array(q, dtype=float64).reshape((-1, 1))
        C = connectivity_matrix(edges, 'csr')
        Ci = C[:, free]
        Cf = C[:, fixed]
        Cit = Ci.transpose()

        # TODO: Change to use function zq_from_qid

        q = q * r
        Q = diags([q.ravel()], [0])

        A = Cit.dot(Q).dot(Ci)
        B = Cit.dot(Q).dot(Cf)

        xyz[free, 2] = spsolve(A, p[free, 2] - B.dot(xyz[fixed, 2]))

        i = 0
        for key in self.vertices():
            self.vertex_attribute(key, 'z', xyz[i, 2])
            i = i + 1

        i = 0
        for u, v in self.edges_where({'_is_edge': True}):
            self.edge_attribute((u, v), 'q', q[i, 0])
            i = i + 1

        return

    def initialise_tna(self, zmax=5.0, method='nodal', plot=False, alpha=100.0, kmax=500, remove_feet=True, display=False):

        corners = list(self.vertices_where({'is_fixed': True}))
        self.vertices_attribute('is_anchor', True, keys=corners)
        self.edges_attribute('fmin', 0.0)
        self.edges_attribute('fmax', 10.0)
        leaves = False
        for u, v in self.edges_on_boundary():
            if self.edge_attribute((u, v), '_is_edge') is False:
                leaves = True
                break
        if leaves is False:
            self.update_boundaries()

        force = ForceDiagram.from_formdiagram(self)
        if plot:
            print('Plot of Primal')
            plotter = MeshPlotter(self, figsize=(10, 10))
            plotter.draw_edges(keys=[key for key in self.edges_where({'_is_edge': True})])
            plotter.draw_vertices(radius=0.05)
            plotter.draw_vertices(keys=[key for key in self.vertices_where({'is_anchor': True})], radius=0.10, facecolor='000000')
            plotter.show()
            print('Plot of Dual')
            force.plot()

        if method == 'nodal':
            horizontal_nodal(self, force, alpha=alpha, kmax=kmax)
        else:
            horizontal(self, force, alpha=alpha, kmax=kmax)

        vertical_from_zmax(self, zmax)
        # if leaves is False and remove_feet is True:
        #     self = self.remove_feet()

        if plot:
            print('Plot of Reciprocal')
            force.plot()

        return self

    def reciprocal_from_form(self, zmax=5.0, method='nodal', plot=False, alpha=100.0, kmax=500, remove_feet=True, display=False):

        corners = list(self.vertices_where({'is_fixed': True}))
        self.vertices_attribute('is_anchor', True, keys=corners)

        for key in self.vertices_where({'is_fixed': True}):
            react = abs(self.vertex_attribute(key, '_rx'))
            break
        leaves = False
        for u, v in self.edges_on_boundary():
            if self.edge_attribute((u, v), '_is_edge') == False:
                leaves = True
                break
        if leaves is False:
            self.update_boundaries(feet=2)

        for u, v in self.edges():
            if self.edge_attribute((u, v), '_is_external') is False:
                qi = self.edge_attribute((u, v), 'q')
                a = self.vertex_coordinates(u)
                b = self.vertex_coordinates(v)
                lh = distance_point_point_xy(a, b)
                self.edge_attribute((u, v), 'fmin', value=qi*lh)
                self.edge_attribute((u, v), 'fmax', value=qi*lh)
                self.edge_attribute((u, v), 'lmin', value=lh)
                self.edge_attribute((u, v), 'lmax', value=lh)
            # else:
            #     self.edge_attribute((u, v), 'fmin', value=react)
            #     self.edge_attribute((u, v), 'fmax', value=react)
            #     self.edge_attribute((u, v), 'lmin', value=1.0)
            #     self.edge_attribute((u, v), 'lmax', value=1.0)
            #     self.edge_attribute((u, v), 'q', value=react)

        if plot:
            plot_form(self).show()

        force = ForceDiagram.from_formdiagram(self)
        if plot:
            print('Plot of Primal')
            plotter = MeshPlotter(self, figsize=(10, 10))
            plotter.draw_edges(keys=[key for key in self.edges_where({'_is_edge': True})])
            plotter.draw_vertices(radius=0.05)
            plotter.draw_vertices(keys=[key for key in self.vertices_where({'is_fixed': True})], radius=0.10, facecolor='000000')
            plotter.show()
            print('Plot of Dual')
            force.plot()

        if method == 'nodal':
            horizontal_nodal(self, force, alpha=alpha, kmax=kmax, display=False)
        else:
            horizontal(self, force, alpha=alpha, kmax=kmax, display=False)

        # Vertical Equilibrium with no updated loads

        if leaves is False and remove_feet is True:
            self = self.remove_feet()

        if plot:
            print('Plot of Reciprocal')
            force.plot()

        return self, force

    def vertex_projected_area(self, key):
        """Compute the projected tributary area of a vertex.

        Parameters
        ----------
        key : int
            The identifier of the vertex.

        Returns
        -------
        float
            The projected tributary area.

        Example
        -------
        >>>

        """

        from compas.geometry import subtract_vectors
        from compas.geometry import length_vector
        from compas.geometry import cross_vectors

        area = 0.

        p0 = self.vertex_coordinates(key)
        p0[2] = 0

        for nbr in self.halfedge[key]:
            p1 = self.vertex_coordinates(nbr)
            p1[2] = 0
            v1 = subtract_vectors(p1, p0)

            fkey = self.halfedge[key][nbr]
            if fkey is not None:
                p2 = self.face_centroid(fkey)
                p2[2] = 0
                v2 = subtract_vectors(p2, p0)
                area += length_vector(cross_vectors(v1, v2))

            fkey = self.halfedge[nbr][key]
            if fkey is not None:
                p3 = self.face_centroid(fkey)
                p3[2] = 0
                v3 = subtract_vectors(p3, p0)
                area += length_vector(cross_vectors(v1, v3))

        return 0.25 * area

    # def adapt_objective(form, zrange=[3.0, 8.0], objective='loadpath', method='nodal', discr=100, plot=False,
    #                     delete_face=False, alpha=100.0, kmax=100, display=False, amax=2.0, rmax=0.01):

    #     form = add_feet(form, delete_face=delete_face, plot=plot)
    #     force = ForceDiagram.from_formdiagram(form)

    #     plot_force(force, form).show()

    #     if method == 'nodal':
    #         horizontal_nodal(form, force, alpha=alpha, kmax=kmax, display=False)
    #     else:
    #         horizontal(form, force, alpha=alpha, kmax=kmax, display=False)

    #     plot_force(force, form).show()

    #     a = evaluate_a(form, plot=plot)
    #     if a > amax:
    #         print('High Angle deviations!')

    #     # Vertical Equilibrium with no updated loads

    #     if objective == 'loadpath':
    #         f = loadpath
    #     if objective == 'target':
    #         f = energy

    #     scale = []
    #     form0 = scale_form(form, 1.0)
    #     z = [form0.vertex_attribute(key, 'z') for key in form0.vertices()]
    #     z_ = max(z)
    #     scale.append(z_ / zrange[1])
    #     scale.append(z_ / zrange[0])
    #     print(scale)

    #     best_scl = evaluate_scale(form0, f, scale, n=discr, plot=plot)

    #     print('Best Scale is: {0}'.format(best_scl))
    #     form = scale_form(form0, best_scl)

    #     r = residual(form, plot=plot)
    #     if r > rmax:
    #         print('High residual forces!')

    #     if plot:
    #         plot_form(form).show()
    #         plot_force(force, form).show()

    #     return form

    # def sym_on_openings(form, xy_span=[[0.0, 10.0], [0.0, 10.0]], exc_length=1.0):

    #     gkey_key = form.gkey_key()

    #     y1 = xy_span[1][1]
    #     y0 = xy_span[1][0]
    #     x1 = xy_span[0][1]
    #     x0 = xy_span[0][0]

    #     bndr = form.vertices_on_boundary()
    #     lines = [form.edge_coordinates(u, v) for u, v in form.edges()]

    #     for key in bndr:
    #         if form.vertex_attribute(key, 'is_fixed') == False:
    #             x, y, _ = form.vertex_coordinates(key)
    #             if x == x1:
    #                 lines.append([[x, y, 0.0], [x + exc_length, y, 0.0]])
    #             if x == x0:
    #                 lines.append([[x, y, 0.0], [x - exc_length, y, 0.0]])
    #             if y == y1:
    #                 lines.append([[x, y, 0.0], [x, y + exc_length, 0.0]])
    #             if y == y0:
    #                 lines.append([[x, y, 0.0], [x, y - exc_length, 0.0]])

    #     form_ = FormDiagram.from_lines(lines, delete_boundary_face=False)
    #     # form.delete_face(0)

    #     form_.update_default_vertex_attributes({'is_roller': False})
    #     form_.update_default_vertex_attributes({'is_fixed': False})
    #     form_.update_default_edge_attributes({'q': 1, 'is_symmetry': False})
    #     form_.attributes['loadpath'] = 0
    #     form_.attributes['indset'] = []

    #     for key_ in form_.vertices():
    #         coord_ = form_.vertex_coordinates(key_)
    #         try:
    #             key = gkey_key[geometric_key(coord_)]
    #             pz = form.vertex_attribute(key, 'pz')
    #             px = form.vertex_attribute(key, 'px')
    #             py = form.vertex_attribute(key, 'py')
    #             fixed = form.vertex_attribute(key, 'is_fixed')
    #             lb = form.vertex_attribute(key, 'lb')
    #             ub = form.vertex_attribute(key, 'ub')
    #             target = form.vertex_attribute(key, 'target')
    #             form_.vertex_attribute(key_, 'pz', value=pz)
    #             form_.vertex_attribute(key_, 'px', value=px)
    #             form_.vertex_attribute(key_, 'py', value=py)
    #             form_.vertex_attribute(key_, 'is_fixed', value=fixed)
    #             form_.vertex_attribute(key_, 'lb', value=lb)
    #             form_.vertex_attribute(key_, 'ub', value=ub)
    #             form_.vertex_attribute(key_, 'target', value=target)
    #         except:
    #             form_.vertex_attribute(key_, 'is_fixed', value=True)
    #             form_.vertex_attribute(key_, 'pz', value=0.0)
    #             form_.vertex_attribute(key_, 'px', value=0.0)
    #             form_.vertex_attribute(key_, 'py', value=0.0)
    #             form_.vertex_attribute(key_, 'lb', value=0.0)
    #             form_.vertex_attribute(key_, 'ub', value=0.0)
    #             form_.vertex_attribute(key_, 'target', value=0.0)
    #             ngb = form_.vertex_neighbors(key_)[0]
    #             form_.edge_attribute((key_, ngb), 'is_symmetry', True)

    #     return form_
