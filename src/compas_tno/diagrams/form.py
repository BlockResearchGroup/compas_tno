
from compas.utilities import geometric_key
from compas.numerical import connectivity_matrix
from random import shuffle
from numpy import max
from numpy import min
from numpy import zeros
from numpy import newaxis
from numpy import array
import math

from compas_tno.diagrams.diagram_arch import create_arch
from compas_tno.diagrams.diagram_rectangular import create_cross_form
from compas_tno.diagrams.diagram_rectangular import create_fan_form
from compas_tno.diagrams.diagram_circular import create_circular_radial_form
from compas_tno.diagrams.diagram_circular import create_circular_radial_spaced_form
from compas_tno.diagrams.diagram_circular import create_circular_spiral_form
from compas_tno.diagrams.diagram_rectangular import create_ortho_form

from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_force

from compas_tno.diagrams import ForceDiagram
from compas_tna.diagrams import FormDiagram

from compas_tna.equilibrium import horizontal
from compas_tna.equilibrium import horizontal_nodal
from compas_tna.equilibrium import vertical_from_zmax



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

        if form_type == 'arch':
            form = cls().create_arch(H=data['H'], L=data['L'], x0=data['x0'], total_nodes=data['total_nodes'])
        if form_type == 'cross_fd':
            form = cls().create_cross_form(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data['fix'])
        if form_type == 'fan_fd':
            form = cls().create_fan_form(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data['fix'])
        if form_type == 'ortho':
            form = cls().create_ortho_form(xy_span=data['xy_span'], discretisation=data['discretisation'], fix=data['fix'])
        if form_type == 'radial_fd':
            form = cls().create_circular_radial_form(center=data['center'], radius=data['radius'], discretisation=data['discretisation'],
                                                     r_oculus=data['r_oculus'], diagonal=data['diagonal'], partial_diagonal=data['partial_diagonal'])
        if form_type == 'radial_spaced_fd':
            form = cls().create_circular_radial_spaced_form(
                center=data['center'], radius=data['radius'], discretisation=data['discretisation'], r_oculus=data['r_oculus'], diagonal=data['diagonal'], partial_diagonal=data['partial_diagonal'])
        if form_type == 'spiral_fd':
            form = cls().create_circular_spiral_form(center=data['center'], radius=data['radius'], discretisation=data['discretisation'], r_oculus=data['r_oculus'])

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


# --------------------------------------------------------------- #
# -----------------------ULTILITIES------------------------------ #
# --------------------------------------------------------------- #

    def overview_forces(self):
        """ Quick overview of the important parameters of a network.
        """

        f = []
        q = []
        z = []
        pz = 0

        lp = 0

        for u, v in self.edges_where({'is_external': False}):
            if self.edge_attribute((u, v), 'is_edge') is True and self.edge_attribute((u, v), 'is_symmetry') is False:
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

        form_ = FormDiagram.from_lines(edges, delete_boundary_face=False)
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
            if vertex.get('is_external') == False:
                z = vertex.get('z')
                s = vertex.get('target')
                w = vertex.get('weight', 1.0)
                f += w * (z - s)**2

        return f

    def loadpath(self):
        """ Compute loadpath in the current configuration.
        """

        lp = 0
        for u, v in self.edges_where({'is_external': False}):
            if self.edge_attribute((u, v), 'is_edge') is True and form.edge_attribute((u, v), 'is_symmetry') is False:
                qi = self.edge_attribute((u, v), 'q')
                li = self.edge_length(u, v)
                lp += qi*li**2

        self.attributes['loadpath'] = lp

        return lp

    def delete_boundary_edges(self):
        """ Delete boundary edges on a diagram.
        """

        for u, v in self.edges_on_boundary():
            self.edge_attribute((u, v), 'is_edge', False)

        return self

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
        q = array([attr['q'] for u, v, attr in self.edges(True)])[:, newaxis]

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
        for u, v, attr in self.edges_where({'is_edge': True}, True):
            a = attr['a']
            a_total += a
            l = self.edge_length(u, v)
            a = a*l
            if a > a_max:
                a_max = a
        if plot is True:
            print('Angle Deviation  Max: {0}'.format(a_max))
            print('Angle Deviation  Total: {0}'.format(a_total))

        return a_max


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

        for u, v in self.edges_where({'is_edge': True, 'is_external': False}):
            s = self.vertex_coordinates(u)
            e = self.vertex_coordinates(v)
            lines.append([s, e])
            qs[geometric_key(self.edge_midpoint(u, v))] = self.edge_attribute((u, v), 'q')

        fixed = [geometric_key(self.vertex_coordinates(key)) for key in self.vertices_where({'is_anchor': True})]
        zs = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_coordinates(key)[2] for key in self.vertices_where({'is_external': False})}
        pz = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'pz') for key in self.vertices()}
        target = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'target') for key in self.vertices()}
        lb = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'lb') for key in self.vertices()}
        ub = {geometric_key(self.vertex_coordinates(key)[:2] + [0]): self.vertex_attribute(key, 'ub') for key in self.vertices()}

        form_ = FormDiagram.from_lines(lines)
        form_.update_default_edge_attributes({'q': 1, 'is_symmetry': False, 'is_edge': True})
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

    def evaluate_scale(self, function, bounds, n = 100, plot = True):

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
        q0 = array([self.edge_attribute((u,v),'q') for u, v in self.edges_where({'is_edge': True})])[:, newaxis]

        form_ = deepcopy(self)

        k_i  = form_.key_index()
        uv_i = form_.uv_index()

        for k in range(n):
            r = r0 + stp * k
            q = q0 * r
            x.append(r)

            for u, v in form_.edges_where({'is_edge': True}):
                i = uv_i[(u,v)]
                [qi] = q[i]
                form_.edge_attribute((u,v),'q',value=qi)

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
        text= "x={:.3f}, y={:.3f}".format(xmin, ymin)
        if not ax:
            ax=plt.gca()
        bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
        arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
        kw = dict(xycoords='data',textcoords="axes fraction",
                arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
        ax.annotate(text, xy=(xmin, ymin), xytext=(0.94,0.96), **kw)

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
        q = array([self.edge_attribute((u,v),'q') for u, v in self.edges_where({'is_edge': True, 'is_external' : False})])[:, newaxis]
        q = q * r

        for u, v in self.edges_where({'is_external': False}):
            if self.edge_attribute((u,v),'is_edge') is True:
                i = uv_i[(u,v)]
                [qi] = q[i]
                self.edge_attribute((u,v),'q',value=qi)

        self = z_from_form(self)

        return self

    def scale_form(self,r):
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

        k_i     = self.key_index()
        uv_i    = self.uv_index()
        vcount  = len(self.vertex)
        anchors = list(self.anchors())
        fixed   = list(self.fixed())
        fixed   = set(anchors + fixed)
        fixed   = [k_i[key] for key in fixed]
        free    = list(set(range(vcount)) - set(fixed))
        edges   = [(k_i[u], k_i[v]) for u, v in self.edges_where({'is_edge': True})]
        xyz     = array(self.get_vertices_attributes('xyz'), dtype=float64)
        p       = array(self.get_vertices_attributes(('px', 'py', 'pz')), dtype=float64)
        q       = [attr.get('q', 1.0) for u, v, attr in self.edges_where({'is_edge': True}, True)]
        q       = array(q, dtype=float64).reshape((-1, 1))
        C       = connectivity_matrix(edges, 'csr')
        Ci      = C[:, free]
        Cf      = C[:, fixed]
        Cit     = Ci.transpose()

        q = q * r
        Q = diags([q.ravel()], [0])

        A       = Cit.dot(Q).dot(Ci)
        B       = Cit.dot(Q).dot(Cf)

        xyz[free, 2] = spsolve(A,p[free, 2] - B.dot(xyz[fixed, 2]))

        for key, attr in self.vertices(True):
            index = k_i[key]
            attr['z']  = xyz[index, 2]

        for u, v, attr in self.edges_where({'is_edge': True}, True):
            index = uv_i[(u, v)]
            attr['q'] = q[index, 0]

        return self


# ------------------------------------------------------------------- #
# -----------------------CHECK IMPORTANCE---------------------------- #
# ------------------------------------------------------------------- #

    def initialise_tna(self, zmax=5.0, method='nodal', plot=False, alpha=100.0, kmax=500, display=False):

        corners = list(self.vertices_where({'is_fixed': True}))
        self.vertices_attribute('is_anchor', True, keys=corners)
        self.edges_attribute('fmin', 0.0)
        leaves = False
        for u, v in self.edges_on_boundary():
            if self.edge_attribute((u,v), 'is_edge') == False:
                leaves = True
                break
        if leaves is False:
            self.update_boundaries(feet=2)
        force = ForceDiagram.from_formdiagram(self)
        if plot:
            print('Plot of Dual')
            force.plot()
            # plot_form(self, show_q=False, fix_width=True).show()

        if method == 'nodal':
            horizontal_nodal(self, force, alpha=alpha, kmax=kmax, display=False)
        else:
            horizontal(self, self, alpha=alpha, kmax=kmax, display=False)

        # Vertical Equilibrium with no updated loads

        vertical_from_zmax(self, zmax)
        if leaves is False:
            self = self.remove_feet()

        # form0 = copy(self)
        # form_ = scale_form(form0, 1.0)
        # z = [form_.vertex_attribute(key, 'z') for key in form_.vertices()]
        # z_ = max(z)
        # scale = (z_ / zmax)
        # form = scale_form(form0, scale)

        if plot:
            print('Plot of Reciprocal')
            force.plot()

        # if plot:
        #     plot_form(self).show()
        #     plot_force(force, self).show()
        #     force.plot()

        return self

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