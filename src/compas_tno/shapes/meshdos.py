from compas.geometry import subtract_vectors
from compas.geometry import length_vector
from compas.geometry import cross_vectors
from compas.geometry import normalize_vector
from compas.geometry import is_point_in_convex_polygon_xy
from compas.geometry import centroid_points
from compas.geometry import normal_polygon
from compas.geometry import norm_vector
from compas.geometry import sum_vectors
from compas.geometry import scale_vector
from compas.geometry import angle_vectors

from compas.datastructures import Mesh

import math


class MeshDos(Mesh):
    """The ``MeshDos`` inherits compas ``Mesh`` object and add some features to create Intrados and Extrados.

    Notes
    -----
    A ``MeshDos`` has the following constructor functions

    *   ``from_library`` : Construct the shape from a dictionary with instructions, the library supports the creation of parametric arches, domes, and vaults.
    *   ``from_rhinomesh`` : Construct Extrados, Intrados and Middle surfaces from RhinoMeshes.
    *   ``from_rhinosurface`` : Construct Extrados, Intrados and Middle surfaces using the U and V isolines.
    *   ``from_pointcloud`` : Construct Extrados, Intrados and Middle surfaces interoplating a pointcloud.

    """

    def __init__(self):  # add 'is_outside': False as default
        super(Mesh, self).__init__()

    @classmethod
    def from_mesh(cls, mesh):
        """ Construct object from an existing COMPAS mesh.

        Parameters
        ----------
        mesh : Mesh
            Base mesh.

        Returns
        -------
        mesh
            A mesh MeshDos mesh.
        """

        vertices, faces = mesh.to_vertices_and_faces()

        return cls().from_vertices_and_faces(vertices, faces)

    @classmethod
    def from_formdiagram_attribute(cls, formdiagram, attribute='lb'):
        """Create a MeshDos object from a base form diagram and an attribute.

        Parameters
        ----------
        formdiagram : :class:`~compas_tno.diagrams.FormDiagram`
            The form diagram to copy
        attribute : str, optional
            The attribute to carry to the ``z``, by default 'lb'

        Returns
        -------
        mesh
            MeshDos
        """

        vertices, faces = formdiagram.to_vertices_and_faces()
        mesh = cls().from_vertices_and_faces(vertices, faces)

        for key in formdiagram.vertices():
            z_ = formdiagram.vertex_attribute(key, attribute)
            mesh.vertex_attribute(key, 'z', z_)

        return mesh

    def offset_mesh(self, n=0.1, direction='up', t=0.0):
        """Offset the mesh upwards considering it as intrados.

        Parameters
        ----------
        n : float, optional
            Normal distance of offset, by default 0.1
        direction : str, optional
            Direction of the offset, by default 'up'
        t : float, optional
            Distance considered in the supports, by default 0.0

        Returns
        -------
        mesh
            MeshDos
        """

        offset_list = []
        mesh_copy = self.copy()

        for key in self.vertices():
            normal = self.vertex_attribute(key, 'n')
            # normal = new_mesh.vertex_normal(key)
            z = self.vertex_coordinates(key)[2]
            if self.vertex_attribute(key, 'is_outside'):
                offset_list.append(t)
            else:
                deviation = 1/math.sqrt(1/(1 + (normal[0]**2 + normal[1]**2)/normal[2]**2))
                # print('deviation:', deviation)
                if direction == 'up':
                    offset_list.append(z + n*deviation*norm_vector(normal))  # Experimenting with this normal norm!
                else:
                    offset_list.append(z - n*deviation*norm_vector(normal))  # Experimenting with this normal norm!

        for i, key in enumerate(mesh_copy.vertices()):
            mesh_copy.vertex_attribute(key, 'z', offset_list[i])

        return mesh_copy

    def offset_up_and_down(self, n=1.0, t=0.0):
        """Offset the mesh up and down considering it as middle surface.

        Parameters
        ----------
        n : float, optional
            Percentage of the normal vector to use, by default 1.0.
        t : float, optional
            Distance considered in the supports, by default 0.0

        Returns
        -------
        mesh_ub, mesh_lb
            Extrados mesh and intrados mesh.

        """

        ub_update = []
        lb_update = []
        mesh_ub = self.copy()
        mesh_lb = self.copy()

        for key in self.vertices():
            nub = self.vertex_attribute(key, 'nub')
            nlb = self.vertex_attribute(key, 'nlb')
            z = self.vertex_coordinates(key)[2]
            dev_ub = 1/math.sqrt(1/(1 + (nub[0]**2 + nub[1]**2)/nub[2]**2))
            if self.vertex_attribute(key, 'is_outside'):
                dev_lb = 0.0
            else:
                dev_lb = 1/math.sqrt(1/(1 + (nlb[0]**2 + nlb[1]**2)/nlb[2]**2))
            ub_update.append(z + n * dev_ub * norm_vector(nub))
            lb_update.append(z - n * dev_lb * norm_vector(nlb))

        for i, key in enumerate(mesh_ub.vertices()):
            mesh_ub.vertex_attribute(key, 'z', ub_update[i])
            mesh_lb.vertex_attribute(key, 'z', lb_update[i])

        return mesh_ub, mesh_lb

    def store_normals(self, correct_creases=False, printout=False, plot=False):
        """Store normals in a MeshDos object

        Parameters
        ----------
        correct_creases : bool, optional
            If creases should be corrected, by default False
        printout : bool, optional
            If prints should appear in the screen, by default False
        plot : bool, optional
            If plots with normals, by default False

        Returns
        -------
        None
            Normals added as attributes
        """

        for key in self.vertices():
            self.vertex_attribute(key, 'n', self.vertex_normal(key))

        if correct_creases:
            self.magnify_normal_at_ribs(printout=printout, plot=plot)

        return

    def identify_creases_by_angle(self, deviation=20.0):
        """ Identify creses in the structure based on a threshold angle limit.

        Parameters
        ----------
        deviation : float, optional
            Angle (deg) deviation of normals along an edge so it is considered a crease. Default is 20.0.

        Returns
        -------
        None
            Identification made in the shape
        """

        self.vertices_attribute('is_crease', False)
        self.edges_attribute('is_crease', False)
        for u, v in self.edges():
            faces = self.edge_faces(u, v)
            if faces[0] is not None and faces[1] is not None:
                normals = [self.face_normal(face) for face in faces]
                dev = angle_vectors(normals[0], normals[1], deg=True)
                if dev > deviation:
                    self.vertex_attribute(u, 'is_crease', True)
                    self.vertex_attribute(v, 'is_crease', True)
                    self.edge_attribute((u, v), 'is_crease', True)

        return

    def identify_creases_at_diagonals(self, xy_span=[[0.0, 10.0], [0.0, 10.0]]):
        """Identify creses in the structure based on the diagonal of a rectangular.

        Parameters
        ----------
        xy_span : list, optional
            Span of the pattern, by default [[0.0, 10.0], [0.0, 10.0]]

        Returns
        -------
        None
            Creases are marked with attributes
        """

        self.vertices_attribute('is_crease', False)
        self.edges_attribute('is_crease', False)
        tol = 10E-3

        y1 = xy_span[1][1]
        y0 = xy_span[1][0]
        x1 = xy_span[0][1]
        x0 = xy_span[0][0]

        for key in self.vertices():
            xi, yi, _ = self.vertex_coordinates(key)
            if abs(yi - (y0 + (y1 - y0)/(x1 - x0) * (xi - x0))) <= tol or abs(yi - (y1 - (y1 - y0)/(x1 - x0) * (xi - x0))) <= tol:
                self.vertex_attribute(key, 'is_crease', True)

        for u, v in self.edges():
            if self.vertex_attribute(u, 'is_crease') and self.vertex_attribute(v, 'is_crease'):
                self.edge_attribute((u, v), 'is_crease', True)

        return

    def magnify_normal_at_ribs(self):
        """Magnify normals at ribs.
        """

        keys = list(self.vertices_where({'is_crease': True}))

        for vkey in keys:
            vfaces = self.vertex_faces(vkey, ordered=True)
            vfaces = vfaces + [vfaces[0]]
            right = []
            left = []
            for i in range(len(vfaces) - 2):
                halfedge = self.face_adjacency_halfedge(vfaces[i], vfaces[i + 1])
                if i == 0:
                    right.append(vfaces[i])
                    if self.edge_attribute(halfedge, 'is_crease'):
                        left.append(vfaces[i + 1])
                    else:
                        right.append(vfaces[i + 1])
                else:
                    if self.edge_attribute(halfedge, 'is_crease'):
                        if vfaces[i] in right:
                            left.append(vfaces[i + 1])
                        else:
                            right.append(vfaces[i + 1])
                    else:
                        if vfaces[i] in right:
                            right.append(vfaces[i + 1])
                        else:
                            left.append(vfaces[i + 1])

            normal_right = normalize_vector(centroid_points([self.face_normal(fkey, False) for fkey in right]))
            normal_left = normalize_vector(centroid_points([self.face_normal(fkey, False) for fkey in left]))
            angle = angle_vectors(normal_right, normal_left)/2
            n = scale_vector(normalize_vector(sum_vectors([normal_right, normal_left])), 1/math.cos(angle))
            self.vertex_attribute(vkey, 'n', n)

        return

    def scale_normals_with_ub_lb(self, zub, zlb, tol=10e-3):
        """Scale the normals on the middle surface to match zub and zlb.

        Parameters
        ----------
        zub : array (n x 1)
            Height of the nodes of extrados.
        zlb : array (n x 1)
            Height of the nodes of intrados.
        tol : float, optional
            The tolerance, by default 10e-3

        """

        self.vertices_attribute('is_outside', False)
        i = 0
        for key in self.vertices():
            n = self.vertex_attribute(key, 'n')
            zti = self.vertex_attribute(key, 'z')
            deviation = 1/math.sqrt(1/(1 + (n[0]**2 + n[1]**2)/n[2]**2))
            fac_ub = (zub[i] - zti)/deviation
            fac_lb = (zlb[i] - zti)/deviation
            nub = scale_vector(n, fac_ub/norm_vector(n))
            nlb = scale_vector(n, fac_lb/norm_vector(n))
            self.vertex_attribute(key, 'nub', nub)
            self.vertex_attribute(key, 'nlb', nlb)
            if abs(fac_lb) < tol:
                self.vertex_attribute(key, 'is_outside', True)
            i += 1

        return

    def get_xy_face_normals(self, XY):
        """ Get normal based on the face normal in which the points are in the plan.

        Parameters
        ----------
        XY : list
            List of Points

        Returns
        -------
        normals : list
            List with normals
        """

        face_coords = [self.face_coordinates(fkey) for fkey in self.faces()]
        normals = []
        for xy in XY:
            normals_for_pt = []
            for face_coord in face_coords:
                i = 0
                if is_point_in_convex_polygon_xy(xy, face_coord):
                    normal = normal_polygon(face_coord, unitized=False)
                    normals_for_pt.append(normal)
                    i = i+1
            if len(normals_for_pt) == 1:
                normal = normalize_vector(normals_for_pt[0])
            else:
                normal = normalize_vector(centroid_points(normals_for_pt))
            normals.append(normal)
        return normals

    def round_xy_boundary(self, span=[[0, 10], [0, 10]], tol=10e-4):
        """Prevent rounding errors by rounding boundary vertices

        Parameters
        ----------
        span : list, optional
            Span of the pattern, by default [[0, 10], [0, 10]]
        tol : float, optional
            Tolerance, by default 10e-4

        Returns
        -------
        None
            Pattern modified in place
        """

        for key in self.vertices():
            x, y, z = self.vertex_coordinates(key)
            if abs(x - span[0][0]) < tol:
                self.vertex_attribute(key, 'x', span[0][0])
            if abs(x - span[1][0]) < tol:
                self.vertex_attribute(key, 'x', span[1][0])
            if abs(y - span[0][1]) < tol:
                self.vertex_attribute(key, 'y', span[0][1])
            if abs(y - span[1][1]) < tol:
                self.vertex_attribute(key, 'y', span[1][1])

        return

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

        """

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
