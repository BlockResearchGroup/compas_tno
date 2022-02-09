from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import Point

from compas_tno.shapes import Shape

from compas_view2 import app
from compas_view2.shapes import Arrow
from compas_view2.shapes import Text


__all__ = ['Viewer']


class Viewer(object):  # CHHANGE view_xxx functions to draw_xxx functions
    """A Class for view 3D thrust networks and shapes.

    Parameters
    ----------
    thrust : FormDiagram, optional
        The FormDiagram you want to plot, by default None
    shape : Shape, optional
        The Shape of masonry to plot, by default None

    Attributes
    ----------

    None

    """

    def __init__(self, thrust=None, shape=None, force=None, **kwargs):

        super().__init__(**kwargs)
        self.title = 'Viewer'
        self.app = None
        self.thrust = thrust
        self.shape = shape
        self.force = force
        self.settings = {
            'show.thrust': True,
            'show.shape': True,
            'show.reactions': True,
            'show.reactionlabels': True,
            'show.cracks': True,
            'show.vertex.outside': True,
            'show.edge.thickness': True,

            'camera.target': [0, 0, 0],
            'camera.distance': 40,
            'camera.rx': -45,
            'camera.rz': 45,
            'camera.fov': 40,
            'camera.show.grid': True,
            'camera.show.axis': True,

            'size.vertex': 15.0,
            'size.edge.max_thickness': 15.0,
            'size.edge.base_thickness': 1.0,
            'size.edge.normals': 0.5,
            'size.reaction.head_width': 0.1,
            'size.reaction.body_width': 0.01,
            'size.reactionlabel': 20,

            'force.scale': 0.05,
            'force.anchor': [0, -10, 0],

            'scale.reactions': 0.01,
            'scale.loads': 0.001,
            'opacity.shapes': 0.5,

            'color.edges.thrust': (255, 0, 0),
            'color.edges.reactions': (125, 125, 125),
            'color.vertex.extrados': (0, 125, 0),
            'color.vertex.intrados': (0, 0, 255),
            'color.vertex.outside': (200, 200, 200),
            'color.mesh.intrados': (125, 125, 125),
            'color.mesh.extrados': (125, 125, 125),
            'color.mesh.middle': (125, 125, 125),
            'color.mesh.general': (0, 0, 0),

            'tol.forces': 1e-3,
        }
        self.initiate_app()

    def initiate_app(self):
        """ Initiate the App with the default camera options

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app = app.App(width=1600, height=900, show_grid=self.settings['camera.show.grid'])

        self.app.view.camera.target = self.settings['camera.target']
        self.app.view.camera.distance = self.settings['camera.distance']
        self.app.view.camera.rx = self.settings['camera.rx']
        self.app.view.camera.rz = self.settings['camera.rz']
        self.app.view.camera.fov = self.settings['camera.fov']

    def show_solution(self):
        """ Show the thrust network, with the shape according to the settings

        Returns
        -------
        None
            The objects are updated in place
        """

        self.view_thrust()
        self.view_cracks()
        self.view_shape()
        self.view_reactions()
        self.show()

    def show(self):
        """ Display the viewer in the screen

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app.show()

    def save(self, path='temp/fig.png'):
        """ Save to a path

        Returns
        -------
        None
            The objects are updated in place
        """

        # self.app.save(path)
        self.app.view.init()
        self.app.record = True
        self.app.view.paint()
        print(self.app.recorded_frames)

    def clear(self):
        """Clear the viewer elements

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app.view.objects = {}

    def view_thrust(self):
        """View thrust network according to the settings

        Returns
        -------
        None
            The viewer is updated in place
        """

        base_thick = self.settings['size.edge.base_thickness']
        max_thick = self.settings['size.edge.max_thickness']
        thickness = self.settings['show.edge.thickness']

        if thickness:
            forces = [self.thrust.edge_attribute((u, v), 'q') * self.thrust.edge_length(u, v) for u, v in self.thrust.edges_where({'_is_edge': True})]
            fmax = max(abs(max(forces)), abs(min(forces)))

        for u, v in self.thrust.edges_where({'_is_edge': True}):
            Xu = self.thrust.vertex_coordinates(u)
            Xv = self.thrust.vertex_coordinates(v)
            line = Line(Xu, Xv)
            if not thickness:
                self.app.add(line, name=str((u, v)), linewidth=base_thick, color=_norm(self.settings['color.edges.thrust']))
                continue
            q = self.thrust.edge_attribute((u, v), 'q')
            length = self.thrust.edge_length(u, v)
            force = abs(q*length)
            thk = force/fmax * max_thick
            if force > self.settings['tol.forces']:
                self.app.add(line, name=str((u, v)), linewidth=thk, color=_norm(self.settings['color.edges.thrust']))

    def view_cracks(self):
        """View cracks according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        intrad = 1
        extrad = 1
        out = 1

        for key in self.thrust.vertices():
            lb = self.thrust.vertex_attribute(key, 'lb')
            ub = self.thrust.vertex_attribute(key, 'ub')
            x, y, z = self.thrust.vertex_coordinates(key)
            if self.settings['show.cracks']:
                if abs(ub - z) < self.settings['tol.forces']:
                    self.app.add(Point(x, y, z), name="Extrados (%s)" % extrad, color=_norm(self.settings['color.vertex.extrados']), size=self.settings['size.vertex'])
                    extrad += 1
                elif abs(lb - z) < self.settings['tol.forces']:
                    self.app.add(Point(x, y, z), name="Intrados (%s)" % intrad, color=_norm(self.settings['color.vertex.intrados']), size=self.settings['size.vertex'])
                    intrad += 1
            if self.settings['show.vertex.outside']:
                if z > ub:
                    self.app.add(Point(x, y, z), name="Outside - Intra (%s)" % out, color=_norm(self.settings['color.vertex.outside']), size=self.settings['size.vertex'])
                    out += 1
                elif z < lb:
                    self.app.add(Point(x, y, z), name="Outside - Extra (%s)" % out, color=_norm(self.settings['color.vertex.outside']), size=self.settings['size.vertex'])
                    out += 1

    def view_shape(self):
        """View the shape (intrados + extrados) according to the settings

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        if self.shape:
            datashape = self.shape.datashape.copy()
            if datashape['type'] in ['dome']:
                datashape['type'] = 'dome_polar'
                shape = Shape.from_library(datashape)
            # elif datashape['type'] == 'arch':
            #     print('WIP = Special Plot for arch')
            else:
                shape = self.shape

        else:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        self.app.add(shape.intrados, name="Intrados", show_edges=False, opacity=self.settings['opacity.shapes'], color=_norm(self.settings['color.mesh.intrados']))
        self.app.add(shape.extrados, name="Extrados", show_edges=False, opacity=self.settings['opacity.shapes'], color=_norm(self.settings['color.mesh.extrados']))

    def view_middle_shape(self):
        """ View the middle of the shape according to the settings

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        shape = self.shape
        if not shape:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        vertices_middle, faces_middle = shape.middle.to_vertices_and_faces()
        mesh_middle = Mesh.from_vertices_and_faces(vertices_middle, faces_middle)
        self.app.add(mesh_middle, name="Middle", show_edges=False, opacity=self.settings['opacity.shapes'], color=_norm(self.settings['color.mesh.middle']))

    def view_shape_normals(self):
        """ View the shape normals at intrados and extrados surfaces

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        shape = self.shape
        if not shape:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        intrados = shape.intrados
        extrados = shape.extrados
        length = self.settings['size.edge.normals']

        for mesh in [intrados, extrados]:
            for key in mesh.vertices():
                n = mesh.vertex_attribute(key, 'n')
                pt0 = mesh.vertex_coordinates(key)
                pt1 = [pt0[0] + n[0]*length, pt0[1] + n[1]*length, pt0[2] + n[2]*length]
                line = Line(pt0, pt1)
                self.app.add(line, name='normal-{}'.format(key))

    def view_mesh(self, mesh=None, show_edges=True, opacity=0.5):
        """Add a mesh to the viewer, if no mesh is given the ``self.thrust`` is taken

        Parameters
        ----------
        mesh : Mesh, optional
            Mesh to plot, by default None
        show_edges : bool, optional
            Whether or not edges are shown, by default True
        opacity : float, optional
            The opacity of the mesh, by default 0.5

        Returns
        -------
        None
            The Viewer is updated in place.
        """

        if not mesh:
            mesh = self.thrust

        vertices_mesh, faces_mesh = mesh.to_vertices_and_faces()
        mesh = Mesh.from_vertices_and_faces(vertices_mesh, faces_mesh)

        self.app.add(mesh, show_edges=show_edges, opacity=opacity)

    def view_reactions(self):
        """View the reaction vectors on the supports according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if self.settings['show.reactions']:
            reaction_scale = self.settings['scale.reactions']

            for key in self.thrust.vertices_where({'is_fixed': True}):
                x, y, z = self.thrust.vertex_coordinates(key)
                rx = self.thrust.vertex_attribute(key, '_rx') * reaction_scale
                ry = self.thrust.vertex_attribute(key, '_ry') * reaction_scale
                rz = self.thrust.vertex_attribute(key, '_rz') * reaction_scale
                arrow = Arrow([x, y, z], [-rx, -ry, -rz], head_width=self.settings['size.reaction.head_width'], body_width=self.settings['size.reaction.body_width'])
                self.app.add(arrow, color=_norm(self.settings['color.edges.reactions']))

    def view_loads(self):
        """View the externally applied loadss as vectors on the applied nodes

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if self.settings['show.reactions']:
            reaction_scale = self.settings['scale.loads']

            for key in self.thrust.vertices():
                x, y, z = self.thrust.vertex_coordinates(key)
                px = self.thrust.vertex_attribute(key, 'px') * reaction_scale
                py = self.thrust.vertex_attribute(key, 'py') * reaction_scale
                pz = self.thrust.vertex_attribute(key, 'pz') * reaction_scale
                arrow = Arrow([x - px, y - py, z - pz], [px, py, pz], head_width=self.settings['size.reaction.head_width'], body_width=self.settings['size.reaction.body_width'])
                self.app.add(arrow, color=_norm(self.settings['color.edges.reactions']))

    def view_reaction_label(self):
        """View the reaction labels (force magnitude) on the supports according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if self.settings['show.reactionlabels']:
            reaction_scale = self.settings['scale.reactions']

            for key in self.thrust.vertices_where({'is_fixed': True}):
                x, y, z = self.thrust.vertex_coordinates(key)
                rx = - self.thrust.vertex_attribute(key, '_rx') * reaction_scale
                ry = - self.thrust.vertex_attribute(key, '_ry') * reaction_scale
                rz = - self.thrust.vertex_attribute(key, '_rz') * reaction_scale
                pt = [x + rx, y + ry, z + rz]
                reaction = '({0:3g}, {1:3g}, {2:3g})'.format(rx, ry, rz)
                text = Text(reaction, pt, height=self.settings['size.reactionlabel'])
                self.app.add(text)

    def view_dome_shape(self, radius=5.0, thickness=0.50, center=[5.0, 5.0, 0.0], spr_angle=0.0):
        """ 3D display of dome with high discretisation

        Parameters
        ----------
        radius : float, optional
            The radius of the dome, by default 5.0
        thickness : float, optional
            The thickness of the dome, by default 0.50
        center : list, optional
            List with the coordinated of the dome, by default [5.0, 5.0, 0.0]
        spr_angle : float, optional
            The springing angle of the dome, by default 0

        Returns
        -------
        [type]
            [description]
        """

        # WIP

        return

    def view_force(self, force=None, show_edges=True, opacity=0.5):
        """Add the ForceDiagram to the viewer, if no force is given the force diagram is recomputed in place.

        Parameters
        ----------
        force : ForceDiagram, optional
            ForceDiagram to plot, by default ``None``
        show_edges : bool, optional
            Whether or not edges are shown, by default ``True``
        opacity : float, optional
            The opacity of the mesh, by default 0.5

        Returns
        -------
        None
            The Viewer is updated in place.
        """

        if not force:
            from compas_tno.algorithms import reciprocal_from_form  # here to avoid double imports
            force = reciprocal_from_form(self.thrust)

        from compas.geometry import Translation
        from compas.geometry import Scale

        scale = self.settings['force.scale']
        translation = self.settings['force.anchor']

        S = Scale.from_factors(3 * [scale])
        T = Translation.from_vector(translation)

        force = force.transformed(S)
        force = force.transformed(T)

        self.app.add(force, show_edges=show_edges, opacity=opacity)


def _norm(rgb):
    """ Normalise the color in RGB dividing each component by 255

    Parameters
    ----------
    rgb : tuple
        Tuple with RBG colors

    Returns
    ----------
    norm_rgb : tuple
        The normalised color
    """
    r, g, b = rgb

    return (r/255.0, g/255.0, b/255.0)
