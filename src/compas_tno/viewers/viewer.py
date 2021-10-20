from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import Point

from compas_tno.shapes import Shape

from compas_view2 import app
from compas_view2.shapes import Arrow
from compas_view2.shapes import Text

from compas_tno.algorithms import reactions


__all__ = ['Viewer']


class Viewer(object):
    """A Class for view 3D thrust networks and shapes.
    """

    def __init__(self, thrust=None, shape=None, **kwargs):
        super().__init__(**kwargs)
        self.title = 'Viewer'
        self.app = None
        self.thrust = thrust
        self.shape = shape
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

            'size.vertex': 15.0,
            'size.edge.max_thickness': 10.0,
            'size.edge.base_thickness': 1.0,
            'size.edge.normals': 0.5,
            'size.reaction.head_width': 0.1,
            'size.reaction.body_width': 0.01,
            'size.reactionlabel': 20,

            'scale.reactions': 0.0005,
            'opacity.shapes': 0.5,

            'color.edges.thrust': (255, 0, 0),
            'color.edges.reactions': (125, 125, 125),
            'color.vertex.extrados': (0, 125, 0),
            'color.vertex.intrados': (0, 0, 255),
            'color.vertex.outside': (200, 200, 200),
            'color.mesh.intrados': (125, 125, 125),
            'color.mesh.extrados': (125, 125, 125),
            'color.mesh.middle': (125, 125, 125),
            'color.mesh.general': (125, 125, 125),

            'tol.forces': 1e-3,
        }
        self.initiate_app()

    def initiate_app(self):
        """ Initiate the App with the default camera options """

        self.app = app.App()

        self.app.view.camera.target = self.settings['camera.target']
        self.app.view.camera.distance = self.settings['camera.distance']
        self.app.view.camera.rx = self.settings['camera.rx']
        self.app.view.camera.rz = self.settings['camera.rz']
        self.app.view.camera.fov = self.settings['camera.fov']

    def show_solution(self):
        """ Show the thrust network, with the shape according to the settings """

        self.view_thrust()
        self.view_cracks()
        self.view_shape()
        self.view_reactions()
        self.show()

    def show(self):
        """ Display the viewer in the screen """

        self.app.show()

    def view_thrust(self):
        """ View thrust network according to the settings """

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
                self.app.add(line, name=str((u, v)), linewidth=base_thick, color=(255, 0, 0))
                continue
            q = self.thrust.edge_attribute((u, v), 'q')
            length = self.thrust.edge_length(u, v)
            force = abs(q*length)
            thk = force/fmax * max_thick
            if force > 10e-4:
                self.app.add(line, name=str((u, v)), linewidth=thk, color=(255, 0, 0))

    def view_cracks(self):
        """ View cracks according to the settings """

        intrad = 1
        extrad = 1
        out = 1

        for key in self.thrust.vertices():
            lb = self.thrust.vertex_attribute(key, 'lb')
            ub = self.thrust.vertex_attribute(key, 'ub')
            x, y, z = self.thrust.vertex_coordinates(key)
            if self.settings['show.cracks']:
                if abs(ub - z) < 10e-4:
                    self.app.add(Point(x, y, z), name="Extrados (%s)" % extrad, color=_norm(self.settings['color.vertex.extrados']), size=self.settings['size.vertex'])
                    extrad += 1
                elif abs(lb - z) < 10e-4:
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
        """ View the shape (intrados + extrados) according to the settings """

        shape = self.shape
        if not shape:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        vertices_intra, faces_intra = shape.intrados.to_vertices_and_faces()
        mesh_intra = Mesh.from_vertices_and_faces(vertices_intra, faces_intra)

        vertices_extra, faces_extra = shape.extrados.to_vertices_and_faces()
        mesh_extra = Mesh.from_vertices_and_faces(vertices_extra, faces_extra)

        self.app.add(mesh_intra, name="Intrados", show_edges=False, opacity=self.settings['opacity.shapes'], color=_norm(self.settings['color.mesh.intrados']))
        self.app.add(mesh_extra, name="Extrados", show_edges=False, opacity=self.settings['opacity.shapes'], color=_norm(self.settings['color.mesh.extrados']))

    def view_middle_shape(self):
        """ View the middle of the shape according to the settings """

        shape = self.shape
        if not shape:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        vertices_middle, faces_middle = shape.middle.to_vertices_and_faces()
        mesh_middle = Mesh.from_vertices_and_faces(vertices_middle, faces_middle)
        self.app.add(mesh_middle, name="Middle", show_edges=False, opacity=0.5, color=(0.5, 0.5, 0.5))

    def view_shape_normals(self):
        """ View the shape normals at intrados and extrados surfaces """

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
                self.app(line, name='normal-{}'.format(key))

    def view_mesh(self, mesh=None):
        """ Add a mesh to the viewer, if no mesh is given the ``self.thrust`` is taken """

        if not mesh:
            mesh = self.thrust

        vertices_mesh, faces_mesh = mesh.to_vertices_and_faces()
        mesh = Mesh.from_vertices_and_faces(vertices_mesh, faces_mesh)

        self.app.add(mesh, show_edges=True, color=_norm(self.settings['color.mesh.general']))

    def view_reactions(self):
        """ View the reaction vectors on the supports according to the settings """

        if self.settings['show.reactions']:
            reaction_scale = self.settings['scale.reactions']

            for key in self.thrust.vertices_where({'is_fixed': True}):
                x, y, z = self.thrust.vertex_coordinates(key)
                rx = self.thrust.vertex_attribute(key, '_rx') * reaction_scale
                ry = self.thrust.vertex_attribute(key, '_ry') * reaction_scale
                rz = self.thrust.vertex_attribute(key, '_rz') * reaction_scale
                arrow = Arrow([x, y, z], [-rx, -ry, -rz], head_width=self.settings['size.reaction.head_width'], body_width=self.settings['size.reaction.body_width'])
                self.app.add(arrow, color=_norm(self.settings['color.edges.reactions']))

    def view_reaction_label(self):
        """ View the reaction labels (force magnitude) on the supports according to the settings """

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

    def view_reaction_label(self):
        """ View the reaction labels (force magnitude) on the supports according to the settings """

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
