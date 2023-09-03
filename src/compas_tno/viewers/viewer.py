from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import Point
from math import radians
from math import sqrt

import json

from compas_tno.shapes import Shape

from compas_view2 import app
from compas_view2.shapes import Arrow
from compas_view2.shapes import Text
from compas_view2.objects import LineObject
from compas_view2.objects import PointObject

from compas.colors import Color


class Viewer(object):
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

    def __init__(self, thrust=None, shape=None, force=None, show_grid=False, **kwargs):

        super().__init__()  # **kwargs
        self.title = 'Viewer'
        self.app = None
        self.thrust = thrust
        self.shape = shape
        self.force = force
        # self.lines = []
        # self.points = []
        # self.vectors = []
        self.settings = {
            'show.thrust': True,
            'show.shape': True,
            'show.reactions': True,
            'show.reactionlabels': True,
            'show.cracks': True,
            'show.vertex.outside': True,

            'camera.target': [5, 5, 0],
            'camera.distance': 35,  # 20,
            'camera.rx': 45,
            'camera.ry': 0,
            'camera.rz': 45,
            'camera.fov': 40,
            'camera.show.grid': show_grid,
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
            'scale.loads': 0.01,
            'scale.edge.thk_absolute': 1.0,
            'opacity.shapes': 0.3,  # before 0.5

            'color.edges.thrust': Color.red(),
            'color.edges.reactions': Color.grey().lightened(80),
            'color.vertex.extrados': Color.from_rgb255(0, 125, 0),
            'color.vertex.intrados': Color.from_rgb255(0, 0, 255),
            'color.vertex.outside': Color.from_rgb255(200, 200, 200),
            'color.mesh.intrados': Color.from_rgb255(125, 125, 125),
            'color.mesh.extrados': Color.from_rgb255(125, 125, 125),
            'color.mesh.middle': Color.from_rgb255(125, 125, 125),
            'color.mesh.general': Color.from_rgb255(0, 0, 0),

            'tol.forces': 1e-2,
            'tol.cracks': 1e-3,
        }

    def initiate_app(self, **kwargs):
        """ Initiate the App with the default camera options

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app = app.App(width=1600, height=900, show_grid=self.settings['camera.show.grid'], **kwargs)

        self.app.view.camera.target = self.settings['camera.target']
        self.app.view.camera.distance = self.settings['camera.distance']
        # self.app.view.camera.rx = radians(self.settings['camera.rx'])
        # self.app.view.camera.rz = radians(self.settings['camera.rz'])
        self.app.view.camera.rotation = [radians(self.settings['camera.rx']),
                                         radians(self.settings['camera.ry']),
                                         radians(self.settings['camera.rz'])]
        self.app.view.camera.fov = self.settings['camera.fov']

    def set_camera(self, target=[5, 5, 0], distance=30, rotation=[45, 0, 45], fov=40, **kwargs):
        """ Set the camera options (can also be set directly modifying the settings dictionary)

        Returns
        -------
        None
            The objects are updated in place
        """

        self.settings['camera.target'] = target
        self.settings['camera.distance'] = distance
        self.settings['camera.rx'] = rotation[0]
        self.settings['camera.ry'] = rotation[1]
        self.settings['camera.rz'] = rotation[2]
        self.settings['camera.fov'] = fov

        self.initiate_app()

    def show_solution(self, **kwargs):
        """ Show the thrust network, with the shape according to the settings

        Returns
        -------
        None
            The objects are updated in place
        """

        if not self.app:
            self.initiate_app()

        self.draw_thrust()
        self.draw_cracks()
        self.draw_shape()
        self.draw_reactions()
        self.show()

    def show(self):
        """ Display the viewer in the screen

        Returns
        -------
        None
            The objects are updated in place
        """

        if not self.app:
            raise ValueError('First Initiate the app')

        self.app.show()

    def save(self, path='temp/fig.png', **kwargs):
        """ Save to a path

        Returns
        -------
        None
            The objects are updated in place
        """

        if not self.app:
            raise ValueError('First Initiate the app')

        # self.app.save(path)
        self.app.view.init()
        self.app.record = True
        self.app.view.paint()
        print(self.app.recorded_frames)

    def add(self, item, **kwargs):
        """Add an item to the viewer.

        Parameters
        ----------
        item : Any
            Item to add

        Returns
        -------
        None
            Viewer updated in place

        """

        if not self.app:
            self.initiate_app()

        self.app.add(item, **kwargs)

    def clear(self):
        """Clear the viewer elements

        Returns
        -------
        None
            The objects are updated in place
        """

        if not self.app:
            raise ValueError('First Initiate the app')

        self.app.view.objects = {}

    def scale_edge_thickness(self, factor=1.0):

        self.settings['scale.edge.thk_absolute'] = self.settings['scale.edge.thk_absolute'] * factor

    def draw_form(self, scale_width=True, absolute_scale=True, cull_negative=False, edges=None, line_color=None, **kwargs):
        """Alias to draw thrust network according to the settings

        Parameters
        ----------
        scale_width : bool, optional
            If the lines of the form diagram should be scaled with regards to the force carried, by default True
        absolute_scale : bool, optional
            If the axial forces should be scaled according to an absoulte scale, by default False

        Returns
        -------
        None
            The viewer is updated in place
        """

        self.draw_thrust(scale_width=scale_width, absolute_scale=absolute_scale, cull_negative=cull_negative, edges=edges, line_color=line_color, **kwargs)

    def draw_thrust(self, scale_width=True, absolute_scale=True, cull_negative=False, edges=None, line_color=None, **kwargs):
        """Draw thrust network according to the settings

        Parameters
        ----------
        scale_width : bool, optional
            If the lines of the form diagram should be scaled with regards to the force carried, by default True
        absolute_scale : bool, optional
            If the axial forces should be scaled according to an absoulte scale, by default False

        Returns
        -------
        None
            The viewer is updated in place
        """

        if not self.app:
            self.initiate_app()
        edges = edges or list(self.thrust.edges_where({'_is_edge': True}))

        base_thick = self.settings['size.edge.base_thickness']
        max_thick = self.settings['size.edge.max_thickness']
        thickness_abs = self.settings['scale.edge.thk_absolute']

        if scale_width:
            forces = [self.thrust.edge_attribute((u, v), 'q') * self.thrust.edge_length(u, v) for u, v in self.thrust.edges_where({'_is_edge': True})]
            fmax = sqrt(max(abs(max(forces)), abs(min(forces))))  # trying sqrt

        if not line_color:
            line_color = {(u, v): self.settings['color.edges.thrust'] for u, v in self.thrust.edges()}

        thks = []
        forces = []

        for u, v in edges:
            Xu = self.thrust.vertex_coordinates(u)
            Xv = self.thrust.vertex_coordinates(v)
            line = Line(Xu, Xv)
            if cull_negative:
                if min(Xu[2], Xv[2]) < 0:
                    if Xu[2] > Xv[2]:
                        Xa, Xb = Xu, Xv
                    else:
                        Xa, Xb = Xv, Xu
                    coef = - Xb[2]/(Xa[2] - Xb[2])
                    Xi = [Xb[0] + (Xa[0] - Xb[0]) * coef, Xb[1] + (Xa[1] - Xb[1]) * coef, 0.0]
                    line = Line(Xa, Xi)
            if not scale_width:
                self.app.add(line, name=str((u, v)), linewidth=base_thick, linecolor=line_color[(u, v)])
                continue
            q = self.thrust.edge_attribute((u, v), 'q')
            length = self.thrust.edge_length(u, v)
            force = sqrt(abs(q*length))
            thk = force/fmax * max_thick
            if absolute_scale:
                thk = force * thickness_abs
            forces.append(force)
            thks.append(thk)
            if force > self.settings['tol.forces'] * 2:
                # print('thk:', thk)
                self.app.add(line, name=str((u, v)), linewidth=thk, linecolor=line_color[(u, v)])

        # print('Min / Max thks:', min(thks), max(thks))
        # print('Min / Max forces:', min(forces), max(forces))

    def draw_cracks(self, cull_negative=False, points=None, **kwargs):
        """Draw cracks according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        intrad = 1
        extrad = 1
        out = 1
        points = points or list(self.thrust.vertices())

        for key in points:
            lb = self.thrust.vertex_attribute(key, 'lb')
            ub = self.thrust.vertex_attribute(key, 'ub')
            x, y, z = self.thrust.vertex_coordinates(key)
            if cull_negative:
                if z < 0.0:
                    continue
            if self.settings['show.cracks']:
                if abs(ub - z) < self.settings['tol.cracks']:
                    self.app.add(Point(x, y, z), name="Extrados (%s)" % extrad, pointcolor=self.settings['color.vertex.extrados'], pointsize=self.settings['size.vertex'])
                    extrad += 1
                elif abs(lb - z) < self.settings['tol.cracks']:
                    self.app.add(Point(x, y, z), name="Intrados (%s)" % intrad, pointcolor=self.settings['color.vertex.intrados'], pointsize=self.settings['size.vertex'])
                    intrad += 1
                if self.thrust.vertex_attribute(key, 'is_fixed'):
                    if z > 0.0:
                        if self.thrust.vertex_attribute(key, 'b'):
                            bx, by = self.thrust.vertex_attribute(key, 'b')
                            rx, ry, rz = self.thrust.vertex_attribute(key, '_rx'), self.thrust.vertex_attribute(key, '_ry'), self.thrust.vertex_attribute(key, '_rz')
                            scale = z/rz
                            if abs(rx) < self.settings['tol.cracks'] and abs(ry) < self.settings['tol.cracks']:
                                continue
                            # print('Reaction', key, rx, ry, rz)
                            crack_Rx = abs(bx) - abs(rx * scale)  # must be >= 0
                            crack_Ry = abs(by) - abs(ry * scale)  # must be >= 0

                            # print('Cracks:', key, crack_Rx, crack_Ry)

                            if crack_Rx > self.settings['tol.cracks'] and crack_Ry > self.settings['tol.cracks']:  # it is not a crack
                                continue
                            crack_point = Point(x - rx * scale, y - ry * scale, 0.0)
                            if abs(abs(bx) - abs(rx * scale)) < self.settings['tol.cracks'] or abs(abs(by) - abs(ry * scale)) < self.settings['tol.cracks']:
                                self.app.add(crack_point, name="Extrados (%s)" % extrad,
                                             pointcolor=self.settings['color.vertex.extrados'], pointsize=self.settings['size.vertex'])
                                extrad += 1
                            else:
                                self.app.add(Point(x, y, z), name="Outside - Intra (%s)" % out, pointcolor=self.settings['color.vertex.outside'],
                                             pointsize=self.settings['size.vertex'])
                                out += 1

            if self.settings['show.vertex.outside']:
                if z > ub:
                    self.app.add(Point(x, y, z), name="Outside - Intra (%s)" % out, pointcolor=self.settings['color.vertex.outside'], pointsize=self.settings['size.vertex'])
                    out += 1
                elif z < lb:
                    self.app.add(Point(x, y, z), name="Outside - Extra (%s)" % out, pointcolor=self.settings['color.vertex.outside'], pointsize=self.settings['size.vertex'])
                    out += 1

    def draw_shape(self, show_lines=False, **kwargs):
        """Draw the shape (intrados + extrados) according to the settings

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        if not self.app:
            self.initiate_app()

        if self.shape:
            datashape = self.shape.datashape.copy()
            if datashape['type'] in ['dome']:
                datashape['type'] = 'dome_polar'
                datashape['discretisation'] = [50, 50]
                shape = Shape.from_library(datashape)
                print('Drawing nicer dome')
            elif datashape['type'] == 'arch':
                datashape['type'] = 'arch_polar'
                shape = Shape.from_library(datashape)
                print('Drawing nicer arch')
            else:
                shape = self.shape

        else:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        self.app.add(shape.intrados, name="Intrados", show_lines=show_lines, opacity=self.settings['opacity.shapes'], facecolor=self.settings['color.mesh.intrados'])
        self.app.add(shape.extrados, name="Extrados", show_lines=show_lines, opacity=self.settings['opacity.shapes'], facecolor=self.settings['color.mesh.extrados'])

    def draw_middle_shape(self, **kwargs):
        """ Draw the middle of the shape according to the settings

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        if not self.app:
            self.initiate_app()

        shape = self.shape
        if not shape:
            shape = Shape.from_formdiagram_and_attributes(self.thrust)

        vertices_middle, faces_middle = shape.middle.to_vertices_and_faces()
        mesh_middle = Mesh.from_vertices_and_faces(vertices_middle, faces_middle)
        self.app.add(mesh_middle, name="Middle", show_lines=True, opacity=self.settings['opacity.shapes'], facecolor=self.settings['color.mesh.middle'])

    def draw_shape_normals(self, **kwargs):
        """ Draw the shape normals at intrados and extrados surfaces

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        if not self.app:
            self.initiate_app()

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

    def draw_mesh(self, mesh=None, show_lines=True, opacity=0.5, color=Color.grey(), **kwargs):
        """Draw a mesh to the viewer, if no mesh is given the ``self.thrust`` is taken

        Parameters
        ----------
        mesh : Mesh, optional
            Mesh to plot, by default None
        show_lines : bool, optional
            Whether or not edges are shown, by default True
        opacity : float, optional
            The opacity of the mesh, by default 0.5

        Returns
        -------
        None
            The Viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        if not mesh:
            mesh = self.thrust

        vertices_mesh, faces_mesh = mesh.to_vertices_and_faces()
        mesh = Mesh.from_vertices_and_faces(vertices_mesh, faces_mesh)

        self.app.add(mesh, show_lines=show_lines, opacity=opacity, facecolor=color)

    def draw_thrustsurface(self, show_edges=False, show_faces=True, opacity=0.2, color=Color.from_rgb255(125, 125, 125), **kwargs):
        """Draw a mesh to the viewer, if no mesh is given the ``self.thrust`` is taken

        Parameters
        ----------
        mesh : Mesh, optional
            Mesh to plot, by default None
        show_faces : bool, optional
            Whether or not faces are shown, by default True
        show_edges : bool, optional
            Whether or not edges are shown, by default False
        opacity : float, optional
            The opacity of the mesh, by default 0.5
        color : Color, optional
            Color of the mesh, by default grey

        Returns
        -------
        None
            The Viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        mesh = self.thrust

        vertices_mesh, faces_mesh = mesh.to_vertices_and_faces()
        mesh = Mesh.from_vertices_and_faces(vertices_mesh, faces_mesh)

        self.app.add(mesh, show_lines=show_edges, show_faces=show_faces, opacity=opacity, facecolor=color)

    def draw_reactions(self, emerging_reactions=False, supports=None, extend_reactions=False, scale=None, **kwargs):
        """Draw the reaction vectors on the supports according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        if self.settings['show.reactions']:
            reaction_scale = scale or self.settings['scale.reactions']
            thickness_scale = self.settings['scale.edge.thk_absolute']

            supports = supports or list(self.thrust.vertices_where({'is_fixed': True}))

            for key in supports:
                x, y, z = self.thrust.vertex_coordinates(key)
                rx = self.thrust.vertex_attribute(key, '_rx') * reaction_scale
                ry = self.thrust.vertex_attribute(key, '_ry') * reaction_scale
                rz = self.thrust.vertex_attribute(key, '_rz') * reaction_scale
                r = 1/reaction_scale * (rx**2 + ry**2 + rz**2)**(1/2)
                if r < 1e-3:
                    continue
                if extend_reactions:
                    if z < 0:
                        continue
                    rs = z / abs(self.thrust.vertex_attribute(key, '_rz'))
                    rx, ry, rz = self.thrust.vertex_attribute(key, '_rx') * rs, self.thrust.vertex_attribute(key, '_ry') * rs, self.thrust.vertex_attribute(key, '_rz') * rs,
                    max_thickness = self.settings['size.edge.max_thickness']
                    max_f = max([abs(self.thrust.edge_attribute(edge, 'q') * self.thrust.edge_length(*edge)) for edge in self.thrust.edges()])
                    thickness = r/max_f*max_thickness

                    thickness = thickness_scale*sqrt(r)

                    line = Line([x, y, z], [x - rx, y - ry, z - rz])
                    self.app.add(line, linewidth=thickness, linecolor=self.settings['color.edges.thrust'])
                    continue
                if emerging_reactions:
                    arrow = Arrow([x, y, z], [-rx, -ry, -rz], head_width=self.settings['size.reaction.head_width'],
                                  body_width=self.settings['size.reaction.body_width'])
                else:
                    arrow = Arrow([x - rx, y - ry, z - rz], [rx, ry, rz], head_width=self.settings['size.reaction.head_width'],
                                  body_width=self.settings['size.reaction.body_width'])
                self.app.add(arrow, linecolor=self.settings['color.edges.reactions'])  # check if this color= works for new compasview2

    def draw_vector(self, vector=None, base=None, **kwargs):
        """Helper to add vectors to the plotter.

        Parameters
        ----------
        vectors : Vector, optional
            Vector to plot, by default None
        bases : Point, optional
            Point with the location of the base of the vector, by default None

        Returns
        -------
        None
            The plotter is updated in place
        """

        if not vector or not base:
            raise ValueError()

        arrow = Arrow(base, vector, head_width=self.settings['size.reaction.head_width'], body_width=self.settings['size.reaction.body_width'])
        self.app.add(arrow, **kwargs)

    def draw_loads(self, **kwargs):
        """Draw the externally applied loadss as vectors on the applied nodes

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        if self.settings['show.reactions']:
            reaction_scale = self.settings['scale.loads']

            for key in self.thrust.vertices():
                x, y, z = self.thrust.vertex_coordinates(key)
                px = self.thrust.vertex_attribute(key, 'px') * reaction_scale
                py = self.thrust.vertex_attribute(key, 'py') * reaction_scale
                pz = self.thrust.vertex_attribute(key, 'pz') * reaction_scale
                arrow = Arrow([x - px, y - py, z - pz], [px, py, pz], head_width=self.settings['size.reaction.head_width'], body_width=self.settings['size.reaction.body_width'])
                self.app.add(arrow, linecolor=self.settings['color.edges.reactions'])

    def draw_reaction_label(self, **kwargs):
        """Draw the reaction labels (force magnitude) on the supports according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

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

    def draw_force(self, force=None, show_edges=True, opacity=0.5, **kwargs):
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

        if not self.app:
            self.initiate_app()

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

        self.app.add(force, show_lines=show_edges, opacity=opacity)

    def draw_assembly(self, assembly, **kwargs):
        """Draw the reaction labels (force magnitude) on the supports according to the settings

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        for block in assembly.blocks():
            self.add(block, **kwargs)

    def draw_points(self, points=[], **kwargs):
        """Draw points in space

        Returns
        -------
        None
            The viewer is updated in place.
        """

        if not self.app:
            self.initiate_app()

        for point in points:
            self.add(Point(*point), **kwargs)

    def draw_b_constraint(self):

        base_thick = self.settings['size.edge.base_thickness']
        for key in self.thrust.vertices_where({'is_fixed': True}):
            x, y, _ = self.thrust.vertex_coordinates(key)
            b = self.thrust.vertex_attribute(key, 'b')
            if b:
                bx, by = b
                if abs(bx) > 1e-4:
                    line = Line([x, y, 0], [x + bx, y + by, 0])
                    self.app.add(line, linewidth=base_thick)
                if abs(by) > 1e-4:
                    line = Line([x, y, 0], [x + bx, y + by, 0])
                    self.app.add(line)

    def to_objects(self):
        """Generate General objects (dicts) from the added elements in the viewer (WIP).

        Returns
        -------
        data
            Raw geometry objects in a .
        """

        objects = self.app.view.objects

        print(objects)

        data = {}
        lines = []
        points = []

        for obj in objects:
            if isinstance(obj, LineObject):
                lines.append({'start': [obj._data[0].x, obj._data[0].y, obj._data[0].z],
                              'end': [obj._data[1].x, obj._data[1].y, obj._data[1].z],
                              'color': (255*obj.linecolor).tolist(),
                              'width': obj.linewidth
                              })
            if isinstance(obj, PointObject):
                points.append({'pos': [obj._data.x, obj._data.y, obj._data.z],
                               'color': (255*obj.pointcolor).tolist()
                               })

        data['lines'] = lines
        data['points'] = points

        return data

    def to_json(self, file="out.json"):

        data = self.to_objects()

        with open(file, "w") as outfile:
            json.dump(data, outfile)
