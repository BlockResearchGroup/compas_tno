from compas_plotters import Plotter
from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import Rotation
from compas.geometry import Translation
from compas.geometry import Scale
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Polygon
from compas_tno.shapes import Shape

from compas.colors import Color

from compas.utilities import rgb_to_hex
import matplotlib.pyplot as plt
import math


class TNOPlotter(object):
    """TNO helper to plot forms and shapes

    Parameters
    ----------
    form : FormDiagram, optional
        The form diagram to be plotted, by default None
    shape : Shape, optional
        The shape of the masonry constraining the solution, by default None
    form_base : FormDiagram, optional
        The base form diagram to be visualised if pattern is allowed to move, by default None
    force : ForceDiagram, optional
        The force diagram associates with the form diagram, by default None
    figsize : tuple, optional
        The size of the plot, by default (8, 8)

    Attributes
    ----------

    None

    """

    def __init__(self, form=None,
                 shape=None,
                 form_base=None,
                 force=None,
                 figsize=(8, 8),
                 *args,
                 **kwargs):

        super().__init__(*args, **kwargs)
        self.title = 'Plotter'
        self.app = None
        self._form = form
        self._shape = shape
        self._form_base = form_base
        self._force = force
        self._formartist = None
        self._otherartists = []
        self.scale_forcediagram = 1.0
        self.translation_forcediagram = [0.0, 0.0, 0.0]
        self.settings = {
            'show.thrust': True,
            'show.shape': True,
            'show.reactions': True,
            'show.reactions.emerging': True,
            'show.reactionlabels': True,
            'show.supports': True,
            'show.edges': True,
            'show.faces': False,
            'show.cracks': True,
            'show.vertex.outside': True,
            'show.edge.thickness': True,
            'show.reactions.extended': False,
            'show.reactions.asarrows': True,

            'camera.target': [0, 0, 0],
            'camera.distance': 40,
            'camera.rx': -45,
            'camera.rz': 45,
            'camera.fov': 40,
            'camera.show.grid': True,
            'camera.show.axis': True,

            'size.vertex': 5.0,
            'size.edge.max_thickness': 10.0,
            'size.edge.base_thickness': 1.0,
            'size.edge.normals': 0.5,
            'size.reaction.head_width': 0.1,
            'size.reaction.body_width': 0.01,
            'size.reactionlabel': 20,

            'scale.reactions': 0.05,
            'scale.loads': 0.001,
            'shape.opacity': 0.5,
            'rotated': False,

            'color.edges.form': Color.from_rgb255(255, 0, 0),
            'color.edges.reactions': Color.from_rgb255(170, 170, 170),
            'color.edges.shape': Color.from_rgb255(200, 200, 200),
            'color.edges.form_base': Color.from_rgb255(200, 200, 200),
            'color.edges.independent': Color.from_rgb255(0, 0, 255),
            'color.vertex.supports': Color.from_rgb255(170, 170, 170),
            'color.vertex.extrados': Color.from_rgb255(0, 125, 0),
            'color.vertex.intrados': Color.from_rgb255(0, 0, 255),
            'color.vertex.outside': Color.from_rgb255(200, 200, 200),
            'color.mesh.intrados': Color.from_rgb255(125, 125, 125),
            'color.mesh.extrados': Color.from_rgb255(125, 125, 125),
            'color.mesh.middle': Color.from_rgb255(125, 125, 125),
            'color.mesh.general': Color.from_rgb255(0, 0, 0),

            'tol.forces': 1e-3,
        }
        self.initiate_app(figsize=figsize)

    @property
    def form(self):
        """The FormDiagram of the Plotter."""
        return self._form

    @form.setter
    def form(self, value):
        self._form = value

    @property
    def shape(self):
        """The Shape of the Plotter."""
        return self._shape

    @shape.setter
    def shape(self, value):
        self._shape = value

    @property
    def form_base(self):
        """The form_base property."""
        return self._form_base

    @form_base.setter
    def form_base(self, value):
        self._form_base = value

    @property
    def force(self):
        """The force property."""
        return self._force

    @force.setter
    def force(self, value):
        self._force = value

    @property
    def formartist(self):
        """The formartist property."""
        return self._formartist

    @formartist.setter
    def formartist(self, value):
        if not self._formartist:
            self._formartist = value

    def initiate_app(self, figsize=(8, 8)):
        """Initiate the Plotter with the default camera options.

        Parameters
        ----------
        figsize : tuple, optional
            Size of the plot, by default (8, 8)

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.app = Plotter(figsize=figsize)
        self.app.fontsize = 10

    def show_solution(self):
        """Show the thrust network, with the shape according to the settings.

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.draw_form()
        self.draw_supports()
        self.draw_cracks()
        self.zoom_extents()
        self.show()

    def zoom_extents(self, padding=0.1):
        """Wrapper to extend the objects in the active view.

        Parameters
        ----------
        padding : int, optional
            Additional white space surrounding the object, by default 0.1.

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.app.zoom_extents(padding=padding)

    def show(self, padding=0.1):
        """Display the plotter in the screen.

        Parameters
        ----------
        padding : int, optional
            Additional white space surrounding the object, by default 0.1.

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.app.zoom_extents(padding=padding)
        self.app.show()

    def clear(self):
        """Clear the plotter elements.

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.app.clear()

    def save(self, filepath):
        """Save the plotter window as an image.

        Parameters
        ----------
        filepath : str
            The file path to save the image.

        Returns
        -------
        None
            Picture is saved

        """

        # self.zoom_extents()

        self.app.viewbox = ([10, 30], [-1, 10])
        self.app.axes.set_xlim([10, 30])
        self.app.axes.set_ylim([-1, 10])
        self.app.axes.autoscale_view()

        self.app.save(filepath)

    def add(self, item, **kwargs):
        """Add an item to the plotter.

        Parameters
        ----------
        item : Any
            Item to add

        Returns
        -------
        Artist
            The artist of the added item

        """

        return self.app.add(item, **kwargs)

    def draw_form(self, scale_width=True, edges=None, color=None, **kwargs):
        """Draw the Form Diagram with or without thicknesses of edges scaled as forces in the edges.

        Parameters
        ----------
        scale_width : bool, optional
            If the lines of the form diagram should be scaled with regards to the force carried, by default True
        color : Color, optional
            Color of the edges of the form diagram, this will overwrite the setting 'color.edges.form', by default None
        edges : [tuple], optional
            Edges to consider in the plot, by default None, in which all edges are plotted

        Returns
        -------
        None
            The plotter is updated in place
        """

        base_thick = self.settings['size.edge.base_thickness']
        max_thick = self.settings['size.edge.max_thickness']
        color = color or self.settings['color.edges.form']
        edgecolor = {(u, v): color for u, v in self.form.edges()}

        if scale_width:
            forcedensities = self.form.edges_attribute('q')
            lengths = [self.form.edge_length(u, v) for u, v in self.form.edges()]
            forces = [abs(forcedensities[i] * lengths[i]) for i in range(len(lengths))]
            fmax = max(forces)
            edgewidths = {edge: forces[i] / fmax * max_thick for i, edge in enumerate(self.form.edges())}
        else:
            edgewidths = base_thick
            edgecolor = {(u, v): Color.black() for u, v in self.form.edges()}

        faces = []

        if self.settings['show.edges']:
            edges = edges or list(self.form.edges_where({'_is_edge': True}))

        if self.settings['show.faces']:
            faces = list(self.form.faces())

        formartist = self.app.add(  # change this to a bunch of lines.... Better to edit
            self.form,              # add cull_negative option to cull the thrust line below a value...
            # vertices=[],
            edges=edges,
            faces=faces,
            edgewidth=edgewidths,
            edgecolor=edgecolor,
            show_vertices=False,
            **kwargs
        )

        self.formartist = formartist

        # artist = self.app.add(
        #     self.form,
        #     edges=edges,
        #     edgewidth=width,
        #     edgetext=text,
        #     edgecolor=color,
        #     show_vertices=False,
        #     show_faces=False
        # )
        # self.formartist = artist

    def draw_cracks(self, points=None, **kwargs):
        """Adds to the basic plot, the cracks which are the points of the mesh that touch intrados or extrados.

        Points
        ------
        points : list, optional
            The list with the cracks to consider

        Returns
        -------
        None
            The plotter is updated in place
        """

        tol = self.settings['tol.forces']
        color_intra = self.settings['color.vertex.intrados']
        color_extra = self.settings['color.vertex.extrados']
        color_outside = self.settings['color.vertex.outside']

        cracks = []
        points = points or list(self.form.vertices())

        for key in points:
            x, y, z = self.form.vertex_coordinates(key)
            if self.settings['rotated']:
                x, z, y = self.form.vertex_coordinates(key)
            lb = self.form.vertex_attribute(key, 'lb')
            ub = self.form.vertex_attribute(key, 'ub')
            if not lb:
                continue
            if not ub:
                continue
            if self.settings['show.cracks']:
                if abs(ub - z) < tol:
                    cracks.append({
                        'key': key,
                        'xyz': [x, y, z],
                        'color': color_extra
                    })
                    continue
                elif abs(lb - z) < tol:
                    cracks.append({
                        'key': key,
                        'xyz': [x, y, z],
                        'color': color_intra
                    })
                    continue
            if self.settings['show.vertex.outside']:
                if z - ub > tol or lb - z > tol:
                    cracks.append({
                        'key': key,
                        'xyz': [x, y, z],
                        'color': color_outside
                    })
                    continue

        if self.settings['show.reactions.extended']:
            for key in self.form.vertices_where({'is_fixed': True}):
                x, y, z = self.form.vertex_coordinates(key)
                if self.settings['rotated']:
                    x, z, y = self.form.vertex_coordinates(key)
                b = self.form.vertex_attribute(key, 'b')
                if z > tol:
                    reaction_scale = z/self.form.vertex_attribute(key, '_rz')
                else:
                    continue
                rx = self.form.vertex_attribute(key, '_rx') * reaction_scale
                ry = self.form.vertex_attribute(key, '_ry') * reaction_scale
                rz = self.form.vertex_attribute(key, '_rz') * reaction_scale
                if abs(abs(rx) - b[0]) < tol:
                    cracks.append({
                        'key': key,
                        'xyz': [x - rx, y - ry, z - rz],
                        'color': color_extra
                    })

        for point in cracks:
            x, y, z = point['xyz']
            if self.settings['rotated']:
                x, z, y = point['xyz']
            color = point['color']
            pt = Point(x, y, z)
            pointartist = self.app.add(pt, facecolor=color, size=self.settings['size.vertex'])
            self._otherartists.append(pointartist)

    def draw_supports(self, color=None, size=None):
        """Add the supports as points in the mesh.

        Parameters
        ----------
        color : Color, optional
            Color of the supports of the form diagram, this will overwrite the setting 'color.vertex.supports', by default None
        size : float, optional
            Size of the supports, this will overwrite the setting 'size.vertex', by default None

        Returns
        -------
        None
            The plotter is updated in place
        """

        if self.settings['show.supports']:
            supportcolor = color or self.settings['color.vertex.supports']  # _norm(self.settings['color.vertex.supports'])
            size = size or self.settings['size.vertex']
            rollercolor = Color.from_rgb255(255, 165, 0)

            for key in self.form.vertices_where({'is_fixed': True}):
                x, y, z = self.form.vertex_coordinates(key)
                pt = Point(x, y, z)
                pointartist = self.app.add(pt, facecolor=supportcolor, size=size)
                self._otherartists.append(pointartist)

            rollers = list(self.form.vertices_where({'rol_x': True})) + list(self.form.vertices_where({'rol_y': True}))
            for key in rollers:
                x, y, z = self.form.vertex_coordinates(key)
                pt = Point(x, y, z)
                pointartist = self.app.add(pt, facecolor=rollercolor, size=self.settings['size.vertex'])
                self._otherartists.append(pointartist)

    def draw_reactions(self, scale_width=True):
        """Add to the plots the vector of the reaction forces.

        Parameters
        ----------
        scale_width : bool, optional
            If edges should be scaled with the force, by default True

        Returns
        -------
        None
            The plotter is updated in place
        """

        if self.settings['show.reactions']:
            reaction_scale = self.settings['scale.reactions']

            for key in self.form.vertices_where({'is_fixed': True}):
                x, y, z = self.form.vertex_coordinates(key)
                if self.settings['show.reactions.extended']:  # assume rotated
                    if y > 0:
                        reaction_scale = y/self.form.vertex_attribute(key, '_rz')
                    else:
                        reaction_scale = 0.0
                rx = self.form.vertex_attribute(key, '_rx') * reaction_scale
                ry = self.form.vertex_attribute(key, '_ry') * reaction_scale
                rz = self.form.vertex_attribute(key, '_rz') * reaction_scale
                norm = math.sqrt(self.form.vertex_attribute(key, '_rx')**2 + self.form.vertex_attribute(key, '_ry')**2 + self.form.vertex_attribute(key, '_rz')**2)
                if self.settings['rotated']:
                    rx, ry, rz = rx, rz, ry
                if self.settings['show.reactions.emerging']:
                    r = Vector(-rx, -ry, -rz)
                    pt = Point(x, y, z)
                else:
                    r = Vector(rx, ry, rz)
                    pt = Point(x - rx, y - ry, z - rz)
                if self.settings['show.reactions.asarrows']:
                    self.app.add(r, point=pt, color=self.settings['color.edges.reactions'])
                else:
                    pt1 = Point(x, y, z)
                    pt2 = Point(x - rx, y - ry, z - rz)
                    line = Line(pt1, pt2)
                    max_f = max([abs(self.form.edge_attribute(edge, 'q') * self.form.edge_length(*edge)) for edge in self.form.edges()])
                    if scale_width:
                        width = max_f/norm*self.settings['size.edge.max_thickness']
                    else:
                        width = self.settings['size.edge.base_thickness']
                    self.app.add(line,
                                 draw_as_segment=True,
                                 color=self.settings['color.edges.form'],
                                 linewidth=width)

    def draw_vector(self, vector=None, base=None):
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

        vectorartist = self.app.add(vector, point=base)
        self._otherartists.append(vectorartist)

    def draw_lines(self, lines=[]):
        """Helper to add vectors to the plotter.

        Parameters
        ----------
        lines : list, optional
            The list of lines to plot, by default []

        Returns
        -------
        None
            The plotter is updated in place
        """

        for line in lines:
            linecompas = Line(line[0], line[1])
            self.app.add(linecompas, draw_as_segment=True)

    def draw_mesh(self, mesh=None, show_edges=True, show_vertices=False, show_faces=False, **kwargs):
        """Initiate the Plotter with the default camera options.

        Parameters
        ----------
        mesh : Mesh, optional
            Mesh to plot, by default None, which means that the form diagram is added
        show_edges : bool, optional
            If edges are shown, by default True
        show_vertices : bool, optional
            If vertices are shown, by default False
        show_faces : bool, optional
            If faces are shown, by default False

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if not mesh:
            mesh = self.form
        self.app.add(mesh,
                     show_edges=show_edges,
                     show_vertices=show_vertices,
                     show_faces=show_faces,
                     **kwargs
                     )

    def draw_form_xz(self, scale_width=True, edges=None, **kwargs):
        """Plot the form diagram rotated 90 degrees

        Parameters
        ----------
        scale_width : bool, optional
            If edges should be scaled with the force, by default True
        edges : [tuple], optional
            Edges to consider in the plot, by default None, in which all edges are plotted

        Returns
        -------
        None
            The plotter is updated in place.
        """

        axis = Vector(1.0, 0, 0)
        self.form = self.form.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))
        self.settings['rotated'] = True

        self.draw_form(scale_width=scale_width, edges=edges, **kwargs)

    def draw_shape(self, update_from_parameters=True, **kwargs):
        """Adds the shape to the plot.

        Parameters
        ----------
        update_from_parameters : bool, optional
            Whether or not should look for shapes having special cases such as "arch" and "dome" to adapt.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if not self.shape:
            self.shape = Shape.from_formdiagram_and_attributes(self.form)

        intrados = self.shape.intrados
        extrados = self.shape.extrados

        if update_from_parameters:
            datashape = self.shape.datashape.copy()
            if datashape['type'] in ['arch']:
                stereotomy = False
                close_bottom = True
                total_nodes = datashape['discretisation']
                for key in kwargs:
                    if key == 'stereotomy':
                        stereotomy = kwargs[key]
                    if key == 'close_bottom':
                        close_bottom = kwargs[key]
                    if key == 'blocks':
                        total_nodes = kwargs[key]
                self.draw_arch_lines(H=datashape['H'],
                                     L=datashape['L'],
                                     x0=datashape['x0'],
                                     thk=datashape['thk'],
                                     total_nodes=total_nodes,
                                     stereotomy=stereotomy,
                                     close_bottom=close_bottom)
                return

        self.app.add(
            intrados,
            opacity=self.settings['shape.opacity'],
            edgecolor={edge: self.settings['color.edges.shape'] for edge in intrados.edges()},
            show_vertices=False
        )

        self.app.add(
            extrados,
            opacity=self.settings['shape.opacity'],
            edgecolor={edge: self.settings['color.edges.shape'] for edge in extrados.edges()},
            show_vertices=False
        )

    def draw_shape_xz(self, update_from_parameters=True, **kwargs):
        """Adds the shape to the plot rotated 90 degrees.

        Parameters
        ----------
        update_from_parameters : bool, optional
            Whether or not should look for shapes having special cases such as "arch" and "dome" to adapt.

        Returns
        -------
        None
            The plotter is updated in place
        """

        if not self.shape:
            self.shape = Shape.from_formdiagram_and_attributes(self.form)

        axis = Vector(1.0, 0, 0)
        self.shape.intrados = self.shape.intrados.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))
        self.shape.extrados = self.shape.extrados.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))
        self.settings['rotated'] = True

        self.draw_shape(update_from_parameters=update_from_parameters, **kwargs)

    def draw_base_form(self, form_base=None):
        """Adds to the plot the base mesh which is the mesh before the nodes moved horizontally.

        Parameters
        ----------
        form_base : :class:`~compas_tno.diagrams.FormDiagram`, optional
            The base to plot, by default None.

        Returns
        -------
        None
            The plotter is updated in place
        """

        if form_base:
            self.form_base = form_base

        if not self.form_base:
            return

        edgecolors = {edge: self.settings['color.edges.form_base'] for edge in self.form_base.edges()}

        if self.settings['show.edges']:
            edges = list(self.form_base.edges_where({'_is_edge': True}))

        self.app.add(
            self.form_base,
            edges=edges,
            opacity=0.5,
            edgecolor=edgecolors,
            show_vertices=False,
            show_faces=False
        )

    def draw_constraints(self):
        """ Draw the constraints applied.

        Returns
        -------
        None
            The Viewer object is modified in place
        """

        draw_xy_constraints = True

        if draw_xy_constraints:
            for key in self.form.vertices():
                x, y, _ = self.form.vertex_coordinates(key)
                xmin = self.form.vertex_attribute(key, 'xmin')
                xmax = self.form.vertex_attribute(key, 'xmax')
                ymin = self.form.vertex_attribute(key, 'ymin')
                ymax = self.form.vertex_attribute(key, 'ymax')
                if xmin:
                    p0 = [xmin, ymin]
                    p1 = [xmax, ymin]
                    p2 = [xmax, ymax]
                    p3 = [xmin, ymax]
                    lines = [Line(p0, p1), Line(p1, p2), Line(p2, p3), Line(p3, p0)]
                    for line in lines:
                        self.app.add(line, draw_as_segment=True, color=self.settings['color.edges.form_base'])

    def draw_form_independents(self, color=None, show_text=True, **kwargs):
        """Draw the form diagram with highlight in the independent edges.

        Parameters
        ----------
        show_text : bool, optional
            If text should be added to the independent edges, by default True
        color : Color, optional
            Color of the independent edges of the form diagram, this will overwrite the setting 'color.edges.independent', by default None

        Returns
        -------
        None
            The plotter is updated in place
        """

        edges = list(self.form.edges_where({'_is_edge': True}))
        base_thick = self.settings['size.edge.base_thickness']
        colorinds = color or self.settings['color.edges.independent']

        color = {}
        width = {}
        text = {}
        i = 0
        for edge in edges:
            if self.form.edge_attribute(edge, 'is_ind'):
                color[edge] = colorinds
                width[edge] = base_thick * 3.0
                if show_text:
                    text[edge] = str(i)
                i = i+1
            else:
                color[edge] = (0, 0, 0)
                width[edge] = base_thick

        artist = self.app.add(
            self.form,
            edges=edges,
            edgewidth=width,
            edgetext=text,
            edgecolor=color,
            show_vertices=False,
            show_faces=False
        )
        self.formartist = artist

    def draw_vertexlabels(self, text=None, label=None, fontsize=12, **kwargs):
        """Draw a specific label in the vertices

        Parameters
        ----------
        label : str, optional
            The vertex label to be shown, by default None, in which case the keys are displayed.
        fontsize : int, optional
            Font size, by default 12

        """

        if not self.formartist:
            return

        print(text)
        if not text:
            text = {}
            for key in self.form.vertices():
                text[key] = key
            if isinstance(label, str):
                text = {}
                for key in self.form.vertices():
                    text[key] = '{:.3}'.format(float(self.form.vertex_attribute(key, label)))
        else:
            pass

        self.app.fontsize = fontsize
        self.formartist.draw_vertexlabels(text=text)

    # def draw_vertexlabels(self, text=None, fontsize=12):
    #     """Draw Labels to the vertices

    #     Parameters
    #     ----------
    #     text : dict, optional
    #         The dictionary to print, if None, the keys are displayed, by default None
    #     fontsize : int, optional
    #         Font size, by default 12
    #     """

    #     if not self.formartist:
    #         return

    #     self.formartist.plotter.fontsize = fontsize
    #     self.formartist.draw_vertexlabels(text=text)

    def highlight_vertices(self, vertices, show_forcepolygon=False):
        """Highlight a vertex in the form diagram and in the force diagram (if they have been previously added to the plotter)

        Parameters
        ----------
        vertices : [int]
            Key of the vertices to highlight
        show_forcepolygon : bool, optional
            If force polygons should be also displayed in the force diagram, by default False
            Note: this requires that the force diagram is already added to the scene.

        Returns
        -------
        None
            The Plotter is updated in place.
        """

        points = [self.form.vertex_coordinates(vertex) for vertex in vertices]
        edges = [self.form.vertex_edges(vertex) for vertex in vertices]
        base_thick = self.settings['size.edge.base_thickness']

        for i, point in enumerate(points):
            pt = Point(point[0], point[1], point[2])
            self.app.add(pt, facecolor=(0.8, 0.8, 0.8), size=self.settings['size.vertex'])

            for edge in edges[i]:
                line = Line(self.form.vertex_coordinates(edge[0]), self.form.vertex_coordinates(edge[1]))
                self.app.add(line, color=(0, 0, 0), linewidth=base_thick * 3, draw_as_segment=True)

        if show_forcepolygon:
            if self.force:
                vertex_index = self.form.vertex_index()
                for vertex in vertices:
                    i = vertex_index[vertex]
                    force_faces = list(self.force.faces())
                    face = self.force.face_coordinates(force_faces[i])
                    polygon = Polygon(face)
                    S = Scale.from_factors(3 * [self.scale_forcediagram])
                    T = Translation.from_vector(self.translation_forcediagram)
                    polygon = polygon.transformed(S)
                    polygon = polygon.transformed(T)
                    self.app.add(polygon, linewidth=base_thick * 3, facecolor=(0.8, 0.8, 0.8))

    def draw_form_sym(self, print_sym=False):
        """Draw the form diagram symmetry edges with the respective colors.

        Parameters
        ----------
        print_sym : bool, optional
            If text should be added to the symmetry edges, by default True

        Returns
        -------
        None
            The plotter is updated in place.
        """

        texts = {}
        colors = {}
        edges = []

        base_thick = self.settings['size.edge.base_thickness'] * 2

        i_sym_max = 0
        for edge in self.form.edges_where({'_is_edge': True}):
            edges.append(edge)
            i_sym = self.form.edge_attribute(edge, 'sym_key')
            if i_sym is None:
                raise NameError('Check if symmetry is applied to to the problem formulation.')
            if i_sym > i_sym_max:
                i_sym_max = i_sym

        colormap = plt.cm.get_cmap('hsv')  # gist_ncar nipy_spectral, Set1, Paired coolwarm
        linspace = [1.0/(i_sym_max) * i for i in range(i_sym_max + 1)]
        colormaps = [rgb_to_hex(colormap(i)[:3]) for i in linspace]

        for edge in edges:
            i_sym = self.form.edge_attribute(edge, 'sym_key')
            colors[edge] = colormaps[i_sym]
            if print_sym:
                texts[edge] = str(i_sym)

        self.formartist = self.app.add(
            self.form,
            edges=edges,
            edgewidth=base_thick,
            # edge_text=texts,  # not working for COMPAS 1.15
            edgecolor=colors,
            show_vertices=False,
            show_faces=False
        )
        if print_sym:
            self.formartist.draw_edgelabels(text=texts)

    def draw_force(self, scale=None, show_edges=True, show_vertices=False, show_faces=False, show_independents=False):
        """Draw the force diagram associated with the form diagram.

        Parameters
        ----------
        scale : float, optional
            If a scale should be defined. If None, a convenient scale is compuuted, by default None
        show_edges : bool, optional
            Whether or not edges are shown in the force diagram, by default True
        show_vertices : bool, optional
            Whether or not vertices are shown in the force diagram, by default False
        show_faces : bool, optional
            Whether or not faces are shown in the force diagram, by default False
        show_independentrs : bool, optional
            Whether or not independent edges are highlighted, by default False

        Returns
        -------
        None
            The plotter is updated in place.
        """
        if not self.force:
            from compas_tno.algorithms import reciprocal_from_form
            self.force = reciprocal_from_form(self.form)

        force = self.force

        base_width = self.settings['size.edge.base_thickness']
        width = {}

        edgecolor = {}
        for edge in force.edges():
            edgecolor[edge] = Color.black()
            width[edge] = base_width
        if show_independents:
            for edge in force.edges_where_dual({'is_ind': True}):
                edgecolor[edge] = self.settings['color.edges.independent']
                width[edge] = base_width * 3

        force_bbox = self.force.bounding_box_xy()
        xmin_force, xmax_force = force_bbox[0][0], force_bbox[1][0]
        ymin_force, ymax_force = force_bbox[0][1], force_bbox[2][1]
        dx_force, dy_force = (xmax_force - xmin_force), (ymax_force - ymin_force)

        if self.form:

            form_bbox = self.form.bounding_box_xy()
            xmin_form, xmax_form = form_bbox[0][0], form_bbox[1][0]
            ymin_form, ymax_form = form_bbox[0][1], form_bbox[2][1]
            dx_form, dy_form = (xmax_form - xmin_form), (ymax_form - ymin_form)

            scale = scale or min(dy_form/dy_force, dx_form/dx_force)
            print('Diagram scale factor:', scale)

            tx = dx_form * 1.2 + (xmax_form - xmin_force) * scale
            ty = (ymin_form - ymin_force) * scale

            frame = Frame([xmin_form, ymin_form, 0.0], (1, 0, 0), (0, 1, 0))

            translation = [tx, ty, 0]

            self.scale_forcediagram = scale
            self.translation_forcediagram = translation

            S = Scale.from_factors(3 * [scale], frame=frame)
            T = Translation.from_vector(translation)

            force = force.transformed(S)
            force = force.transformed(T)

        self.app.add(force,
                     show_edges=show_edges,
                     show_vertices=show_vertices,
                     show_faces=show_faces,
                     edgecolor=edgecolor,
                     edgewidth=width,
                     )

    def draw_edgelabels(self, text=None):
        """Draw Labels to the vertices

        Parameters
        ----------
        text : dict, optional
            The dictionary to print, if None, the keys are displayed, by default None
        """

        if not self.formartist:
            return

        self.formartist.draw_edgelabels(text=text)

    def draw_arch_lines(self, H=1.00, L=2.0, x0=0.0, thk=0.20, total_nodes=50, stereotomy=False, close_bottom=True):
        """Helper to draw the lines of intrados and extrados of an arch for given parameters

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
        total_nodes : int, optional
            Density of the shape, equals to the number of blocks, by default 50
        stereotomy : bool, optional
            Whether or not interfaces of the stereotomy should be drawn, by default False
        close_bottom : bool, optional
            Whether or not the last interfaces of the arch should be drawn, by default True

        Returns
        -------
        None
            Lines are added to the plotter
        """

        lines_intrados = []
        lines_extrados = []

        radius = H / 2 + (L**2 / (8 * H))
        ri = radius - thk/2
        re = radius + thk/2
        spr = math.atan2((L/2), (radius - H))
        tot_angle = 2*spr
        angle_init = (math.pi - tot_angle)/2
        an = tot_angle / (total_nodes)
        zc = radius - H
        xc = L/2 + x0
        i = 0

        for i in range(total_nodes):
            angle_i = angle_init + i * an
            angle_f = angle_init + (i + 1) * an

            xii = xc - ri * math.cos(angle_i)
            xif = xc - ri * math.cos(angle_f)
            zii = ri * math.sin(angle_i) - zc
            zif = ri * math.sin(angle_f) - zc

            xei = xc - re * math.cos(angle_i)
            xef = xc - re * math.cos(angle_f)
            zei = re * math.sin(angle_i) - zc
            zef = re * math.sin(angle_f) - zc

            line_i = Line([xii, zii, 0.0], [xif, zif, 0.0])
            lines_intrados.append(line_i)

            line_e = Line([xei, zei, 0.0], [xef, zef, 0.0])
            lines_extrados.append(line_e)

        interfaces = []

        if close_bottom:
            interfaces = [0, total_nodes]
        if stereotomy:
            interfaces = list(range(total_nodes + 1))

        for i in interfaces:
            angle = angle_init + i * an

            xi = xc - ri * math.cos(angle)
            xf = xc - re * math.cos(angle)
            zi = ri * math.sin(angle) - zc
            zf = re * math.sin(angle) - zc

            line = Line([xi, zi, 0.0], [xf, zf, 0.0])
            lines_intrados.append(line)

        # update this when collections are available

        for le in lines_extrados:
            self.app.add(le, draw_as_segment=True, opacity=self.settings['shape.opacity'], color=self.settings['color.edges.shape'])
        for li in lines_intrados:
            self.app.add(li, draw_as_segment=True, opacity=self.settings['shape.opacity'], color=self.settings['color.edges.shape'])

        return

    def make_gif(self, images_file, filepath, interval=10):
        """Make a gif from a list of image files

        Parameters
        ----------
        images_file : str
            List with the files to make the image. Note that list should be ordered
        filepath : str
            Pathh to save the .gif file
        interval : int, optional
            Interval between images, by default 10

        Returns
        -------
        None
            GIF is saved in the address provided
        """

        from PIL import Image

        images = []
        for path in images_file:
            images.append(Image.open(path))

        images[0].save(filepath, save_all=True, append_images=images[1:], optimize=False, duration=interval * 10, loop=0)
