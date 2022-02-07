from compas_plotters import Plotter
from compas.geometry import Vector
from compas.geometry import Point
from compas.geometry import Rotation
from compas_tno.shapes import Shape
import math


__all__ = ['TNOPlotter']


class TNOPlotter(object):
    """A Class plot 2D forms and shapes.

    Parameters
    ----------
    form : FormDiagram, optional
        The FormDiagram you want to plot, by default None
    form_base : FormDiagram, optional
        The base FormDiagram if the pattern can move, by default None
    shape : Shape, optional
        The Shape of masonry to plot, by default None

    Attributes
    ----------

    None

    """

    def __init__(self, form=None, form_base=None, shape=None, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.title = 'Plotter'
        self.app = None
        self._form = form
        self._shape = shape
        self._form_base = form_base
        self._formartist = None
        self._otherartists = []
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

            'camera.target': [0, 0, 0],
            'camera.distance': 40,
            'camera.rx': -45,
            'camera.rz': 45,
            'camera.fov': 40,
            'camera.show.grid': True,
            'camera.show.axis': True,
            'camera.figsize': (8, 8),

            'size.vertex': 5.0,
            'size.edge.max_thickness': 10.0,
            'size.edge.base_thickness': 1.0,
            'size.edge.normals': 0.5,
            'size.reaction.head_width': 0.1,
            'size.reaction.body_width': 0.01,
            'size.reactionlabel': 20,

            'scale.reactions': 0.05,
            'scale.loads': 0.001,
            'opacity.shapes': 0.5,

            'color.edges.form': (255, 0, 0),
            'color.edges.reactions': (170, 170, 170),
            'color.edges.shape': (200, 200, 200),
            'color.edges.form_base': (200, 200, 200),
            'color.edges.independent': (255, 0, 120),
            'color.vertex.supports': (170, 170, 170),
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
    def formartist(self):
        """The formartist property."""
        return self._formartist

    @formartist.setter
    def formartist(self, value):
        if not self._formartist:
            self._formartist = value

    def initiate_app(self):
        """ Initiate the Plotter with the default camera options.

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app = Plotter(figsize=self.settings['camera.figsize'])

    def show_solution(self):
        """ Show the thrust network, with the shape according to the settings.

        Returns
        -------
        None
            The objects are updated in place
        """

        self.draw_form()
        self.draw_supports()
        self.draw_cracks()
        self.zoom_extends()
        self.show()

    def zoom_extends(self):
        """ Wrapper to extend the objects in the active view.

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app.zoom_extents()

    def show(self):
        """ Display the plotter in the screen.

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.app.zoom_extents()
        self.app.show()

    def clear(self):
        """Clear the plotter elements.

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app.clear()

    def draw_form(self):
        """Draw the Form Diagram with or without thicknesses of edges scaled as forces in the edges.

        Returns
        -------
        None
            The plotter is updated in place
        """

        base_thick = self.settings['size.edge.base_thickness']
        max_thick = self.settings['size.edge.max_thickness']

        if self.settings['show.edge.thickness']:
            forcedensities = self.form.edges_attribute('q')
            lengths = [self.form.edge_length(u, v) for u, v in self.form.edges()]
            forces = [abs(forcedensities[i] * lengths[i]) for i in range(len(lengths))]
            fmax = max(forces)
            edgewidths = {edge: forces[i] / fmax * max_thick for i, edge in enumerate(self.form.edges())}
        else:
            edgewidths = base_thick

        edges = []
        faces = []

        if self.settings['show.edges']:
            edges = list(self.form.edges_where({'_is_edge': True}))

        if self.settings['show.faces']:
            faces = list(self.form.faces())

        formartist = self.app.add(
            self.form,
            vertices=[],
            edges=edges,
            faces=faces,
            edgewidth=edgewidths,
            edgecolor=_norm(self.settings['color.edges.form'])
        )

        self.formartist = formartist

    def draw_cracks(self):
        """Adds to the basic plot, the cracks which are the points of the mesh that touch intrados or extrados.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        tol = self.settings['tol.forces']
        color_intra = _norm(self.settings['color.vertex.intrados'])
        color_extra = _norm(self.settings['color.vertex.extrados'])
        color_outside = _norm(self.settings['color.vertex.outside'])

        cracks = []

        for key in self.form.vertices():
            x, y, z = self.form.vertex_coordinates(key)
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

        for point in cracks:
            x, y, z = point['xyz']
            color = point['color']
            pt = Point(x, y, z)
            pointartist = self.app.add(pt, facecolor=color, size=self.settings['size.vertex'])
            self._otherartists.append(pointartist)

    def draw_supports(self):
        """Add the supports as points in the mesh.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if self.settings['show.supports']:
            supportcolor = _norm(self.settings['color.vertex.supports'])

            for key in self.form.vertices_where({'is_fixed': True}):
                x, y, z = self.form.vertex_coordinates(key)
                pt = Point(x, y, z)
                pointartist = self.app.add(pt, facecolor=supportcolor, size=self.settings['size.vertex'])
                self._otherartists.append(pointartist)

    def draw_reactions(self):
        """Add to the plots the vector of the reaction forces.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if self.settings['show.reactions']:
            reaction_scale = self.settings['scale.reactions']

            for key in self.form.vertices_where({'is_fixed': True}):
                x, y, z = self.form.vertex_coordinates(key)
                rx = self.form.vertex_attribute(key, '_rx') * reaction_scale
                ry = self.form.vertex_attribute(key, '_ry') * reaction_scale
                rz = self.form.vertex_attribute(key, '_rz') * reaction_scale
                if self.settings['show.reactions.emerging']:
                    r = Vector(-rx, -ry, -rz)
                    pt = Point(x, y, z)
                else:
                    r = Vector(rx, ry, rz)
                    pt = Point(x - rx, y - ry, z - rz)
                vectorartist = self.app.add(r, point=pt, color=_norm(self.settings['color.edges.reactions']))
                self._otherartists.append(vectorartist)

    def draw_vectors(self, vectors=[], bases=[]):
        """Helper to add vectors to the plotter.

        Parameters
        ----------
        vectors : list, optional
            The list of vectors to plot, by default []
        bases : list, optional
            The list with the location of the base of the vectors, by default []

        Returns
        -------
        None
            Plotter app is updated.

        """

        if not(len(vectors) == len(bases)):
            raise ValueError('Please provide the vector and the bases')

        for i in range(len(vectors)):
            vector = vectors[i]
            point = bases[i]
            vectorartist = self.app.add(vector, point=point)
            self._otherartists.append(vectorartist)

    def draw_mesh(self, mesh=None, show_edges=True, show_vertices=False, show_faces=False):
        """ Initiate the Plotter with the default camera options.

        Returns
        -------
        None
            The objects are updated in place
        """

        if not mesh:
            mesh = self.form
        self.app.add(mesh,
                     show_edges = show_edges,
                     show_vertices = show_vertices,
                     show_faces = show_faces
                     )

    def draw_form_xz(self):
        """Plot the form diagram rotated 90 degrees.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        axis = Vector(1.0, 0, 0)
        self.form = self.form.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))

        self.draw_form()

    def draw_shape(self):
        """Adds the shape to the plot.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if not self.shape:
            self.shape = Shape.from_formdiagram_and_attributes(self.form)

        intrados = self.shape.intrados
        extrados = self.shape.extrados

        self.app.add(
            intrados,
            opacity=0.5,
            edgecolor=_norm(self.settings['color.edges.shape']),
            show_vertices=False
        )

        self.app.add(
            extrados,
            opacity=0.5,
            edgecolor=_norm(self.settings['color.edges.shape']),
            show_vertices=False
        )

    def draw_shape_xz(self):
        """Plot the shape rotated 90 degrees.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if not self.shape:
            self.shape = Shape.from_formdiagram_and_attributes(self.form)

        axis = Vector(1.0, 0, 0)
        self.shape.intrados = self.shape.intrados.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))
        self.shape.extrados = self.shape.extrados.transformed(Rotation.from_axis_and_angle(axis, -math.pi/2))

        self.draw_shape()

    def draw_base_form(self):
        """Adds to the plot the base mesh which is the mesh before the nodes moved horizontally.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        if not self.form_base:
            return

        if self.settings['show.edges']:
            edges = list(self.form_base.edges_where({'_is_edge': True}))

        self.app.add(
            self.form_base,
            edges=edges,
            opacity=0.5,
            edgecolor=_norm(self.settings['color.edges.form_base']),
            show_vertices=False,
            show_faces=False
        )

    def draw_form_independents(self, show_text=True):
        """Draw the form diagram with highlight in the independent edgess.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        edges = list(self.form.edges_where({'_is_edge': True}))
        base_thick = self.settings['size.edge.base_thickness']

        color = {}
        width = {}
        text = {}
        i = 0
        for edge in edges:
            if self.form.edge_attribute(edge, 'is_ind'):
                color[edge] = _norm(self.settings['color.edges.independent'])
                width[edge] = base_thick * 3.0
                if show_text:
                    text[edge] = str(i)
                i = i+1
            else:
                color[edge] = (0, 0, 0)
                width[edge] = base_thick

        self.app.add(
                    self.form,
                    edges=edges,
                    edgewidth=width,
                    edge_text=text,
                    edgecolor=color,
                    show_vertices=False,
                    show_faces=False
                    )

    def draw_form_sym(self):
        """Draw the form diagram symmetry edges with the respective colors.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        raise NotImplementedError

    def draw_force(self):
        """Draw the force diagram associated with the form diagram.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        raise NotImplementedError


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
