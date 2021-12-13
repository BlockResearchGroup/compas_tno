from compas_plotters import Plotter
from compas.geometry import Vector
from compas.geometry import Point
from compas_plotters.artists import VectorArtist
from compas.geometry import Rotation
import math


__all__ = ['FormPlotter']


class FormPlotter(object):
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
        self._formbase = form_base
        self._formartist = None
        self._otherartists = []
        self.settings = {
            'show.thrust': True,
            'show.shape': True,
            'show.reactions': True,
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

            'size.vertex': 10.0,
            'size.edge.max_thickness': 10.0,
            'size.edge.base_thickness': 1.0,
            'size.edge.normals': 0.5,
            'size.reaction.head_width': 0.1,
            'size.reaction.body_width': 0.01,
            'size.reactionlabel': 20,

            'scale.reactions': 0.005,
            'scale.loads': 0.001,
            'opacity.shapes': 0.5,

            'color.edges.form': (255, 0, 0),
            'color.edges.reactions': (125, 125, 125),
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
        """ Initiate the Plotter with the default camera options

        Returns
        -------
        None
            The objects are updated in place
        """

        self.app = Plotter(figsize=self.settings['camera.figsize'])

    def show_solution(self):
        """ Show the thrust network, with the shape according to the settings

        Returns
        -------
        None
            The objects are updated in place
        """

        self.draw_form()
        self.draw_cracks()
        self.show()

    def show(self):
        """ Display the plotter in the screen

        Returns
        -------
        None
            The plotter is updated in place
        """

        self.app.zoom_extents()
        self.app.show()

    def clear(self):
        """Clear the plotter elements

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
        vertices = []
        faces = []

        if self.settings['show.edges']:
            edges = list(self.form.edges_where({'_is_edge': True}))

        if self.settings['show.faces']:
            faces = list(self.form.faces())

        if self.settings['show.supports']:
            vertices = list(self.form.fixed())

        formartist = self.app.add(
            self.form,
            vertices=vertices,
            edges=edges,
            faces=faces,
            edgewidth=edgewidths,
            edgecolor=_norm(self.settings['color.edges.form']),
            vertexcolor=_norm(self.settings['color.vertex.supports']),
            vertexsize=self.settings['size.vertex']
            )

        self.formartist = formartist

    def draw_cracks(self):
        """Adds to the basic plot, the cracks which are the points of the mesh that touch intrados or extrados

        Returns
        -------
        None
            The plotter is updated in place.
        """

        intrad = 1
        extrad = 1
        out = 1
        tol = self.settings['tol.forces']
        color_intra = _norm(self.settings['color.vertex.intrados'])
        color_extra = _norm(self.settings['color.vertex.extrados'])
        color_outside = _norm(self.settings['color.vertex.outside'])

        vertices = []
        color = {}

        for key in self.form.vertices():
            lb = self.form.vertex_attribute(key, 'lb')
            ub = self.form.vertex_attribute(key, 'ub')
            if not lb:
                continue
            if not ub:
                continue
            x, y, z = self.form.vertex_coordinates(key)
            if self.settings['show.cracks']:
                if abs(ub - z) < tol:
                    vertices.append(key)
                    color[key] = color_extra
                    extrad += 1
                    continue
                elif abs(lb - z) < tol:
                    vertices.append(key)
                    color[key] = color_intra
                    intrad += 1
                    continue
            if self.settings['show.vertex.outside']:
                if z - ub > tol:
                    vertices.append(key)
                    color[key] = color_outside
                    out += 1
                    continue
                elif lb - z > tol:
                    vertices.append(key)
                    color[key] = color_outside
                    out += 1
                    continue

        artist = self.formartist

        if not artist:
            self.app.add(
                self.form,
                show_edges=False,
                show_faces=False,
                show_vertices=True,
                vertices=vertices,
                vertexcolor=color,
                vertexsize=self.settings['size.vertex']
                )
        else:  # WARNING: this will delete the vertices added before as the base mesh - not good
            artist.vertex_size = self.settings['size.vertex']  # I think it does not make to have to set this here...
            artist.draw_vertices(
                vertices=vertices,
                color=color,
                # vertexsize=self.settings['size.vertex']  # THIS SHOULD BE PASSABLE HERE
            )

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
                r = Vector(-rx, -ry, -rz)
                pt = Point(x, y, z)
                vectorartist = VectorArtist(r, point=pt)
                vectorartist.draw()
                self._otherartists.append(vectorartist)

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

    def draw_shape_xz(self):
        """Plot the shape rotated 90 degrees.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        self.draw_form_xz()
        # then add the mesh self.shape.intrados and self.shape.extrados to the plot

    def draw_base_form(self):
        """Adds to the plot the base mesh which is the mesh before the nodes moved horizontally.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        self.draw_form()
        self.draw_cracks()


        if not self.form_base:
            return

        # add here the self.form_base as a light gray mesh behind the form...

    def draw_form_independents(self):
        """Draw the form diagram with highlight in the independent edgess.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        # add here the form diagram withh highlight in the inependent edges


    def draw_form_sym(self):
        """Draw the form diagram symmetry edges with the respective colors.

        Returns
        -------
        None
            The plotter is updated in place.
        """

        # add here the form diagram withh highlight in the inependent edges


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
