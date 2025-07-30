from compas_viewer.config import Config
from compas_viewer.viewer import Viewer
from compas.geometry import Cylinder, Point
from compas.colors import Color
import math


class TNOViewer:
    """
    Wrapper for a COMPAS Viewer to visualize:
    - A form-diagram as grouped pipes scaled by force ('Form').
    - Optional shape meshes for intrados and extrados in respective groups ('Intrados', 'Extrados').
    - Cracks at vertices touching intrados/extrados ('Cracks').
    """
    def __init__(self, form, shape=None, config=None):
        """
        Parameters
        ----------
        form : FormDiagram
            A COMPAS form-diagram with edge attribute 'q' and vertex attributes 'lb', 'ub'.
        shape : object, optional
            An object with attributes `intrados` and `extrados`, each a COMPAS mesh.
        config : Config, optional
            A compas_viewer Config object. If None, a default is used.
        """
        if config is None:
            config = Config()
        self.viewer = Viewer(config=config)
        self.form = form
        self.shape = shape
        self.groups = {}
        self.settings = {
            'form_max_thk': 0.1,  # Maximum pipe radius for the largest detected force
            'intrados_color': Color.blue().darkened(30),
            'extrados_color': Color.green().darkened(30),
            'form_color': Color.red(),
            'shape_color': Color.grey().lightened(30),
            'shape_opacity': 0.4,
            'thrust_opacity': 0.9,
            'crack_opacity': 0.9,
            'crack_size': 20,  # Size of crack points
            'cracks_tol': 1e-3,  # Tolerance for comparing vertex z to lb/ub values
        }

    def clear(self):
        """Clear all objects and reset all groups."""
        self.viewer.scene.clear()
        self.groups.clear()

    def add_shape(self):
        """Create and add groups for intrados and extrados meshes."""
        if not self.shape:
            return
        if getattr(self.shape, 'intrados', None):
            grp = self.viewer.scene.add_group(name='Intrados', show=True)
            grp.add(
                self.shape.intrados,
                show_faces=True,
                show_lines=False,
                facecolor=self.settings['shape_color'],
                opacity=self.settings['shape_opacity'],
                name='intrados'
            )
            self.groups['intrados'] = grp
        if getattr(self.shape, 'extrados', None):
            grp = self.viewer.scene.add_group(name='Extrados', show=True)
            grp.add(
                self.shape.extrados,
                show_faces=True,
                show_lines=False,
                facecolor=self.settings['shape_color'],
                opacity=self.settings['shape_opacity'],
                name='extrados'
            )
            self.groups['extrados'] = grp

    def add_form(self):
        """
        Create and add cylindrical pipes for each form-diagram edge into 'Form' group.

        Parameters
        ----------
        max_thick : float, optional
            Maximum pipe radius for the largest detected force.
        """
        max_thick = self.settings['form_max_thk']
        grp = self.viewer.scene.add_group(name='Form', show=True)

        edges = list(self.form.edges_where({'_is_edge': True}))
        forces = [self.form.edge_attribute(e, 'q') * self.form.edge_length(e) for e in edges]
        f_max = math.sqrt(max(abs(max(forces)), abs(min(forces)))) or 1e-6
        for edge in edges:
            q = self.form.edge_attribute(edge, 'q')
            line = self.form.edge_line(edge)
            length = line.length
            force = math.sqrt(abs(q * length))
            if force < 1e-3:
                continue
            radius = (force / f_max) * max_thick
            cyl = Cylinder.from_line_and_radius(line, radius)
            grp.add(cyl, 
                    name=f"thrust_{edge}", 
                    color=self.settings['form_color'], 
                    opacity=self.settings['thrust_opacity'],)
        self.groups['form'] = grp

    def add_cracks(self):
        """
        Identify vertices where the form touches intrados/extrados and add them to 'Cracks' group.

        Parameters
        ----------
        tol : float, optional
            Tolerance for comparing vertex z to lb/ub values.
        """
        tol = self.settings['cracks_tol']
        grp = self.viewer.scene.add_group(name='Cracks', show=True)

        for key in self.form.vertices():
            x, y, z = self.form.vertex_coordinates(key)
            lb = self.form.vertex_attribute(key, 'lb')
            ub = self.form.vertex_attribute(key, 'ub')
            if lb is not None and abs(z - lb) < tol:
                grp.add(Point(x, y, z), 
                        name=f"intrados_crack_{key}", 
                        pointsize=self.settings['crack_size'], 
                        pointcolor=self.settings['intrados_color'],
                        opacity=self.settings['crack_opacity'])
            if ub is not None and abs(ub - z) < tol:
                grp.add(Point(x, y, z), 
                        name=f"extrados_crack_{key}", 
                        pointsize=self.settings['crack_size'], 
                        pointcolor=self.settings['extrados_color'],
                        opacity=self.settings['crack_opacity'])
        self.groups['cracks'] = grp

    def show(self):
        """
        Render intrados/extrados (if provided), form-diagram pipes, and cracks by default.
        This method clears the viewer, adds the shape, form, and cracks, and then displays the viewer.

        """
        self.clear()
        self.add_shape()
        self.add_form()
        self.add_cracks()
        self.viewer.show()
