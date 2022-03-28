import math

from compas.data import Data

from compas_tno.problems import set_up_convex_optimisation
from compas_tno.problems import set_up_general_optimisation

from compas_tno.solvers import run_optimisation_scipy
from compas_tno.solvers import run_optimisation_MATLAB
from compas_tno.solvers import run_optimisation_CVXPY
from compas_tno.solvers import run_optimisation_MMA
from compas_tno.solvers import run_optimisation_ipopt

from compas_tno.utilities import apply_selfweight_from_shape
from compas_tno.utilities import apply_selfweight_from_pattern
from compas_tno.utilities import get_shape_ub
from compas_tno.utilities import get_shape_lb

from compas_tno.utilities import apply_horizontal_multiplier
from compas_tno.utilities import apply_envelope_from_shape
from compas_tno.utilities import apply_bounds_on_q
from compas_tno.utilities import apply_bounds_reactions
from compas_tno.utilities import apply_envelope_on_xy

from compas_tno.diagrams import FormDiagram
from compas_tno.shapes import Shape
from compas_tno.optimisers import Optimiser

from numpy import array
from scipy import interpolate


class Analysis(Data):
    """The ``Analysis`` class unites the ``FormDiagram``, the ``Shape`` and the ``Optimiser`` to compute the optimisarion.

    Its metohds include assignment of loads, partial supports, cracks, etc.

    Parameters
    ----------
    name: str, optional
        The name of the analysis.
        Defaults to "Analysis".

    Attributes
    ----------
    node : dict
        The node dictionary. Each key in the node dictionary
        represents a node of the network and maps to a dictionary of
        node attributes.
    edge : dict of dict
        The edge dictionary. Each key in the edge dictionary
        corresponds to a key in the node dictionary, and maps to a dictionary
        with connected nodes. In the latter, the keys are again references
        to items in the node dictionary, and the values are dictionaries
        of edge attributes. For example, an edge between node 1 and node 2 is represented as follows
        ``Graph.edge[1][2] -> {...}``
    adjacency : dict of dict
        The edges of the graph are directed.
        The undirected connectivity information is represented in the adjacency dict.
    attributes : dict
        A dictionary of miscellaneous information about the graph.
    default_node_attributes : dict
        A dictionary mapping node attribute names to their default values.
    default_edge_attributes : dict
        A dictionary mapping edge attribute names to their default values.
    data : dict
        A dictionary representing the essential data of a graph that can be used in serialization
        processes.

    Examples
    --------
    >>>

    """

    def __init__(self, name="Analysis"):
        self.settings = {}
        self.form = None
        self.shape = None
        self.optimiser = None
        self.name = name

    @property
    def data(self):
        """dict : A data dict representing the shape data structure for serialization.
        """
        data = {
            'form': self.form.data,
            'optimiser': self.optimiser.data,
            'shape': self.shape.data,
            'name': self.name,
            'settings': self.settings
        }
        return data

    @data.setter
    def data(self, data):
        if 'data' in data:
            data = data['data']
        self.settings = data.get('settings') or {}
        self.name = data.get('name')

        formdata = data.get('form', None)
        shapedata = data.get('shape', None)
        optimiserdata = data.get('optimiser', None)

        self.form = None
        self.shape = None
        self.optimiser = None

        if formdata:
            self.form = FormDiagram.from_data(formdata)
        if shapedata:
            self.shape = Shape.from_data(shapedata)
        if optimiserdata:
            self.optimiser = Optimiser.from_data(optimiserdata)

    @classmethod
    def from_elements(cls, shape, form, optimiser):
        """Create a analysis from the elements of the problem (form, shape and optimiser) """

        analysis = cls()
        analysis.shape = shape
        analysis.form = form
        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def from_form_and_optimiser(cls, form, optimiser):
        """Create a analysis from the elements of the problem (form and optimiser)"""

        analysis = cls()
        analysis.shape = None
        analysis.form = form
        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def from_form_and_shape(cls, form, shape):
        """Create a analysis from the elements of the problem (form and shape) """

        analysis = cls()
        analysis.shape = shape
        analysis.form = form
        analysis.optimiser = None

        return analysis

    def apply_selfweight(self):
        """Invoke method to apply selfweight to the nodes of the form diagram based on the shape"""

        norm = self.optimiser.settings.get('normalize_loads', True)
        apply_selfweight_from_shape(self.form, self.shape, normalize=norm)

        return

    def apply_selfweight_from_pattern(self, pattern, plot=False):
        """Apply selfweight to the nodes considering a different Form Diagram to locate loads.

            Warning, the base pattern has to coincide with nodes from the original form diagram.
        """

        apply_selfweight_from_pattern(self.form, pattern, plot=plot)

        return

    def apply_hor_multiplier(self, multiplier, component):
        """Apply a multiplier on the selfweight to the nodes of the form diagram based"""

        apply_horizontal_multiplier(self.form, lambd=multiplier, direction=component)

        return

    def apply_envelope(self):
        """Invoke method to apply ub and lb to the nodes based on the shape's intrados and extrados"""

        apply_envelope_from_shape(self.form, self.shape)

        return

    def apply_bounds_on_q(self, qmax=0.0, qmin=-10000.0):
        """Invoke method to apply bounds on the force densities of the pattern (qmax, qmin)"""

        apply_bounds_on_q(self.form, qmax=qmax, qmin=qmin)

        return

    def apply_envelope_with_damage(self):
        """Apply ub and lb to the nodes based on the shape's intrados and extrados and in the intra/extra damaged"""

        form = self.form
        shape = self.shape
        extrados_damage = array(shape.extrados_damage.vertices_attributes('xyz'))
        intrados_damage = array(shape.intrados_damage.vertices_attributes('xyz'))

        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            ub_ = get_shape_ub(shape, x, y)
            lb_ = get_shape_lb(shape, x, y)
            lb_damage = float(interpolate.griddata(intrados_damage[:, :2], intrados_damage[:, 2], [x, y]))
            ub_damage = float(interpolate.griddata(extrados_damage[:, :2], extrados_damage[:, 2], [x, y]))
            if math.isnan(lb_damage) or math.isnan(lb_):
                lb_ = -1 * shape.datashape['t']
                # print('Shape interpolation got NaN, check results (x,y): ({0:.3f}, {1:.3f})'.format(x, y))
            else:
                lb_ = max(lb_, lb_damage)
            ub_ = min(ub_, ub_damage)
            form.vertex_attribute(key, 'ub', value=ub_)
            form.vertex_attribute(key, 'lb', value=lb_)

        self.form = form  # With correct forces

        return

    def apply_target(self):
        """Apply target to the nodes based on the shape's target surface"""

        form = self.form
        shape = self.shape

        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            form.vertex_attribute(key, 'target', value=shape.get_middle(x, y))

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form  # With correct forces

        return

    def apply_envelope_on_xy(self, c=0.5):
        """_summary_

        Parameters
        ----------
        c : float, optional
            Distance in (x, y) to constraint the nodes limiting the hor. movement, by default 0.5

        """

        apply_envelope_on_xy(self.form, c=c)

        return

    def apply_cracks(self, key, position):
        """Apply cracks on the nodes (key) and in the positions up / down"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form  # With correct forces

        return

    def apply_reaction_bounds(self, assume_shape=None):
        """Apply limit thk to be respected by the anchor points"""

        apply_bounds_reactions(self.form, self.shape, assume_shape)

        return

    def set_up_optimiser(self):
        """With the data from the elements of the problem compute the matrices for the optimisation"""

        objective = self.optimiser.settings['objective']
        features = self.optimiser.settings['features']

        if objective == 'loadpath' and 'fixed' in features:
            self.optimiser.settings['type'] = 'convex'
            set_up_convex_optimisation(self)
        else:
            self.optimiser.settings['type'] = 'nonlinear'
            self = set_up_general_optimisation(self)

        return

    def run(self):
        """With the data from the elements of the problem compute the matrices for the optimisation"""

        opt_type = self.optimiser.settings['type']
        solver = self.optimiser.settings['solver']
        library = self.optimiser.settings['library']

        if opt_type == 'convex':
            if solver == 'MATLAB':
                run_optimisation_MATLAB(self)
            elif solver == 'CVXPY':
                run_optimisation_CVXPY(self)
            else:
                raise NotImplementedError('Only <CVXPY> and <MATLAB> are suitable for this optimisation')
        elif library == 'pyOpt':
            self = run_optimisation_MMA(self)
        elif solver == 'MMA':
            self = run_optimisation_MMA(self)
        elif solver == 'IPOPT':
            self = run_optimisation_ipopt(self)
        else:
            self = run_optimisation_scipy(self)

        return
