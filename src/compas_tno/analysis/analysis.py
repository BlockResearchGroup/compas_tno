import copy
import math

from compas_tno.algorithms.problems import initialise_form
from compas_tno.algorithms.problems import initialise_problem
from compas_tno.algorithms.general_solver import set_up_optimisation

__all__ = ['Analysis']

class Analysis(object):

    """The ``Analysis`` class puts together the FormDiagram, the Shape and the Optimiser, assign loads, partial supports, cracks and much more to come...

    Notes
    -----
    A ``Analysis`` has the following constructor functions

    *   ``initialise`` : Construct the Analysis for a given FormDiagram, Shape and Optimiser.

    A ``Shape`` has the following elements:

    *   ``FormDiagram``

        *   ``FormDiagram``  : Layout for the forces to follow

    *   ``Optimiser``

        *   ``Optimiser``    : Information about the solver and other specificities.

    *   ``Shape``

        *   ``Shape``    : Intrados, Extrados and Middle surface of the structure.


    The main methods available with the ``Analysis`` previous to runing the optimisation are:

        *   ``apply_selfweight``        : Apply deadload to the structure based on the shape.
        *   ``apply_pointed_load``      : Apply pointed load to the structure in a given key, in a direction 'px', 'py', 'pz'.
        *   ``apply_hor_multiplier``    : Apply horizontal multiplier in a direction 'px', 'py', 'pz'.
        *   ``apply_envelope``          : Apply Envelope to the nodes' attributes.
        *   ``aplpy_cracks``            : Create aditional constraints on the cracks.
        *   ``apply_reaction_bounds``   : Create bounds on the position of the reactions.
        *   ``apply_partial_reactions`` : Allow for partial reactions to form on the open edges.


    To run the optimisation and see results one can use:

        *   ``run``    : Method that will run the specified optimiser in the problem
        *  ``results_summary``:

    """

    __module__ = 'compas_tna.shapes'


    def __init__(self):
        self.data = { }
        self.shape = None
        self.form = None
        self.optimiser = None


    @classmethod
    def from_elements(cls, shape, form, optimiser):
        """Create a analysis from the elements of the problem """

        analysis = cls()
        analysis.shape = shape
        analysis.form = form
        analysis.optimiser = optimiser

        return analysis


    def apply_selfweight(self):
        """Apply selfweight to the nodes of the form diagram based on the shape"""

        form = self.form
        form_ = form.copy()
        shape = self.shape
        total_selfweight = shape.compute_selfweight()

        for key in form_.vertices():
            x, y, _ = form.vertex_coordinates(key)
            z = shape.get_middle(x, y)
            form_.set_vertex_attribute(key, 'z', value=z)
            form.set_vertex_attribute(key, 'target', value=z)

        pzt = 0
        for key in form.vertices():
            pz = form_.vertex_area(key)
            form.set_vertex_attribute(key, 'pz', value = pz)
            pzt += pz

        factor = total_selfweight/pzt

        for key in form.vertices():
            pzi = factor * form.get_vertex_attribute(key, 'pz')
            form.set_vertex_attribute(key, 'pz', value = pzi)

        self.form = form

        return


    def apply_pointed_load(self, key, direction, magnitude):
        """Apply pointed load a node of the form diagram based on the shape"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form # With correct forces

        return


    def apply_hor_multiplier(self, multiplier, direction):
        """Apply a multiplier on the selfweight to the nodes of the form diagram based"""

        form = self.form
        shape = self.shape

        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            form.set_vertex_attribute(key, 'ub', value = shape.get_ub(x, y))
            form.set_vertex_attribute(key, 'lb', value = shape.get_lb(x, y))

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form # With correct forces

        return


    def apply_envelope(self):
        """Apply ub and lb to the nodes based on the shape's intrados and extrados"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form # With correct forces

        return

    def apply_cracks(self, key, position):
        """Apply cracks on the nodes (key) and in the positions up / down"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form # With correct forces

        return

    def apply_reaction_bounds(self):  # Make one infividual...
        """Apply limit thk to be respected by the anchor points"""

        form = self.form
        shape = self.shape
        thk = shape.data['thk']

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        if shape.data['type'] == 'dome':
            [x0, y0] = shape.data['center']
            for key in form.vertices_where({'is_fixed': True}):
                x, y, _ = form.vertex_coordinates(key)
                theta = math.atan2((y - y0), (x - x0))
                x_ = thk/2*math.cos(theta)
                y_ = thk/2*math.sin(theta)
                form.set_vertex_attribute(key, 'b', [x_, y_])

        self.form = form # With correct forces

        return


    def apply_partial_reactions(self, key, direction, magnitude):
        """Apply ub and lb to the nodes based on the shape's intrados and extrados"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form # With correct forces

        return


    def set_up_optimiser(self):
        """With the data from the elements of the problem compute the matrices for the optimisation"""

        self = set_up_optimisation(self)

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        # self.form = form # With correct forces


        return

