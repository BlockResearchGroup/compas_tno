import copy
import math

from compas_tno.shapes import Shape

from compas_plotters import MeshPlotter

from compas_tno.algorithms.problems import initialise_form
from compas_tno.algorithms.problems import initialise_problem
from compas_tno.algorithms import set_up_nonlinear_optimisation
from compas_tno.algorithms import set_up_convex_optimisation

from compas_tno.algorithms import run_optimisation_scipy
from compas_tno.algorithms import run_optimisation_MATLAB
from compas_tno.algorithms import run_optimisation_MMA
from compas_tno.algorithms import run_optimisation_ipopt

__all__ = ['Analysis']

class Analysis(object):

    """The ``Analysis`` class puts together the FormDiagram, the Shape and the Optimiser, assign loads, partial supports, cracks and much more to come...

    Notes
    -----
    A ``Analysis`` has the following constructor functions

    *   ``initialise`` : Construct the Analysis for a given FormDiagram, Shape and Optimiser.

    A ``Analysis`` is composed by the following elements:

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
        *  ``results_summary``: Work in progress

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

    @classmethod
    def from_form_and_optimier(cls, form, optimiser):
        """Create a analysis from the elements of the problem """

        analysis = cls()
        analysis.shape = None
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
            form_.vertex_attribute(key, 'z', value=z)
            form.vertex_attribute(key, 'target', value=z)

        pzt = 0
        for key in form.vertices():
            pz = form_.vertex_area(key)
            form.vertex_attribute(key, 'pz', value = pz)
            pzt += pz

        if shape.data['type'] == 'arch':
            pzt = 0
            for key in form.vertices():
                form.vertex_attribute(key, 'pz', value = 1.0)
                if form.vertex_attribute(key, 'is_fixed') == True:
                    form.vertex_attribute(key, 'pz', value = 0.5)
                pzt += form.vertex_attribute(key, 'pz')

        factor = total_selfweight/pzt

        for key in form.vertices():
            pzi = factor * form.vertex_attribute(key, 'pz')
            form.vertex_attribute(key, 'pz', value = pzi)

        self.form = form

        return


    def apply_fill_load(self, plot=False):

        shape = self.shape
        form = self.form
        vol_calc = 0.
        vertices_under_fill = []
        if shape.fill:
            for key in form.vertices():
                x, y, z = form.vertex_coordinates(key)
                z_fill = shape.get_ub(x, y)
                z_brick = shape.get_ub_bricks(x, y)
                if z_fill > z_brick:
                    dz = z_fill - z_brick
                    proj_area = form.vertex_projected_area(key)
                    form.vertex_attribute(key, 'fill_volume', proj_area*dz)
                    vol_calc += proj_area*dz
                    vertices_under_fill.append(key)
            fill_load = shape.compute_fill_weight()
            fill_vol = shape.fill_volume
            print('Volume of Fill in this Form-Diagram:', vol_calc, 'against a shape fill vol and load:', fill_vol, fill_load)
            for key in vertices_under_fill:
                load_fraction = form.vertex_attribute(key, 'fill_volume')/vol_calc * fill_load
                pz0 = form.vertex_attribute(key, 'pz')
                form.vertex_attribute(key, 'pz', pz0 + load_fraction)
                form.vertex_attribute(key, 'fill_load', load_fraction)
            if plot:
                plotter = MeshPlotter(form, figsize=(10, 10))
                plotter.draw_edges()
                plotter.draw_vertices(keys=vertices_under_fill, text={key: str(round(form.vertex_attribute(key, 'fill_load'), 2)) for key in vertices_under_fill})
                plotter.show()
        else:
            print('Warning, no Fill was assigned to Shape')


        return

    def apply_pointed_load(self, keys, magnitudes, proportional=True, component='pz', component_get='pz'):
        """Apply pointed load a node of the form diagram based on the shape"""

        if isinstance(keys, list) == False:
            keys = [keys]
        if isinstance(magnitudes, list) == False:
            magnitudes = [magnitudes]
        if len(keys) == len(magnitudes):
            pass
        elif len(magnitudes) is 1:
            magnitudes = magnitudes * len(keys)
        else:
            print('Error, check the number of nodes to apply loads and magnitude of forces!')
        print('magnitudes',magnitudes)
        form = self.form
        i = 0
        for key in keys:
            p0 = form.vertex_attribute(key, component_get)
            if proportional:
                pnew = (p0 + magnitudes[i] * p0)
            else:
                pnew = magnitudes[i]
            form.vertex_attribute(key, component, pnew)
            print('Load:', pnew, 'applied to key:', key, 'in direction:', component)
            i += 1

        self.form = form # With correct forces

        return

    def apply_hor_multiplier(self, multiplier, component):
        """Apply a multiplier on the selfweight to the nodes of the form diagram based"""

        form = self.form
        shape = self.shape

        pzt = 0
        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            pz0 = form.vertex_attribute(key, 'pz')
            pzi = multiplier * pz0
            form.vertex_attribute(key, component, pzi)
            pzt += pzi

        print('Load applied in ', component, 'total: ', pzt)

        self.form = form

        return


    def apply_envelope(self):
        """Apply ub and lb to the nodes based on the shape's intrados and extrados"""


        form = self.form
        shape = self.shape

        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            ub_ = shape.get_ub(x, y)
            lb_ = shape.get_lb(x, y)
            if math.isnan(lb_):
                lb_ = -1 * shape.data['t']
                print('Shape interpolation got NaN, check results')
            form.vertex_attribute(key, 'ub', value = ub_)
            form.vertex_attribute(key, 'lb', value = lb_)

        self.form = form # With correct forces

        return

    def apply_target(self):
        """Apply target to the nodes based on the shape's target surface"""


        form = self.form
        shape = self.shape

        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            form.vertex_attribute(key, 'target', value = shape.get_middle(x, y))

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form # With correct forces

        return

    def apply_cracks(self, key, position):
        """Apply cracks on the nodes (key) and in the positions up / down"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form  # With correct forces

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
                form.vertex_attribute(key, 'b', [x_, y_])

        if shape.data['type'] == 'arch':
            H = shape.data['H']
            L = shape.data['L']
            thk = shape.data['thk']
            radius = H / 2 + (L**2 / (8 * H))
            zc = radius - H
            re = radius + thk/2
            x = math.sqrt(re**2 - zc**2)
            for key in form.vertices_where({'is_fixed': True}):
                form.vertex_attribute(key, 'b', [x - L/2, 0.0])

        if shape.data['type'] == 'dome_spr':
            # print(shape.data['center'])
            x0 = shape.data['center'][0]
            y0 = shape.data['center'][1]
            radius = shape.data['radius']
            [_, theta_f] = shape.data['theta']
            r_proj_e = (radius + thk/2) * math.sin(theta_f)
            r_proj_m = (radius) * math.sin(theta_f)
            delt = r_proj_e - r_proj_m
            for key in form.vertices_where({'is_fixed': True}):
                x, y, _ = form.vertex_coordinates(key)
                theta = math.atan2((y - y0), (x - x0))
                x_ = delt*math.cos(theta)
                y_ = delt*math.sin(theta)
                form.vertex_attribute(key, 'b', [x_, y_])
                # print(x,y,x_,y_)

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

        if self.optimiser.data['library'] == 'MATLAB':
            self = set_up_convex_optimisation(self)
        else:
            self = set_up_nonlinear_optimisation(self)

        return

    def run(self):
        """With the data from the elements of the problem compute the matrices for the optimisation"""

        if self.optimiser.data['library'] == 'MATLAB':
            self = run_optimisation_MATLAB(self)
        elif self.optimiser.data['library'] == 'MMA':
            self = run_optimisation_MMA(self)
        elif self.optimiser.data['library'] == 'IPOPT':
            self = run_optimisation_ipopt(self)
        else:
            self = run_optimisation_scipy(self)

        return

    def limit_analysis_GSF(self, thk, thk_reduction, R, fill=None, rollers_ratio=None):

        solutions_min = []  # empty lists to keep track of the solutions for min thrust
        solutions_max = []  # empty lists to keep track of the solutions for max thrust
        size_parameters = []  # empty lists to keep track of  the parameters
        thicknesses = []
        t0 = thk

        data_diagram = self.form.parameters
        data_shape = self.shape.data

        print('Limit Analysis - GSF: For ', data_shape['type'], 'with diagram ', data_diagram['type'])

        for self.optimiser.data['objective'] in ['min', 'max']:

            exitflag = 0
            thk = t0
            while exitflag == 0:

                exitflag = 0
                t_over_R = thk/R
                print('\n----- Starting the' , data_diagram['type'], 'problem for thk/R:', t_over_R, '\n')
                data_shape['thk'] = thk
                self.shape = Shape.from_library(data_shape)
                if fill:
                    self.shape.add_fill_with_height(fill)
                swt = self.shape.compute_selfweight()
                if rollers_ratio:
                    self.form.set_boundary_rollers(total_rx=[rollers_ratio[0]*swt]*2, total_ry=[rollers_ratio[1]*swt]*2)

                self.apply_selfweight()
                self.apply_envelope()
                if fill:
                    self.apply_fill_load()
                self.set_up_optimiser()
                self.run()

                exitflag = self.optimiser.exitflag  # get info if optimisation was succeded ot not
                fopt = self.optimiser.fopt  # objective function optimum value
                fopt_over_weight = fopt/swt  # divide by selfweight

                if exitflag == 0:
                    if self.optimiser.data['objective'] == 'min':
                        solutions_min.append(fopt_over_weight)
                        size_parameters.append(t_over_R)
                        thicknesses.append(thk)
                    else:
                        solutions_max.append(fopt_over_weight)
                    thk = round(thk - thk_reduction, 4)
                    print('Reduce thickness to', thk, '\n')

        print('\n SUMMARY')
        print('Thicknesses calculated:')
        print(thicknesses)
        print('Thk/R calculated:')
        print(size_parameters)
        print('Solutions Found:')
        print(solutions_min)
        print('Solutions Found:')
        print(solutions_max)

        return thicknesses, size_parameters, solutions_min, solutions_max


    def limit_analysis_load_mult(self, load0, load_increase, load_direction, percentual=True):
        """"" W.I.P. """""

        return







