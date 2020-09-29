import copy
import math

from compas_tno.shapes import Shape
import compas_tno
import os

from compas_plotters import MeshPlotter

from copy import deepcopy

from compas_tno.problems import set_up_nonlinear_optimisation
from compas_tno.problems import set_up_convex_optimisation

from compas_tno.solvers import run_optimisation_scipy
from compas_tno.solvers import run_optimisation_MATLAB
from compas_tno.solvers import run_optimisation_MMA
from compas_tno.solvers import run_optimisation_ipopt

from compas_tno.diagrams import FormDiagram

from compas_tno.plotters import plot_form
from compas_tno.plotters import plot_independents

from numpy import array
from scipy import interpolate

import time

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
        self.data = {}
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
    def from_form_and_optimiser(cls, form, optimiser):
        """Create a analysis from the elements of the problem """

        analysis = cls()
        analysis.shape = None
        analysis.form = form
        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def from_form_and_shape(cls, form, shape):
        """Create a analysis from the elements of the problem """

        analysis = cls()
        analysis.shape = shape
        analysis.form = form
        analysis.optimiser = None

        return analysis

    def apply_selfweight(self):
        """Invoke method to apply selfweight to the nodes of the form diagram based on the shape"""

        self.form.selfweight_from_shape(self.shape)

        return

    def apply_selfweight_from_pattern(self, pattern, plot=False):
        """Apply selfweight to the nodes considering a different Form Diagram to locate loads. Warning, the base pattern has to coincide with nodes from the original form diagram"""

        form = self.form
        form_ = pattern
        shape = self.shape
        total_selfweight = shape.compute_selfweight()
        tol = 10e-10

        form.vertices_attribute('pz', 0.0)
        key_real_to_key = {}

        for key in form_.vertices():
            x, y, _ = form_.vertex_coordinates(key)
            for key_real in form.vertices():
                x_real, y_real, _ = form.vertex_coordinates(key_real)
                if x - tol < x_real < x + tol and y - tol < y_real < y + tol:
                    key_real_to_key[key_real] = key
                    break
            z = shape.get_middle(x, y)
            form_.vertex_attribute(key, 'z', value=z)
            form.vertex_attribute(key_real, 'target', value=z)

        pzt = 0
        for key in key_real_to_key:
            pz = form_.vertex_area(key_real_to_key[key])
            form.vertex_attribute(key, 'pz', value=pz)
            pzt += pz

        if plot:
            plotter = MeshPlotter(form, figsize=(10, 10))
            plotter.draw_edges()
            plotter.draw_vertices(text=key_real_to_key)
            plotter.show()

            plotter = MeshPlotter(form_, figsize=(10, 10))
            plotter.draw_edges()
            plotter.draw_vertices(text={key: key for key in form_.vertices()})
            plotter.show()

            plotter = MeshPlotter(form, figsize=(10, 10))
            plotter.draw_edges()
            plotter.draw_vertices(text={key: round(form.vertex_attribute(key, 'pz'), 1) for key in form.vertices()})
            plotter.show()

        if shape.data['type'] == 'arch':
            pzt = 0
            for key in form.vertices():
                form.vertex_attribute(key, 'pz', value=1.0)
                if form.vertex_attribute(key, 'is_fixed') == True:
                    form.vertex_attribute(key, 'pz', value=0.5)
                pzt += form.vertex_attribute(key, 'pz')

        factor = total_selfweight/pzt

        for key in form.vertices():
            pzi = factor * form.vertex_attribute(key, 'pz')
            form.vertex_attribute(key, 'pz', value=pzi)

        self.form = form

        return

    def apply_symmetry(self, center_point=[5.0, 5.0, 0.0]):
        """Apply symmetry to the pattern based on the center point of the Form Diagram"""

        self.form.apply_symmetry(center_point=center_point)

        return

    def apply_fill_load(self, plot=False):

        shape = self.shape
        form = self.form
        vol_calc = 0.
        vertices_under_fill = []
        swt = shape.compute_selfweight()
        proj_area_under = 0.
        if shape.fill:
            for key in form.vertices():
                x, y, z = form.vertex_coordinates(key)
                z_fill = shape.get_ub_fill(x, y)
                z_brick = shape.get_ub(x, y)
                if z_fill > z_brick:
                    dz = z_fill - z_brick
                    proj_area = form.vertex_projected_area(key)
                    form.vertex_attribute(key, 'fill_volume', proj_area*dz)
                    vol_calc += proj_area*dz
                    proj_area_under += proj_area
                    vertices_under_fill.append(key)
            fill_load = shape.compute_fill_weight()
            fill_vol = shape.fill_volume
            print('Volume of Fill in this Form-Diagram:', vol_calc, 'against a shape fill vol of:', fill_vol, ' difference (%):', 100*(vol_calc-fill_vol)/fill_vol)
            print('Load that will be added to the Structure:', fill_load, 'against a initial SWT:', swt, ' Addition of (%):', 100*(fill_load)/swt)
            print('Proj area under:', proj_area_under)
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
        print('magnitudes', magnitudes)
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

        self.form = form  # With correct forces

        return

    def apply_hor_multiplier(self, multiplier, component):
        """Apply a multiplier on the selfweight to the nodes of the form diagram based"""

        form = self.form

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
        """Invoke method to apply ub and lb to the nodes based on the shape's intrados and extrados"""

        self.form.envelope_from_shape(self.shape)

        return

    def apply_envelope_with_damage(self):
        """Apply ub and lb to the nodes based on the shape's intrados and extrados and in the intra/extra damaged"""

        form = self.form
        shape = self.shape
        extrados_damage = array(shape.extrados_damage.vertices_attributes('xyz'))
        intrados_damage = array(shape.intrados_damage.vertices_attributes('xyz'))

        for key in form.vertices():
            x, y, _ = form.vertex_coordinates(key)
            ub_ = shape.get_ub(x, y)
            lb_ = shape.get_lb(x, y)
            lb_damage = float(interpolate.griddata(intrados_damage[:, :2], intrados_damage[:, 2], [x, y]))
            ub_damage = float(interpolate.griddata(extrados_damage[:, :2], extrados_damage[:, 2], [x, y]))
            if math.isnan(lb_damage) or math.isnan(lb_):
                lb_ = -1 * shape.data['t']
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

    def apply_cracks(self, key, position):
        """Apply cracks on the nodes (key) and in the positions up / down"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form  # With correct forces

        return

    def apply_reaction_bounds(self):  # TODO: Move this to an appropriate place
        """Apply limit thk to be respected by the anchor points"""

        form = self.form
        shape = self.shape
        thk = shape.data['thk']

        if shape.data['type'] == 'dome' or shape.data['type'] == 'dome_polar':
            [x0, y0] = shape.data['center'][:2]
            for key in form.vertices_where({'is_fixed': True}):
                x, y, _ = form.vertex_coordinates(key)
                theta = math.atan2((y - y0), (x - x0))
                x_ = thk/2*math.cos(theta)
                y_ = thk/2*math.sin(theta)
                form.vertex_attribute(key, 'b', [x_, y_])

        b_manual = shape.data.get('b_manual', None)
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
                if b_manual:
                    form.vertex_attribute(key, 'b', [b_manual, 0.0])
                    print('Applied b manual')

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

        if shape.data['type'] == 'pavillionvault':
            x0, x1 = shape.data['xy_span'][0]
            y0, y1 = shape.data['xy_span'][1]
            for key in form.vertices_where({'is_fixed': True}):
                x, y, _ = form.vertex_coordinates(key)
                if x == x0:
                    form.vertex_attribute(key, 'b', [-thk/2, 0])
                elif x == x1:
                    form.vertex_attribute(key, 'b', [+thk/2, 0])
                if y == y0:
                    if form.vertex_attribute(key, 'b'):
                        b = form.vertex_attribute(key, 'b')
                        b[1] = -thk/2
                        form.vertex_attribute(key, 'b', b)
                    else:
                        form.vertex_attribute(key, 'b', [0, -thk/2])
                elif y == y1:
                    if form.vertex_attribute(key, 'b'):
                        b = form.vertex_attribute(key, 'b')
                        b[1] = +thk/2
                        form.vertex_attribute(key, 'b', b)
                    else:
                        form.vertex_attribute(key, 'b', [0, +thk/2])
                # print('x = {0:.1f} | y = {1:.1f} | b = {2}'.format(x, y, form.vertex_attribute(key, 'b')))
            # print('Pavillion limits')

        self.form = form  # With correct forces

        return

    def apply_partial_reactions(self, key, direction, magnitude):
        """Apply ub and lb to the nodes based on the shape's intrados and extrados"""

        form = self.form

        # Go over nodes and find node = key and apply the pointed load pz += magnitude

        self.form = form  # With correct forces

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

    def limit_analysis_GSF(self, thk, thk_reduction, span, thk_refined=None, limit_equal=0.01, fill_percentage=None, rollers_ratio=None, rollers_absolute=None, printout=True, plot=False, save_forms=None):

        solutions_min = []  # empty lists to keep track of the solutions for min thrust
        solutions_max = []  # empty lists to keep track of the solutions for max thrust
        thicknesses_min = []
        thicknesses_max = []
        t0 = thk
        form0 = self.form.copy()
        swt0 = self.shape.compute_selfweight()
        thk_reduction0 = thk_reduction
        data_diagram = self.form.parameters
        data_shape = self.shape.data
        ro = self.shape.ro
        last_min = 0
        last_max = 100
        last_thk_min = t0
        last_thk_max = t0
        objectives = ['min', 'max']
        exitflag = 0
        address = None

        print('Limit Analysis - GSF: For ({0}) with diagram ({1})'.format(data_shape['type'], data_diagram['type']))

        while (last_max - last_min) > limit_equal:

            if exitflag == 0:
                pass
            else:
                objectives = objectives[::-1]
                print('Warning: Did not Achieve precision required. Stopping Anyway.\n')
                break

            for self.optimiser.data['objective'] in objectives:

                self.form = form0
                exitflag = 0
                thk = t0
                thk_reduction = thk_reduction0
                count = 0

                print('\n----- Starting the [', self.optimiser.data['objective'], '] problem for intial thk:', thk)
                print('THK  |   Solved  |   Opt.Val |   Opt/W   |   THK red.  |   Setup time  |   Run time')

                while exitflag == 0 and thk > 0:
                    # if self.optimiser.data['objective'] == 'max':
                    #     if count < len(thicknesses_min):
                    #         thk = thicknesses_min[count]
                    #     else:
                    #         break
                    # try:
                    #     address_load = save_forms + '_' + 'min' + '_thk_' + str(100*thk) + '.json'
                    #     self.form = FormDiagram.from_json(address_load)
                    # except:
                    #     pass
                    data_shape['thk'] = thk
                    time0 = time.time()
                    self.shape = Shape.from_library(data_shape)
                    self.shape.ro = ro
                    swt = self.shape.compute_selfweight()
                    if fill_percentage:
                        self.shape.add_fill_with_height(max([point[2] for point in self.shape.extrados.bounding_box()]) * fill_percentage)
                        swt += self.shape.compute_fill_weight()
                    if rollers_ratio:
                        self.form.set_boundary_rollers(total_rx=[rollers_ratio[0]*swt0]*2, total_ry=[rollers_ratio[1]*swt0]*2)
                    elif rollers_absolute:
                        self.form.set_boundary_rollers(total_rx=[rollers_absolute[0]]*2, total_ry=[rollers_absolute[1]]*2)
                    self.apply_selfweight()
                    self.apply_envelope()
                    self.apply_reaction_bounds()
                    if fill_percentage:
                        self.apply_fill_load()
                    self.set_up_optimiser()
                    setup_time = time.time() - time0
                    self.run()
                    run_time = time.time() - time0 - setup_time
                    if plot:
                        plot_form(self.form, show_q=False, cracks=True).show()
                    exitflag = self.optimiser.exitflag  # get info if optimisation was succeded ot not
                    fopt = self.optimiser.fopt  # objective function optimum value
                    fopt_over_weight = fopt/swt  # divide by selfweight

                    if exitflag == 0:
                        if self.optimiser.data['objective'] == 'min':
                            solutions_min.append(fopt_over_weight)
                            thicknesses_min.append(thk)
                            last_min = fopt_over_weight
                            last_thk_min = thk
                        else:
                            solutions_max.append(fopt_over_weight)
                            thicknesses_max.append(thk)
                            last_max = abs(fopt_over_weight)
                            last_thk_max = thk
                        if save_forms:
                            address = save_forms + '_' + self.optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'
                        else:
                            address = os.path.join(compas_tno.get('/temp/'), 'form_' + self.optimiser.data['objective'] + '.json')
                        self.form.to_json(address)
                        print('{0}  |   True  |   {1:.1f} |   {2:.6f} |   {3}   | {4:.2f}s  |   {5:.2f}s'.format(thk, fopt, fopt_over_weight, thk_reduction, setup_time, run_time))
                        thk = round(thk - thk_reduction, 5)
                    else:
                        if not address:
                            if self.optimiser.data['objective'] == objectives[0]:
                                print('Failed in Initial Thickness and objective ({0})'.format(self.optimiser.data['objective']))
                                objectives = objectives[::-1]
                        elif self.optimiser.data['objective'] == 'max' and last_thk_min <= thk:
                            print('{0}  |   False  |   XXXX |   XXXX |   Load next [min]  | {1:.2f}s   |  {2:.2f}s'.format(thk, setup_time, run_time))
                            # print('Try loading next [min] optimisation - 1')
                            next_min = None
                            for i in range(len(thicknesses_min)):
                                if thicknesses_min[i] < thk or (thicknesses_min[i] == thk and count == 0):
                                    next_min = thicknesses_min[i]
                                    break
                            if next_min:
                                address = save_forms + '_' + 'min' + '_thk_' + str(100*next_min) + '.json'
                                self.form = FormDiagram.from_json(address)
                                thk = next_min
                                exitflag = 0
                            else:
                                print('Tried loading all [min] optimisation')
                                print('---------- End of process -----------', '\n')
                        elif thk_reduction > thk_refined:  # or (self.optimiser.data['objective'] == 'max' and last_thk_min < last_thk_max and thk_reduction > thk_refined/2):
                            self.form = FormDiagram.from_json(address)  # reload last solved
                            thk_ = thk
                            thk = round(thk + thk_reduction, 5)
                            thk_reduction = thk_reduction / 2
                            print('{0}  |   False  |   XXXX |   XXXX |   {1}    | {2:.2f}s   |  {3:.2f}s'.format(thk_, thk_reduction, setup_time, run_time))
                            thk = round(thk - thk_reduction, 5)
                            exitflag = 0
                        else:
                            print('---------- End of process -----------', '\n')
                    count += 1

        if printout:
            print('------------------ SUMMARY ------------------ ')
            print('ANALYSIS FOR THICKNESS t0:', t0)
            print('thicknesses min ({0}):'.format(len(thicknesses_min)))
            print(thicknesses_min)
            print('thicknesses max ({0}):'.format(len(thicknesses_max)))
            print(thicknesses_max)
            print('min solutions ({0}):'.format(len(solutions_min)))
            print(solutions_min)
            print('max solutions ({0}):'.format(len(solutions_max)))
            print(solutions_max)

        return [thicknesses_min, thicknesses_max],  [solutions_min, solutions_max]


    def thk_minmax_GSF(self, thk_max, thk_step=0.05, fill_percentage=None, rollers_ratio=None, rollers_absolute=None, printout=True, plot=False, save_forms=None):

        solutions_min = []  # empty lists to keep track of the solutions for min thrust
        solutions_max = []  # empty lists to keep track of the solutions for max thrust
        thicknesses_min = []
        thicknesses_max = []
        objectives = ['min', 'max']
        ro = self.shape.ro
        data_shape = self.shape.data

        # Find extreme (min thickness) solution:

        print('\n----- Starting the min thk optimisation for starting thk: {0:.4f}'.format(self.shape.data['thk']))

        time0 = time.time()
        self.optimiser.data['variables'] = ['ind', 'zb', 't']
        self.optimiser.data['objective'] = 't'
        self.apply_selfweight()
        self.apply_envelope()
        self.apply_reaction_bounds()
        self.set_up_optimiser()
        setup_time = time.time() - time0
        self.run()
        run_time = time.time() - time0 - setup_time

        exitflag = self.optimiser.exitflag

        if exitflag == 0:
            thk_min = self.form.attributes['thk']
            T = self.form.thrust()
            swt = self.shape.compute_selfweight()
            T_over_swt = T/swt
            thk_increase = 5.0*math.ceil(thk_min*100/5.0)/100 - thk_min
            print('Min THK  |   Solved  |   Thrust |   T/W   |  THK incr. |  Setup time  |   Run time')
            print('{0:.4f}  |   True  |   {1:.1f} |   {2:.6f} |   {3}   | {4:.2f}s  |   {5:.2f}s'.format(thk_min, T, T_over_swt, thk_increase, setup_time, run_time))

            if plot:
                plot_form(self.form, show_q=False, cracks=True).show()

            # STORE
            form0 = self.form.copy()
            solutions_min.append(T_over_swt)
            solutions_max.append(T_over_swt)
            thicknesses_min.append(thk_min)
            thicknesses_max.append(thk_max)

            # SAVE
            if save_forms:
                address_min = save_forms + '_' + 'min' + '_thk_' + str(100*thk_min) + '.json'
                address_max = save_forms + '_' + 'max' + '_thk_' + str(100*thk_min) + '.json'
                self.form.to_json(address_min)
                self.form.to_json(address_max)
            else:
                address = os.path.join(compas_tno.get('/temp/'), 'form_' + self.optimiser.data['objective'] + '.json')
                self.form.to_json(address)

            thk = thk_min + thk_increase
            thk_increase = thk_step
            thk0 = thk

            if thk_min == self.shape.data['thk']:
                print('Warning: Minimum THK Optimisation found optimum at start. Try rescalling the problem.')
                return
        else:
            print('Error: Minimum THK Optimisation did not find a solution: Try scalling the optimisation, or starting in a different thickness value')
            return

        # Start inverse loop for min/max:

        self.optimiser.data['variables'] = ['ind', 'zb']
        for self.optimiser.data['objective'] in objectives:

            self.form = form0
            exitflag = 0
            thk = thk0
            count = 0

            print('\n----- Starting the inverse [', self.optimiser.data['objective'], '] problem for minimum thk:', thk)
            print('THK  |   Solved  |   Opt.Val |   Opt/W   |   THK incr.  |   Setup time  |   Run time')

            while exitflag == 0 and thk <= thk_max:

                data_shape['thk'] = thk

                time0 = time.time()
                self.shape = Shape.from_library(data_shape)
                self.shape.ro = ro
                swt = self.shape.compute_selfweight()
                lumped_swt = self.form.lumped_swt()
                if fill_percentage:
                    self.shape.add_fill_with_height(max([point[2] for point in self.shape.extrados.bounding_box()]) * fill_percentage)
                    swt += self.shape.compute_fill_weight()
                if rollers_ratio:
                    self.form.set_boundary_rollers(total_rx=[rollers_ratio[0]*swt0]*2, total_ry=[rollers_ratio[1]*swt0]*2)
                elif rollers_absolute:
                    self.form.set_boundary_rollers(total_rx=[rollers_absolute[0]]*2, total_ry=[rollers_absolute[1]]*2)
                self.apply_selfweight()
                self.form.scale_form(swt/lumped_swt)
                self.apply_envelope()
                self.apply_reaction_bounds()
                if fill_percentage:
                    self.apply_fill_load()
                self.set_up_optimiser()
                setup_time = time.time() - time0
                self.run()
                run_time = time.time() - time0 - setup_time
                if plot:
                    plot_form(self.form, show_q=False, cracks=True).show()
                exitflag = self.optimiser.exitflag  # get info if optimisation was succeded ot not
                fopt = self.optimiser.fopt  # objective function optimum value
                fopt_over_weight = fopt/swt  # divide by selfweight

                if exitflag == 0:
                    if self.optimiser.data['objective'] == 'min':
                        solutions_min.append(fopt_over_weight)
                        thicknesses_min.append(thk)
                        last_min = fopt_over_weight
                        last_thk_min = thk
                    else:
                        solutions_max.append(fopt_over_weight)
                        thicknesses_max.append(thk)
                        last_max = abs(fopt_over_weight)
                        last_thk_max = thk
                    if save_forms:
                        address = save_forms + '_' + self.optimiser.data['objective'] + '_thk_' + str(100*thk) + '.json'
                    else:
                        address = os.path.join(compas_tno.get('/temp/'), 'form_' + self.optimiser.data['objective'] + '.json')
                    self.form.to_json(address)
                    print('{0}  |   True  |   {1:.1f} |   {2:.6f} |   {3}   | {4:.2f}s  |   {5:.2f}s'.format(thk, fopt, fopt_over_weight, thk_increase, setup_time, run_time))
                    thk = thk + thk_increase
                else:
                    # here put some conditions to try to make it work ....

                    print('---------- End of process -----------', '\n')
                count += 1

        if printout:
            print('------------------ SUMMARY ------------------ ')
            print('ANALYSIS FOUND MIN THK:', thk_min)
            print('thicknesses min ({0}):'.format(len(thicknesses_min)))
            print(thicknesses_min)
            print('thicknesses max ({0}):'.format(len(thicknesses_max)))
            print(thicknesses_max)
            print('min solutions ({0}):'.format(len(solutions_min)))
            print(solutions_min)
            print('max solutions ({0}):'.format(len(solutions_max)))
            print(solutions_max)

        return [thicknesses_min, thicknesses_max],  [solutions_min, solutions_max]


    def limit_analysis_load_mult(self, load0, load_increase, load_direction, percentual=True):
        """"" W.I.P. """""

        return
