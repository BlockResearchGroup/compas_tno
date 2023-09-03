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
    """The ``Analysis`` class unites the :class:`~compas_tno.diagrams.FormDiagram`, the :class:`~compas_tno.shapes.Shape`
    and the :class:`~compas_tno.optimisers.Optimiser` to compute the optimisation.

    Examples
    --------
    >>> from compas_tno.analysis import Analysis
    >>> from compas_tno.shapes import Shape
    >>> from compas_tno.diagrams import FormDiagram
    >>> arch = Shape.create_arch()
    >>> form = FormDiagram.create_arch()
    >>> analysis = Analysis.create_minthrust_analysis(form, arch)

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
        """Create a analysis from the elements of the problem (form, shape and optimiser)

        Parameters
        ----------
        shape : :class:`~compas_tno.shapes.Shape`
            The vaulted structure constraining the problem
        form : :class:`~compas_tno.diagrams.FormDiagram`
            The form diagram to analyse the structure
        optimiser : :class:`~compas_tno.optimisers.Optimiser`
            The optimiser with information about the problem

        Returns
        -------
        analysis: Analysis
            The analysis object

        """

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

    @classmethod
    def create_minthk_analysis(cls, form, shape, printout=False, plot=False, max_iter=500, starting_point='loadpath', solver='SLSQP', derivatives=True):
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problem
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        optimiser = Optimiser.create_minthk_optimiser(printout=printout,
                                                      plot=plot,
                                                      max_iter=max_iter,
                                                      starting_point=starting_point,
                                                      solver=solver,
                                                      derivatives=derivatives)

        if printout:
            print('-'*20)
            print('Minimum thickness analysis created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_bestfit_analysis(cls, form, shape, printout=False, plot=False, max_iter=500, starting_point='loadpath', solver='SLSQP', derivatives=True):
        """Create a bestfit analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problem
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        derivatives : bool, optional
            Whether or not derivatives should be considered, by default False

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        optimiser = Optimiser.create_bestfit_optimiser(printout=printout,
                                                       plot=plot,
                                                       max_iter=max_iter,
                                                       starting_point=starting_point,
                                                       solver=solver,
                                                       derivatives=derivatives)

        if printout:
            print('-'*20)
            print('Minimum thickness analysis created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_minthrust_analysis(cls, form, shape, printout=False, plot=False, max_iter=500, starting_point='loadpath', solver='SLSQP'):
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problemf
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        optimiser = Optimiser.create_minthrust_optimiser(printout=printout,
                                                         plot=plot,
                                                         max_iter=max_iter,
                                                         starting_point=starting_point,
                                                         solver=solver)

        if printout:
            print('-'*20)
            print('Minimum thrust analysis created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_maxthrust_analysis(cls, form, shape, printout=False, plot=False, max_iter=500, starting_point='loadpath', solver='SLSQP'):
        """Create a maximum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problemf
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        optimiser = Optimiser.create_maxthrust_optimiser(printout=printout,
                                                         plot=plot,
                                                         max_iter=max_iter,
                                                         starting_point=starting_point,
                                                         solver=solver)

        if printout:
            print('-'*20)
            print('Maximium thrust analysis created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_max_load_analysis(cls, form, shape, printout=False, plot=False, horizontal=False, max_iter=500, starting_point='loadpath',
                                 solver='IPOPT', derivatives=True, load_direction=None, max_lambd=1.0):
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problem
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        plot : bool, optional
            Whether or load applied is horizontal, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        if horizontal:
            optimiser = Optimiser.create_max_horload_optimiser(printout=printout,
                                                               plot=plot,
                                                               max_iter=max_iter,
                                                               starting_point=starting_point,
                                                               solver=solver,
                                                               derivatives=derivatives,
                                                               load_direction=load_direction,
                                                               max_lambd=max_lambd)

            if printout:
                print('-'*20)
                print('Max horizontal load analysis created')
                # print(optimiser)

        else:
            optimiser = Optimiser.create_max_vertload_optimiser(printout=printout,
                                                                plot=plot,
                                                                max_iter=max_iter,
                                                                starting_point=starting_point,
                                                                solver=solver,
                                                                derivatives=derivatives,
                                                                load_direction=load_direction,
                                                                max_lambd=max_lambd)

            if printout:
                print('-'*20)
                print('Max vertical load analysis created')
                # print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_compl_energy_analysis(cls, form, shape, printout=False, solver='IPOPT', plot=False, max_iter=500, starting_point='loadpath',
                                     support_displacement=None, Emethod='simplified'):
        """Create a complementary energy analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problem
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        support_displacement : array [nb x 3], optional
            Vector with the displacement applied to the supports, by default None
        Emethod : str, optional
            Whether or not internal deformations should be considered, by default 'simplified' which considers the material rigid

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        optimiser = Optimiser.create_compl_energy_optimiser(printout=printout,
                                                            plot=plot,
                                                            max_iter=max_iter,
                                                            starting_point=starting_point,
                                                            support_displacement=support_displacement,
                                                            Emethod=Emethod,
                                                            solver=solver)

        if printout:
            print('-'*20)
            print('Complementary energy created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_quad_compl_energy_analysis(cls, form, shape, printout=False, solver='IPOPT', plot=False, max_iter=500,
                                          starting_point='loadpath', support_displacement=None, Emethod='simplified'):
        """Create a complementary energy analysis including a quadratic term from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`
            The shape constraining the problemf
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        analysis = cls().from_form_and_shape(form, shape)

        optimiser = Optimiser.create_compl_energy_optimiser(printout=printout,
                                                            plot=plot,
                                                            max_iter=max_iter,
                                                            starting_point=starting_point,
                                                            support_displacement=support_displacement,
                                                            Emethod=Emethod,
                                                            solver=solver)

        if printout:
            print('-'*20)
            print('Complementary energy analysis created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    @classmethod
    def create_lp_analysis(cls, form, shape=None, solver='CVXPY', printout=False, plot=False, max_iter=500):
        """Create a minimum thickness analysis from the elements of the problem (form and shape)

        Parameters
        ----------
        form : :class:`~compas_tno.diagrams.FormDiagram`
            _description_
        shape : :class:`~compas_tno.shapes.Shape`, optional
            The shape constraining the problem, by default None
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysis
            The Analysis object

        """

        optimiser = Optimiser.create_lp_optimiser(solver=solver,
                                                  printout=printout,
                                                  plot=plot,
                                                  max_iter=max_iter)

        if shape:
            analysis = cls().from_elements(shape, form, optimiser)
        else:
            analysis = cls().from_form_and_optimiser(form, optimiser)

        if printout:
            print('-'*20)
            print('Load path Analysis created')
            print(optimiser)

        analysis.optimiser = optimiser

        return analysis

    def is_convex(self):
        """Check if the analysis problem is convex."""

        if not self.optimiser:
            raise ValueError('Define the Optimiser for the problem')

        objective = self.optimiser.settings['objective']
        features = self.optimiser.settings['features']

        if objective == 'loadpath' and 'fixed' in features:
            return True
        else:
            return False

    def set_optimiser_options(self, **kwargs):
        """Set the additional options of the optimisation.
        """

        if kwargs:
            self.optimiser.set_additional_options(**kwargs)

    def clear_previous_results(self):
        """Clear Previous results stored in the Analysis object. Necessary to perform sequential optimisation with the same analysis object
        """

        self.optimiser.clear_optimiser()

        return

    def apply_selfweight(self, normalize_loads=True):
        """Invoke method to apply selfweight to the nodes of the form diagram based on the shape"""
        apply_selfweight_from_shape(self.form, self.shape, normalize=normalize_loads)

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

        if self.is_convex():
            set_up_convex_optimisation(self)
        else:
            self = set_up_general_optimisation(self)

        return

    def run(self):
        """With the data from the elements of the problem compute the matrices for the optimisation"""

        solver = self.optimiser.settings.get('solver', 'SLSQP')

        if not isinstance(solver, str):
            raise ValueError('Please provide the name of the solver')

        if self.is_convex():
            if solver == 'MATLAB':
                run_optimisation_MATLAB(self)
            elif solver == 'CVXPY':
                run_optimisation_CVXPY(self)
            else:
                raise NotImplementedError('Only <CVXPY> and <MATLAB> are suitable for this optimisation')
        else:
            if solver.split('-') == 'pyOpt':
                self = run_optimisation_MMA(self)  # change to PyOpt
            elif solver == 'MMA':
                self = run_optimisation_MMA(self)
            elif solver == 'IPOPT':
                self = run_optimisation_ipopt(self)
            else:
                self = run_optimisation_scipy(self)

        return

    def __str__(self):
        tpl = "<Analysis with parameters: {} >".format(self.data)
        return tpl
