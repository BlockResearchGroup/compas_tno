from compas.datastructures import Datastructure

__all__ = ['Optimiser']


class Optimiser(Datastructure):

    """The ``Optimiser`` sets the parameters of the optimisation.

    Parameters
    ----------
    x0 : array
        The starting point of the optimisation.
    xopt : array
        The value of the variables that solves thhe optimisation problem.
    fopt : float
        The value of the objective function at the end of the optimisation process.
    message : str
        The message provided by the solver.
    M : :class:`~compas_tno.problems.Problem`
        The problem with the relevat matrices.
    time : float
        The time consumption for solving.
    niter : int
        The number of iterations.
    exitflag : int
        Whether or not the optimisation problem was solved successfully.
        Value ``0``represents a solved problem.
    log : str
        Log of the optimisation, available for some solves.

    Settings
    --------
    The main settings and default values are shown below:

    *  'library'           : ['Scipy', 'MATLAB', 'MMA', 'IPOPT', 'Scipy', ...]
    *  'solver'            : ['SLSQP', 'IPOPT', 'MMA', ...],
    *  'objective'         : ['min', 'max', 't', 'loadpath', 'bestfit', ...],
    *  'constraints'       : ['funicular', 'envelope', 'reac_bounds', ...],
    *  'variables'         : ['q', 'zb', 'xyb', 't', ...],
    *  'features'          : ['fixed', 'sym', ...],
    *  'find_inds'         : True,
    *  'exitflag'          : 1,
    *  'printout'          : True,
    *  'starting_point'    : 'current',
    *  'support_displacement': None,
    *  'gradient'          : True,
    *  'jacobian'          : True,
    *  'solver_options'    : {},
    *  'max_iter'          : 500,
    *  'qmin'              : -1e+4,
    *  'qmax'              : 1e-8,


    """

    def __init__(self):
        self.settings = {
            'solver': 'SLSQP',
            'objective': 'min',
            'constraints': ['funicular', 'envelope'],
            'features': ['fixed'],
            'variables': ['q', 'zb'],
            'find_inds': False,
            'printout': False,
            'plot': False,
            'starting_point': 'current',
            'solver_convex': 'CVXPY',
            'support_displacement': None,
            'gradient': True,
            'jacobian': True,
            'solver_options': {},
            'max_iter': 500,
            'qmin': -1e+4,
            'qmax': 1e-8,
        }
        self.clear_optimiser()

    @property
    def data(self):
        """dict : A data dict representing the shape data structure for serialization.
        """
        data = {
            'settings': self.settings,
            'x0': self.x0,
            'xopt': self.xopt,
            'fopt': self.fopt,
            'message': self.message,
            'niter': self.niter,
            'exitflag': self.exitflag,
            # 'log': self.log,
        }
        return data

    @data.setter
    def data(self, data):
        if 'data' in data:
            data = data['data']
        self.settings = data.get('settings') or {}

        self.x0 = data.get('x0', None)
        self.xopt = data.get('xopt', None)
        self.fopt = data.get('fopt', None)
        self.message = data.get('message', None)
        self.niter = data.get('niter', None)
        self.exitflag = data.get('exitflag', None)
        # self.log = data.get('log', None)

    @classmethod
    def create_minthk_optimiser(cls, solver='SLSQP', max_iter=500, printout=False, plot=False, starting_point='loadpath', derivatives=True):
        """Create a minimum thickness optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        derivatives : bool, optional
            Whether or not derivatives are computed by hand, by default True

        Returns
        -------
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb', 't'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('t')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=derivatives, jacobian=derivatives)
        optimiser.set_starting_point(starting_point=starting_point)

        return optimiser

    @classmethod
    def create_minthrust_optimiser(cls, solver='SLSQP', max_iter=500, printout=False, plot=False, starting_point='loadpath'):
        """Create a minimum thickness optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
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
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('min')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=True, jacobian=True)
        optimiser.set_starting_point(starting_point=starting_point)

        return optimiser

    @classmethod
    def create_maxthrust_optimiser(cls, solver='SLSQP', max_iter=500, printout=False, plot=False, starting_point='loadpath'):
        """Create a maximum thickness optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
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
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls().create_minthrust_optimiser(solver=solver,
                                                     max_iter=max_iter,
                                                     printout=printout,
                                                     plot=plot,
                                                     starting_point=starting_point)
        optimiser.set_objective('max')

        return optimiser

    @classmethod
    def create_max_horload_optimiser(cls, solver='IPOPT', max_iter=500, printout=False, plot=False, starting_point='loadpath', max_lambd=1.0,
                                     load_direction=None, derivatives=True):
        """Create a minimum thickness optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        max_lambd : float, optional
            Maximum load multiplier, by default 1.0
        load_direction : array [2n x 1], optional
            Direction of the horizontal load, by default None
        derivatives : bool, optional
            Whether or not derivatives are computed by hand, by default True

        Returns
        -------
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb', 'lambdh'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('max_load')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=derivatives, jacobian=derivatives)
        optimiser.set_starting_point(starting_point=starting_point)
        optimiser.set_additional_options(max_lambd=max_lambd, load_direction=load_direction)

        return optimiser

    @classmethod
    def create_bestfit_optimiser(cls, solver='SLSQP', max_iter=500, printout=False, plot=False, starting_point='loadpath', derivatives=True):
        """Create a bestfit optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
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
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('bestfit')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=derivatives, jacobian=derivatives)
        optimiser.set_starting_point(starting_point=starting_point)

        return optimiser

    @classmethod
    def create_max_vertload_optimiser(cls, solver='IPOPT', max_iter=500, printout=False, plot=False, starting_point='loadpath', max_lambd=1.0,
                                      load_direction=None, derivatives=True):
        """Create a minimum thickness optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500
        starting_point : str, optional
            Which starting point use, by default 'loadpath'
        max_lambd : float, optional
            Maximum load multiplier, by default 1.0
        load_direction : array [n x 1], optional
            Direction of the vertical external load, by default None
        derivatives : bool, optional
            Whether or not derivatives are computed by hand, by default True

        Returns
        -------
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb', 'lambdv'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('max_load')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=derivatives, jacobian=derivatives)
        optimiser.set_starting_point(starting_point=starting_point)
        optimiser.set_additional_options(max_lambd=max_lambd, load_direction=load_direction)

        return optimiser

    @classmethod
    def create_lp_optimiser(cls, solver='MATLAB', printout=False, plot=False, max_iter=500):
        """Create a loadpath optimisation optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        form : FormDiagram
            _description_
        shape : Shape
            The shape cconstraining the problem
        printout : bool, optional
            Whether or not prints appear in the creen, by default False
        plot : bool, optional
            Whether or not plots showing intermediate states appear, by default False
        max_iter : int, optional
            Maximum number of itetations, by default 500

        Returns
        -------
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        if solver not in ['MATLAB', 'CVXPY']:
            raise ValueError('For loadpath optimisation only MATLAB or CVXPY are possible solvers. See solvers page.')

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular'])
        optimiser.set_variables(['q'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('loadpath')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=True, jacobian=True)

        return optimiser

    @classmethod
    def create_compl_energy_optimiser(cls, solver='SLSQP', max_iter=500, printout=False, plot=False, starting_point='loadpath',
                                      support_displacement=None, Emethod='simplified'):
        """Create a linear complementary energy optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
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
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('Ecomp-linear')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=True, jacobian=True)
        optimiser.set_starting_point(starting_point=starting_point)
        optimiser.set_additional_options(support_displacement=support_displacement, Ecomp_method=Emethod)

        return optimiser

    @classmethod
    def create_quad_compl_energy_optimiser(cls, solver='SLSQP', max_iter=500, printout=False, plot=False,
                                           starting_point='loadpath', support_displacement=None, Emethod='simplified'):
        """Create a quadratic complementary energy optimiser to be sent with instructions to the Analysis.

        Parameters
        ----------
        solver : str, optional
            Which solver to use, by default 'SLSQP'. See Solvers page for more information.
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
        :class:`~compas_tno.optimisers.Optimiser`
            The Optimiser object

        """

        optimiser = cls()
        optimiser.set_solver(solver)
        optimiser.set_constraints(['funicular', 'envelope'])
        optimiser.set_variables(['q', 'zb'])
        optimiser.set_features(['fixed'])
        optimiser.set_objective('Ecomp-nonlinear')
        optimiser.set_display_options(plot=plot, printout=printout)
        optimiser.set_max_iterations(max_iter=max_iter)
        optimiser.set_gradient_options(gradient=True, jacobian=True)
        optimiser.set_starting_point(starting_point=starting_point)
        optimiser.set_additional_options(support_displacement=support_displacement, Ecomp_method=Emethod)

        return optimiser

    def set_objective(self, objective='min'):
        """Set the objective function for the problem.

        Parameters
        ----------
        objective : str, optional
            Name of the objective function to call, by default 'min'
            For a list of possible objective functions see [link].

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(objective, str):
            raise ValueError('Please provide the name of the objective')

        self.settings['objective'] = objective

    def set_constraints(self, constraints=[]):
        """Set the constraints for the problem.

        Parameters
        ----------
        constraints : [str], optional
            List wihth the names of the constraints to activate, by default []
            For a list of possible constraints functions see [link].

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(constraints, list):
            raise ValueError('Please provide a list with the name of constraints')

        self.settings['constraints'] = constraints

    def set_variables(self, variables=[]):
        """Set the variables for the problem.

        Parameters
        ----------
        variables : [str], optional
            List wihth the names of the variables to activate, by default []
            For a list of possible variables functions see [link].

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(variables, list):
            raise ValueError('Please provide a list with the name of variables')

        self.settings['variables'] = variables

    def set_features(self, features=[]):
        """Set the features for the problem.

        Parameters
        ----------
        features : [str], optional
            List with the names of the features to activate, by default []
            For a list of possible features functions see [link].

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(features, list):
            raise ValueError('Please provide a list with the name of features')

        self.settings['features'] = features

    def set_solver(self, solver='SLSQP'):
        """Set the features for the problem.

        Parameters
        ----------
        solver : str, optional
            Set the solver to use.
            For a list of possible solvers see [link].

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(solver, str):
            raise ValueError('Please provide the name of the solver')

        self.settings['solver'] = solver

    def set_max_iterations(self, max_iter=500):
        """Set the features for the problem.

        Parameters
        ----------
        max_iter : int, optional
            Number of iterations, by default 500.

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(max_iter, int):
            raise ValueError('Please provide a int as max iteration')

        self.settings['max_iter'] = max_iter

    def set_starting_point(self, starting_point='loadpath'):
        """Set the starating point of the problem.

        Parameters
        ----------
        starting_point : str, optional
            Set the starting point to use.
            For a list of possible starting points see [link].

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(starting_point, str):
            raise ValueError('Please provide the name of the starting point')

        self.settings['starting_point'] = starting_point

    def set_display_options(self, plot=False, printout=True):
        """Set the display options of the optimisation

        Parameters
        ----------
        plot : bool, optional
            If plots with partial states should be evokes, by default False
        printout : bool, optional
            If prints with iterations should appear in the screen, by default True

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(plot, bool):
            raise ValueError('Please provide a bool')
        if not isinstance(printout, bool):
            raise ValueError('Please provide a bool')

        self.settings['plot'] = plot
        self.settings['printout'] = printout

    def set_gradient_options(self, gradient=True, jacobian=True):
        """Set the display options of the optimisation

        Parameters
        ----------
        gradient : bool, optional
            If gradient should be provided analytically, by default True
        jacobian : bool, optional
            If jacobian should be provided analytically, by default True

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(gradient, bool):
            raise ValueError('Please provide a bool')
        if not isinstance(jacobian, bool):
            raise ValueError('Please provide a bool')

        self.settings['gradient'] = gradient
        self.settings['jacobian'] = jacobian

    def set_axis_symmetry(self, axis_sym=[]):
        """Set the axis of symmetry.

        Parameters
        ----------
        axis_sym : list, optional
            List with the lines of symmetry

        Returns
        -------
        None
            Optimiser is modified in place
        """

        if not isinstance(axis_sym, list):
            raise ValueError('Please provide a list with the name of features')

        self.settings['axis_sym'] = axis_sym

    def set_additional_options(self, **kwargs):
        """Set the additional options of the optimisation.
        """

        if kwargs:
            for key in kwargs:
                self.settings[key] = kwargs[key]

    def clear_optimiser(self):
        """Initiate an empty Optimiser objects, or clear Previous results stored in the Optimiser object.
        Necessary to perform sequential optimisation with the same analysis object
        """

        self.x0 = None
        self.xopt = None
        self.fopt = None
        self.M = None
        self.message = None
        self.time = None
        self.niter = None
        self.exitflag = None
        self.callback = None

        return

    def __str__(self):
        tpl = "-"*20 + '\nOptimiser with settings:\n'
        for key in self.settings:
            tpl = tpl + key + ' : ' + str(self.settings[key]) + '\n'
        tpl = tpl + "-"*20 + '\n'
        # tpl = "<Optimiser with settings {}>".format(self.settings)
        return tpl
