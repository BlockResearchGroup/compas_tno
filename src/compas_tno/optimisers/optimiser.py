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
    M : problem
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
            'library': 'Scipy',
            'solver': 'slsqp',
            'objective': 'min',
            'constraints': ['funicular', 'envelope'],
            'features': ['fixed'],
            'variables': ['q', 'zb'],
            'find_inds': False,
            'printout': False,
            'plot': False,
            'starting_point': 'current',
            'solver-convex': 'CVX',
            'support_displacement': None,
            'gradient': True,
            'jacobian': True,
            'solver_options': {},
            'max_iter': 500,
            'qmin': -1e+4,
            'qmax': 1e-8,
        }
        self.x0 = None
        self.xopt = None
        self.fopt = None
        self.M = None
        self.message = None
        self.time = None
        self.niter = None
        self.exitflag = None
        # self.log = None
        self.callback = None

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
    def create_minthk_optimiser(cls, printout=False, plot=False , max_iter=500, starting_point='loadpath'):
        """Create a minimum thickness analysus from the elements of the problem (form and shape)

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
        starting_point : str, optional
            Which starting point use, by default 'loadpath'

        Returns
        -------
        analysis: Analysiss
            The Anallysis object

        """

        optimiser = cls()
        optimiser.settings['library'] = 'IPOPT'
        optimiser.settings['solver'] = 'IPOPT'
        optimiser.settings['constraints'] = ['funicular', 'envelope']
        optimiser.settings['variables'] = ['q', 'zb', 't']
        optimiser.settings['features'] = ['fixed']
        optimiser.settings['objective'] = 't'
        optimiser.settings['plot'] = plot
        optimiser.settings['find_inds'] = False
        optimiser.settings['max_iter'] = max_iter
        optimiser.settings['gradient'] = True
        optimiser.settings['jacobian'] = True
        optimiser.settings['printout'] = printout
        optimiser.settings['starting_point'] = starting_point

        return optimiser

    def __str__(self):
        tpl = "<Optimiser with settings {}>".format(self.settings)
        return tpl
