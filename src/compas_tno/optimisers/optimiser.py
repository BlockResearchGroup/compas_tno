from compas.datastructures import Datastructure

__all__ = ['Optimiser']


class Optimiser(Datastructure):

    """The ``Optimiser`` sets the parameters of the optimisation.

    Notes
    -----
    An ``Optimiser`` keeps track in the following information in the dictionary stored in optimiser.settings:

        *  'library'           : ['Scipy','MATLAB','MMA','IPOPT','Scipy', ... others to come]
        *  'solver'            : ['slsqp','SDPT3','MMA', ... others to come],
        *  'objective'         : ['min','max','t','loadpath','target'],
        *  'constraints'       : ['funicular', 'envelope', 'reac_bounds', 'partial_reactions', 'cracks'],
        *  'variables'         : ['ind', 'zb', 'all-qs'],
        *  'find_inds'         : True,
        *  'exitflag'          : 1,
        *  'fopt'              : None,
        *  'xopt'              : None,

    """

    # __module__ = 'compas_tno.optimisers'

    def __init__(self):
        self.settings = {
            'library': 'Scipy',
            'solver': 'slsqp',
            'objective': 'min',
            'constraints': ['funicular', 'envelope'],
            'features': ['fixed'],
            'variables': ['q', 'zb'],
            'find_inds': True,
            'printout': False,
            'starting_point': 'current',
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
        self.mesage = None
        self.time = None
        self.niter = None
        self.exitflag = None
        self.log = None

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
            'log': self.log,
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
        self.log = data.get('log', None)
