from compas.datastructures import Datastructure

__all__ = ['Optimiser']


class Optimiser(Datastructure):

    """The ``Optimiser`` sets the parameters of the optimisation.

    Notes
    -----
    An ``Optimiser`` keeps track in the following information in the dictionary stored in optimiser.settings:

        *  'library'           : ['Scipy','MATLAB','MMA','IPOPT','Scipy', ... others to come]
        *  'solver'            : ['slsqp','SDPT3','MMA', ... others to come],
        *  'objective'         : ['min','max','loadpath','target'],
        *  'constraints'       : ['funicular', 'envelope', 'reac_bounds', 'partial_reactions', 'cracks'],
        *  'variables'         : ['ind', 'zb', 'all-qs'],
        *  'use_indset'        : True,
        *  'solver_options'    : {},
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
            'variables': ['q', 'zb'],
            'use_indset': True,
            'solver_options': {},
            'qmin': -1e-6,
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
