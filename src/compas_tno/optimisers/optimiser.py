


__all__ = ['Optimiser']

class Optimiser(object):

    """The ``Optimiser`` sets the parameters of the optimisation.

    Notes
    -----
    An ``Optimiser`` keeps track in the following information in the dictionary stored in Optimiser.data:

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

    __module__ = 'compas_tna.optimiser'


    def __init__(self):
        self.data = {
        'library'           : 'Scipy',
        'solver'            : 'slsqp',
        'objective'         : 'min',
        'constraints'       : ['funicular', 'envelope', 'reac_bounds', 'partial_reactions', 'cracks'],
        'variables'         : ['ind', 'zb', 'all-qs'],
        'use_indset'        : True,
        'solver_options'    : {},
        'qmin'              : -1e-6,
        }
        self.fobj = None
        self.fconstr = None
        self.x0 = None
        self.xopt = None
        self.fopt = None

    # This class must separate the functions objective functions and etc...
